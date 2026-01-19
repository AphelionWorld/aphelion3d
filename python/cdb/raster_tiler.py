from __future__ import annotations
from pathlib import Path
from typing import Callable
import numpy as np
import rasterio
from rasterio.warp import transform_bounds, reproject, Resampling
from rasterio.crs import CRS as RioCRS
from rasterio.transform import from_bounds as rio_from_bounds
from .tiling import tile_bounds, MAX_LOD, get_zone_width, geocell_for_lonlat
from .paths import CDBPaths
from .datasets import Dataset, ELEVATION


def tile_raster_to_cdb(
    input_path: Path,
    cdb: CDBPaths,
    min_lod: int = 0,
    max_lod: int = MAX_LOD,
    tile_size: int = 256,
    dataset: Dataset = ELEVATION,
    progress_callback: Callable[[str], None] | None = None,
) -> list[Path]:
    """
    Tile an input elevation raster into CDB geocell tiles at multiple LODs.
    
    CDB 1.2 structure:
    Tiles/{lat_cell}/{lon_cell}/{code}_{name}/L{lod:02d}/U{u}/{lat_cell}{lon_cell}_D{code:03d}_S{S}_T{T}_L{lod:02d}_U{u}_R{r}.tif
    
    Example:
    Tiles/N32/W118/001_Elevation/L00/U0/N32W118_D001_S001_T001_L00_U0_R0.tif
    
    Args:
        input_path: Path to input raster (GeoTIFF, etc.)
        cdb: CDB paths configuration
        min_lod: Minimum LOD to generate (default 0)
        max_lod: Maximum LOD to generate (default MAX_LOD=5)
        tile_size: Output tile size in pixels (default 256)
        dataset: Dataset to write to (default ELEVATION)
        progress_callback: Optional callback for progress messages
    
    Returns:
        List of paths to generated tiles
    """
    def log(msg: str):
        if progress_callback:
            progress_callback(msg)
    
    written = []
    
    # Ensure input path is absolute
    input_path = Path(input_path).resolve()
    cdb_root = Path(cdb.root).resolve()
    
    log(f"CDB Root: {cdb_root}")
    log(f"Input file: {input_path}")
    
    with rasterio.open(input_path) as src:
        # Get bounds in WGS84
        src_crs = src.crs
        wgs84 = RioCRS.from_epsg(4326)
        
        if src_crs != wgs84:
            log(f"Transforming from {src_crs} to WGS84")
            west, south, east, north = transform_bounds(src_crs, wgs84, *src.bounds)
        else:
            west, south, east, north = src.bounds
        
        log(f"Input bounds (WGS84): {west:.6f}, {south:.6f}, {east:.6f}, {north:.6f}")
        
        # Determine geocells covered using CDB zone-aware logic
        geocells_set = set()
        # Sample points across the bounding box to find all geocells
        lat_samples = np.linspace(south, north, max(int((north - south) * 10), 10))
        lon_samples = np.linspace(west, east, max(int((east - west) * 10), 10))
        
        for lat_sample in lat_samples:
            for lon_sample in lon_samples:
                lon_cell, lat_cell, zone_width = geocell_for_lonlat(lon_sample, lat_sample)
                geocells_set.add((lon_cell, lat_cell, zone_width))
        
        geocells = list(geocells_set)
        
        log(f"Processing {len(geocells)} geocell(s) at LODs {min_lod}-{max_lod}")
        
        # Count total tiles to process
        total_tiles = 0
        for lod in range(min_lod, max_lod + 1):
            n_tiles = 1 << lod
            for lon_cell, lat_cell, zone_width in geocells:
                for u in range(n_tiles):
                    for v in range(n_tiles):
                        tile_w, tile_s, tile_e, tile_n = tile_bounds(
                            lon_cell, lat_cell, lod, u, v, zone_width
                        )
                        if not (tile_e < west or tile_w > east or
                                tile_n < south or tile_s > north):
                            total_tiles += 1
        
        log(f"Total tiles to generate: {total_tiles}")
        
        processed = 0
        last_pct = -1
        
        # Process each LOD
        for lod in range(min_lod, max_lod + 1):
            n_tiles = 1 << lod  # tiles per side in a geocell
            
            for lon_cell, lat_cell, zone_width in geocells:
                for u in range(n_tiles):
                    for v in range(n_tiles):
                        tile_w, tile_s, tile_e, tile_n = tile_bounds(
                            lon_cell, lat_cell, lod, u, v, zone_width
                        )
                        
                        # Check if tile overlaps input raster bounds
                        if (tile_e < west or tile_w > east or
                            tile_n < south or tile_s > north):
                            continue
                        
                        try:
                            # Format geocell coordinates with N/S and E/W
                            if lat_cell >= 0:
                                lat_str = f"N{lat_cell:02d}"
                            else:
                                lat_str = f"S{abs(lat_cell):02d}"
                            
                            if lon_cell >= 0:
                                lon_str = f"E{lon_cell:03d}"
                            else:
                                lon_str = f"W{abs(lon_cell):03d}"
                            
                            # Dataset folder format: {code}_{name}
                            dataset_folder = f"{dataset.code:03d}_{dataset.name}"
                            
                            # CDB 1.2 compliant path structure
                            # Tiles/{lat_cell}/{lon_cell}/{code}_{name}/L{lod:02d}/U{v}/
                            # U folder groups tiles by row (v = latitude index)
                            tile_dir = (
                                cdb_root / "Tiles" / lat_str / lon_str / dataset_folder / 
                                f"L{lod:02d}" / f"U{v}"
                            )
                            
                            # Create directory
                            tile_dir.mkdir(parents=True, exist_ok=True)
                            
                            # Filename format: {lat_cell}{lon_cell}_D{code:03d}_S{S}_T{T}_L{lod:02d}_U{v}_R{u}.tif
                            # U is the row (v = latitude tile index within geocell)
                            # R is the column (u = longitude tile index within geocell)
                            # S and T are selector indices (typically 001 for single series/tile)
                            filename = f"{lat_str}{lon_str}_D{dataset.code:03d}_S001_T001_L{lod:02d}_U{v}_R{u}.{dataset.extension}"
                            output_path = tile_dir / filename
                            
                            log(f"Writing: {output_path.relative_to(cdb_root)}")
                            
                            # Create optimized sampler with windowed read
                            sampler = create_windowed_sampler(
                                src, tile_w, tile_s, tile_e, tile_n, tile_size, src_crs, wgs84
                            )
                            
                            # Write directly to the output path
                            write_tile_direct(
                                output_path=output_path,
                                west=tile_w,
                                south=tile_s,
                                east=tile_e,
                                north=tile_n,
                                size=tile_size,
                                sampler=sampler,
                            )
                            
                            written.append(output_path)
                            
                        except Exception as e:
                            log(f"  Warning: Failed to write tile LOD{lod} {lat_str}{lon_str} U{u} V{v}: {e}")
                            import traceback
                            traceback.print_exc()
                        
                        processed += 1
                        pct = int((processed / total_tiles) * 100)
                        if pct != last_pct:
                            log(f"Progress: {pct}% ({processed}/{total_tiles} tiles)")
                            last_pct = pct
        
        log(f"Complete! Generated {len(written)} tiles.")
    
    return written


def write_tile_direct(
    output_path: Path,
    west: float,
    south: float,
    east: float,
    north: float,
    size: int,
    sampler: Callable[[np.ndarray, np.ndarray], np.ndarray],
    band_count: int = 1,
    dtype: str = 'float32',
) -> None:
    """
    Write a tile directly to the specified path.
    
    Args:
        output_path: Full path including filename
        west, south, east, north: Tile bounds in WGS84
        size: Tile size in pixels
        sampler: Function to get raster data (returns shape (bands, height, width) or (height, width))
        band_count: Number of bands to write
        dtype: Data type for output
    """
    from rasterio.transform import from_bounds
    from rasterio.crs import CRS
    
    # Ensure output path is absolute
    output_path = Path(output_path).resolve()
    
    # Create coordinate grids
    lons = np.linspace(west, east, size)
    lats = np.linspace(south, north, size)
    xx, yy = np.meshgrid(lons, lats[::-1])  # Flip lats for north-up
    
    # Sample raster data
    raster_data = sampler(xx, yy)
    
    # Handle both single-band and multi-band data
    if raster_data.ndim == 2:
        # Single band: shape is (height, width)
        raster_data = raster_data[np.newaxis, :, :]  # Add band dimension
    
    # Create geotransform
    transform = from_bounds(west, south, east, north, size, size)
    
    # Write GeoTIFF
    with rasterio.open(
        output_path,
        'w',
        driver='GTiff',
        height=size,
        width=size,
        count=raster_data.shape[0],
        dtype=raster_data.dtype,
        crs=CRS.from_epsg(4326),
        transform=transform,
        compress='deflate',
    ) as dst:
        # Write all bands
        for band_idx in range(raster_data.shape[0]):
            dst.write(raster_data[band_idx], band_idx + 1)


def create_windowed_sampler(
    src: rasterio.DatasetReader,
    west: float,
    south: float,
    east: float,
    north: float,
    tile_size: int,
    src_crs: RioCRS,
    target_crs: RioCRS,
) -> Callable[[np.ndarray, np.ndarray], np.ndarray]:
    """
    Create an optimized sampler that reads a window from the source raster once
    and uses it for all queries. Much faster than pixel-by-pixel access.
    Handles both single-band and multi-band rasters.
    
    Args:
        src: Open rasterio dataset
        west, south, east, north: Tile bounds in target CRS (WGS84)
        tile_size: Output tile resolution
        src_crs: Source raster CRS
        target_crs: Target CRS (WGS84)
    
    Returns:
        Sampler function that takes (xx, yy) arrays and returns raster data
        with shape (bands, height, width) or (height, width) for single band
    """
    # Get number of bands
    band_count = src.count
    
    # Pre-allocate output array for the entire tile
    if band_count == 1:
        out_data = np.zeros((tile_size, tile_size), dtype=np.float32)
    else:
        out_data = np.zeros((band_count, tile_size, tile_size), dtype=np.float32)
    
    # Create transform for output tile in WGS84
    dst_transform = rio_from_bounds(west, south, east, north, tile_size, tile_size)
    
    # Reproject all bands
    try:
        if band_count == 1:
            # Single band
            reproject(
                source=rasterio.band(src, 1),
                destination=out_data,
                src_transform=src.transform,
                src_crs=src_crs,
                dst_transform=dst_transform,
                dst_crs=target_crs,
                resampling=Resampling.bilinear,
                src_nodata=src.nodata,
                dst_nodata=0.0,
            )
        else:
            # Multiple bands
            for band_idx in range(1, band_count + 1):
                reproject(
                    source=rasterio.band(src, band_idx),
                    destination=out_data[band_idx - 1],
                    src_transform=src.transform,
                    src_crs=src_crs,
                    dst_transform=dst_transform,
                    dst_crs=target_crs,
                    resampling=Resampling.bilinear,
                    src_nodata=src.nodata,
                    dst_nodata=0.0,
                )
    except Exception as e:
        # If reprojection fails, return zeros
        print(f"Warning: Reprojection failed: {e}")
        out_data[:] = 0.0
    
    # Return a sampler that uses the pre-computed data
    def sampler(xx: np.ndarray, yy: np.ndarray) -> np.ndarray:
        """
        Sample from the pre-loaded tile data.
        xx, yy are in WGS84 lon/lat and should match the tile bounds.
        Returns shape (bands, height, width) or (height, width) for single band.
        """
        # Map lon/lat to pixel indices in the pre-loaded tile
        # Normalize to [0, 1] within tile bounds
        u = (xx - west) / (east - west)
        v = (yy - south) / (north - south)
        
        # Convert to pixel indices
        col = (u * tile_size).astype(int)
        row = ((1.0 - v) * tile_size).astype(int)  # Flip Y (north-up to row-down)
        
        # Clamp to valid range
        col = np.clip(col, 0, tile_size - 1)
        row = np.clip(row, 0, tile_size - 1)
        
        # Lookup values
        if band_count == 1:
            result = out_data[row, col]
        else:
            result = out_data[:, row, col]
        
        return result
    
    return sampler