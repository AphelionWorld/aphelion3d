from __future__ import annotations
from pathlib import Path
from typing import Callable
import numpy as np
import geopandas as gpd
from shapely.geometry import box
from .tiling import tile_bounds, MAX_LOD, geocell_for_lonlat
from .paths import CDBPaths
from .datasets import Dataset, VECTOR


def tile_vector_to_cdb(
    input_path: Path,
    cdb: CDBPaths,
    min_lod: int = 0,
    max_lod: int = MAX_LOD,
    dataset: Dataset = VECTOR,
    output_format: str = "GPKG",
    progress_callback: Callable[[str], None] | None = None,
) -> list[Path]:
    """
    Tile an input vector dataset into CDB geocell tiles at multiple LODs.
    Vectors are clipped to tile bounds.
    
    Tiles/{lat_cell}/{lon_cell}/{code}_{name}/L{lod:02d}/U{u}/{lat_cell}{lon_cell}_D{code:03d}_S{S}_T{T}_L{lod:02d}_U{u}_R{r}.{ext}
    
    Example:
    Tiles/N32/W118/200_VectorData/L00/U0/N32W118_D200_S001_T001_L00_U0_R0.gpkg
    
    Args:
        input_path: Path to input vector file (GeoPackage, Shapefile, etc.)
        cdb: CDB paths configuration
        min_lod: Minimum LOD to generate (default 0)
        max_lod: Maximum LOD to generate (default MAX_LOD=5)
        dataset: Dataset to write to (default VECTOR)
        output_format: Output format "GPKG" or "ESRI Shapefile" (default GPKG)
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
    
    # Load vector data
    gdf = gpd.read_file(input_path)
    
    # Ensure WGS84
    if gdf.crs != "EPSG:4326":
        log(f"Transforming from {gdf.crs} to WGS84")
        gdf = gdf.to_crs("EPSG:4326")
    
    # Get bounds
    bounds = gdf.total_bounds  # [minx, miny, maxx, maxy]
    west, south, east, north = bounds
    
    log(f"Input bounds (WGS84): {west:.6f}, {south:.6f}, {east:.6f}, {north:.6f}")
    log(f"Features: {len(gdf)}")
    
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
    
    # Determine output extension based on format
    if output_format.upper() == "GPKG":
        ext = "gpkg"
        driver = "GPKG"
    elif output_format.upper() in ("SHAPEFILE", "ESRI SHAPEFILE"):
        ext = "shp"
        driver = "ESRI Shapefile"
    else:
        raise ValueError(f"Unsupported output format: {output_format}")
    
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
        n_tiles = 1 << lod
        
        for lon_cell, lat_cell, zone_width in geocells:
            for u in range(n_tiles):
                for v in range(n_tiles):
                    tile_w, tile_s, tile_e, tile_n = tile_bounds(
                        lon_cell, lat_cell, lod, u, v, zone_width
                    )
                    
                    if (tile_e < west or tile_w > east or
                        tile_n < south or tile_s > north):
                        continue
                    
                    try:
                        # Create tile clip box
                        tile_box = box(tile_w, tile_s, tile_e, tile_n)
                        
                        # Clip vector data to tile bounds
                        clipped_gdf = gpd.clip(gdf, tile_box)
                        
                        # Skip empty tiles
                        if len(clipped_gdf) == 0:
                            processed += 1
                            pct = int((processed / total_tiles) * 100)
                            if pct != last_pct:
                                log(f"Progress: {pct}% ({processed}/{total_tiles} tiles)")
                                last_pct = pct
                            continue
                        
                        # Format geocell coordinates
                        if lat_cell >= 0:
                            lat_str = f"N{lat_cell:02d}"
                        else:
                            lat_str = f"S{abs(lat_cell):02d}"
                        
                        if lon_cell >= 0:
                            lon_str = f"E{lon_cell:03d}"
                        else:
                            lon_str = f"W{abs(lon_cell):03d}"
                        
                        dataset_folder = f"{dataset.code:03d}_{dataset.name}"
                        
                        # U folder groups tiles by row (v = latitude index)
                        tile_dir = (
                            cdb_root / "Tiles" / lat_str / lon_str / dataset_folder / 
                            f"L{lod:02d}" / f"U{v}"
                        )
                        
                        tile_dir.mkdir(parents=True, exist_ok=True)
                        
                        # Filename with U (row=v) and R (column=u)
                        filename = f"{lat_str}{lon_str}_D{dataset.code:03d}_S001_T001_L{lod:02d}_U{v}_R{u}.{ext}"
                        output_path = tile_dir / filename
                        
                        log(f"Writing: {output_path.relative_to(cdb_root)} ({len(clipped_gdf)} features)")
                        
                        # Write clipped vector
                        clipped_gdf.to_file(output_path, driver=driver)
                        
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