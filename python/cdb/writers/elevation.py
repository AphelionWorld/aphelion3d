from __future__ import annotations
from pathlib import Path
from typing import Callable
import numpy as np
import rasterio
from rasterio.transform import from_bounds
from rasterio.crs import CRS
from ..paths import CDBPaths
from ..datasets import Dataset, ELEVATION


def write_elevation_tile(
    cdb: CDBPaths,
    lon: float,
    lat: float,
    lod: int,
    size: int,
    sampler: Callable[[np.ndarray, np.ndarray], np.ndarray],
    dataset: Dataset = ELEVATION,
    output_path: Path | None = None,
) -> Path:
    """
    Write an elevation tile to CDB structure.
    
    Args:
        cdb: CDB paths configuration
        lon: Longitude of tile center (degrees)
        lat: Latitude of tile center (degrees)
        lod: Level of detail (0-5)
        size: Tile size in pixels (typically 256)
        sampler: Function that takes (xx, yy) arrays and returns elevations
        dataset: Dataset to write to (default ELEVATION)
        output_path: Optional explicit output path; if None, uses CDB structure
    
    Returns:
        Path to written tile file
    """
    # Create output grid in WGS84
    tile_size_deg = 1.0 / (2 ** lod)
    west = lon - tile_size_deg / 2.0
    south = lat - tile_size_deg / 2.0
    east = lon + tile_size_deg / 2.0
    north = lat + tile_size_deg / 2.0
    
    # Create coordinate grids
    lons = np.linspace(west, east, size)
    lats = np.linspace(south, north, size)
    xx, yy = np.meshgrid(lons, lats[::-1])  # Flip lats for north-up
    
    # Sample elevation data
    elevation_data = sampler(xx, yy)
    
    # Determine output path if not provided
    if output_path is None:
        lon_cell = int(np.floor(lon))
        lat_cell = int(np.floor(lat))
        u = int((lon - lon_cell) * (2 ** lod))
        v = int((lat - lat_cell) * (2 ** lod))
        
        output_path = cdb.tile_path(
            dataset=dataset,
            lod=lod,
            lat_cell=lat_cell,
            lon_cell=lon_cell,
            u=u,
            v=v
        )
    
    # Ensure parent directory exists
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Create geotransform
    transform = from_bounds(west, south, east, north, size, size)
    
    # Write GeoTIFF
    with rasterio.open(
        output_path,
        'w',
        driver='GTiff',
        height=size,
        width=size,
        count=1,
        dtype=elevation_data.dtype,
        crs=CRS.from_epsg(4326),
        transform=transform,
        compress='deflate',
    ) as dst:
        dst.write(elevation_data, 1)
    
    return output_path