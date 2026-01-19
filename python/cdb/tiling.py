from __future__ import annotations
from typing import Tuple
import math

MAX_LOD = 5  # LOD=5 => 32x32 tiles per geocell

def get_zone_width(lat: float) -> int:
    """
    Return the longitude slice width (in degrees) for the CDB zone containing the latitude.
    
    CDB 1.2 Zones (11 total: 0-10):
    - Zone 0:  +89° ≤ lat < +90° → 1° × 12°
    - Zone 1:  +80° ≤ lat < +89° → 1° × 6°
    - Zone 2:  +75° ≤ lat < +80° → 1° × 4°
    - Zone 3:  +70° ≤ lat < +75° → 1° × 3°
    - Zone 4:  +50° ≤ lat < +70° → 1° × 2°
    - Zone 5:  -50° ≤ lat < +50° → 1° × 1°
    - Zone 6:  -70° ≤ lat < -50° → 1° × 2°
    - Zone 7:  -75° ≤ lat < -70° → 1° × 3°
    - Zone 8:  -80° ≤ lat < -75° → 1° × 4°
    - Zone 9:  -89° ≤ lat < -80° → 1° × 6°
    - Zone 10: -90° ≤ lat < -89° → 1° × 12°
    """
    abs_lat = abs(lat)
    
    if abs_lat >= 89.0:
        return 12  # Zone 0 (north) / Zone 10 (south)
    elif abs_lat >= 80.0:
        return 6   # Zone 1 (north) / Zone 9 (south)
    elif abs_lat >= 75.0:
        return 4   # Zone 2 (north) / Zone 8 (south)
    elif abs_lat >= 70.0:
        return 3   # Zone 3 (north) / Zone 7 (south)
    elif abs_lat >= 50.0:
        return 2   # Zone 4 (north) / Zone 6 (south)
    else:
        return 1   # Zone 5 (equatorial, ±50°)

def geocell_for_lonlat(lon: float, lat: float) -> Tuple[int, int, int]:
    """
    Return the geocell (lon_cell, lat_cell, zone_width) containing the point.
    
    Returns:
        lon_cell: Western edge of the geocell (degrees)
        lat_cell: Southern edge of the geocell (degrees)
        zone_width: Width of the zone in longitude degrees (1, 2, 3, 4, 6, or 12)
    """
    lon = (lon + 180.0) % 360.0 - 180.0
    lat = max(min(lat, 89.999999), -89.999999)
    
    zone_width = get_zone_width(lat)
    
    # Use floor to get the southern/western edge of the cell
    lat_cell = math.floor(lat)
    lon_cell = math.floor(lon / zone_width) * zone_width
    
    return lon_cell, lat_cell, zone_width

def tile_index_in_cell(lon: float, lat: float, lod: int, zone_width: int = 1) -> Tuple[int, int]:
    """
    Return tile indices (u, v) within the geocell at the given LOD.
    
    Uses CDB specification formula adapted for zone-aware geocells:
    - dlat = (lat + 90) % 1  (fractional latitude within 1-degree cell)
    - dlon = ((lon + 180) % zone_width) / zone_width (fractional position within zone-width cell)
    - u is the column index (Right/longitude direction): Uref = int(dlon * pow(2, lod))
    - v is the row index (Up/latitude direction): Vref = int(dlat * pow(2, lod))
    
    Args:
        lon: Longitude in degrees
        lat: Latitude in degrees  
        lod: Level of detail
        zone_width: Width of the zone in longitude degrees (1, 2, 3, 4, 6, or 12)
    
    Returns:
        (u, v) where u is column index and v is row index
    """
    if lod < 0 or lod > MAX_LOD:
        raise ValueError(f"lod must be 0..{MAX_LOD}")
    
    n = 1 << lod  # 2^lod
    
    # Latitude is always within 1-degree cells
    dlat = (lat + 90.0) % 1.0
    
    # Longitude is within zone_width-degree cells
    # Normalize to [0, zone_width) then divide by zone_width to get [0, 1)
    dlon = ((lon + 180.0) % zone_width) / zone_width
    
    # u is column (longitude/Right), v is row (latitude/Up)
    u = int(dlon * n)
    v = int(dlat * n)
    
    # Clamp to valid range
    u = min(max(u, 0), n - 1)
    v = min(max(v, 0), n - 1)
    
    return u, v

def tile_bounds(lon_cell: int, lat_cell: int, lod: int, u: int, v: int, zone_width: int = 1) -> Tuple[float, float, float, float]:
    """
    Return (west, south, east, north) bounds of the tile at (u, v) within the geocell.
    
    Args:
        lon_cell: Western edge of geocell
        lat_cell: Southern edge of geocell
        lod: Level of detail
        u: Tile column index (0 to 2^lod - 1)
        v: Tile row index (0 to 2^lod - 1)
        zone_width: Longitude width of the zone (1, 2, 3, 4, 6, or 12 degrees)
    
    Returns:
        (west, south, east, north) in degrees
    """
    n = 1 << lod
    if not (0 <= u < n and 0 <= v < n):
        raise ValueError("u/v out of range for LOD")
    
    # Longitude spacing depends on zone width
    dlon = zone_width / n
    # Latitude is always 1 degree per geocell
    dlat = 1.0 / n
    
    west  = lon_cell + u * dlon
    east  = west + dlon
    south = lat_cell + v * dlat
    north = south + dlat
    
    return west, south, east, north