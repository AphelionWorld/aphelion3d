from __future__ import annotations
from pathlib import Path
from typing import Iterable, Mapping, Any
import fiona
from fiona.crs import CRS
from shapely.geometry import mapping
from shapely.geometry.base import BaseGeometry
from ..datasets import VECTOR, Dataset
from ..tiling import geocell_for_lonlat, tile_index_in_cell
from ..paths import CDBPaths

def write_vector_features(
    cdb: CDBPaths,
    lod: int,
    features: Iterable[tuple[BaseGeometry, Mapping[str, Any]]],
    dataset: Dataset = VECTOR,
    driver: str = "GPKG",
    layer_name: str = "features",
    crs_epsg: int = 4326,
) -> list[Path]:
    buckets: dict[tuple[int,int,int,int,int], list[tuple[BaseGeometry, Mapping[str, Any]]]] = {}
    for geom, props in features:
        c = geom.representative_point()
        lon, lat = c.x, c.y
        lon_cell, lat_cell, zone_width = geocell_for_lonlat(lon, lat)
        u, v = tile_index_in_cell(lon, lat, lod, zone_width)
        key = (lod, lat_cell, lon_cell, u, v)
        buckets.setdefault(key, []).append((geom, props))

    written: list[Path] = []
    for (lod, lat_cell, lon_cell, u, v), feats in buckets.items():
        path = cdb.tile_path(dataset, lod, lat_cell, lon_cell, u, v)
        schema = {
            "geometry": feats[0][0].geom_type,
            "properties": {k: "str" for k in feats[0][1].keys()},
        }
        path.parent.mkdir(parents=True, exist_ok=True)
        with fiona.open(
            path,
            "w",
            driver=driver,
            schema=schema,
            crs=CRS.from_epsg(crs_epsg),
            layer=layer_name if driver.upper() == "GPKG" else None,
        ) as dst:
            for geom, props in feats:
                dst.write({"geometry": mapping(geom), "properties": dict(props)})
        written.append(path)
    return written