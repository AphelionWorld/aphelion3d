from __future__ import annotations
from pathlib import Path
import click
import numpy as np
from shapely.geometry import Point
from collections import defaultdict
import re
import json
import xml.etree.ElementTree as ET
from datetime import datetime

from .paths import CDBPaths
from .raster_tiler import tile_raster_to_cdb, MAX_LOD
from .vector_tiler import tile_vector_to_cdb
from .obj_converter import convert_obj_directory_to_cdb, convert_obj_to_flt_single_lod, convert_obj_to_flt_osg
from .datasets import ELEVATION, IMAGERY, VECTOR, BUILDINGS, ROADS, GSFEATURE, GTFEATURE, get_feature_info
from .metadata import initialize_cdb_metadata
from .writers.elevation import write_elevation_tile
from .writers.vector import write_vector_features


@click.group()
def main():
    """CDB (OGC) command-line utilities for creating and managing CDB datasets."""
    pass


@main.command("init")
@click.option("--root", type=click.Path(path_type=Path), required=True, help="CDB root directory")
@click.option("--version", type=str, default="1.2", show_default=True, help="CDB specification version")
@click.option("--with-metadata", is_flag=True, help="Also generate metadata files")
def init(root: Path, version: str, with_metadata: bool):
    """
    Initialize a CDB root directory structure (OGC CDB compliant).
    
    Supports CDB 1.2 specification with proper folder hierarchy for:
    - Elevation (Dataset 1, 2)
    - Imagery (Dataset 101, 102)
    - Vector Data (Dataset 200-208)
    - Point Features (Dataset 300-302)
    - 3D Models (Dataset 400, 500)
    """
    #(root / "Configuration").mkdir(parents=True, exist_ok=True)
    (root / "Metadata").mkdir(parents=True, exist_ok=True)
    (root / "Tiles").mkdir(parents=True, exist_ok=True)
    click.echo(f"Initialized CDB {version} at {root}")
    
    if with_metadata:
        written = initialize_cdb_metadata(root, version=version)
        click.echo(f"Generated {len(written)} metadata files:")
        for p in written:
            click.echo(f"  {p.relative_to(root)}")


@main.command("tile-raster")
@click.option("--input", "input_path", type=click.Path(exists=True, path_type=Path), required=True, help="Input raster (GeoTIFF, etc.)")
@click.option("--root", type=click.Path(path_type=Path), required=True, help="CDB root directory")
@click.option("--min-lod", type=int, default=0, show_default=True, help="Minimum LOD")
@click.option("--max-lod", type=int, default=5, show_default=True, help="Maximum LOD")
@click.option("--tile-size", type=int, default=256, show_default=True, help="Output tile size in pixels")
def tile_raster_cmd(input_path: Path, root: Path, min_lod: int, max_lod: int, tile_size: int):
    """Tile an elevation raster into CDB tiles at multiple LODs."""
    cdb = CDBPaths(root)
    
    def progress(msg: str):
        click.echo(msg)
    
    try:
        written = tile_raster_to_cdb(
            input_path=input_path,
            cdb=cdb,
            min_lod=min_lod,
            max_lod=max_lod,
            tile_size=tile_size,
            progress_callback=progress,
        )
        click.echo(f"Success! Generated {len(written)} tiles.")
    except Exception as e:
        raise click.ClickException(f"Failed to tile raster: {e}")


@main.command("tile-vector")
@click.option("--input", "input_path", type=click.Path(exists=True, path_type=Path), required=True, help="Input vector file (GeoPackage, Shapefile, etc.)")
@click.option("--root", type=click.Path(path_type=Path), required=True, help="CDB root directory")
@click.option("--dataset", type=click.Choice(["vector", "buildings", "roads"], case_sensitive=False), default="vector", show_default=True, help="Dataset type")
@click.option("--format", type=click.Choice(["gpkg", "shp"], case_sensitive=False), default="gpkg", show_default=True, help="Output format")
@click.option("--min-lod", type=int, default=0, show_default=True, help="Minimum LOD")
@click.option("--max-lod", type=int, default=5, show_default=True, help="Maximum LOD")
def tile_vector_cmd(input_path: Path, root: Path, dataset: str, format: str, min_lod: int, max_lod: int):
    """Tile a vector dataset into CDB tiles at multiple LODs with clipping."""
    
    dataset_map = {
        "vector": VECTOR,
        "buildings": BUILDINGS,
        "roads": ROADS,
    }
    
    dataset_obj = dataset_map[dataset.lower()]
    
    format_map = {
        "gpkg": "GPKG",
        "shp": "ESRI Shapefile",
    }
    output_format = format_map[format.lower()]
    cdb = CDBPaths(root)
    
    def progress(msg: str):
        click.echo(msg)
    
    try:
        written = tile_vector_to_cdb(
            input_path=input_path,
            cdb=cdb,
            min_lod=min_lod,
            max_lod=max_lod,
            dataset=dataset_obj,
            output_format=output_format,
            progress_callback=progress,
        )
        click.echo(f"Success! Generated {len(written)} tiles.")
    except Exception as e:
        raise click.ClickException(f"Failed to tile vector data: {e}")


@main.command("write-elev")
@click.option("--root", type=click.Path(path_type=Path), required=True)
@click.option("--lon", type=float, required=True)
@click.option("--lat", type=float, required=True)
@click.option("--lod", type=int, default=5, show_default=True)
@click.option("--size", type=int, default=256, show_default=True)
def write_elev_cmd(root: Path, lon: float, lat: float, lod: int, size: int):
    """Write one elevation tile using a simple demo sampler."""
    cdb = CDBPaths(root)
    
    def sampler(xx: np.ndarray, yy: np.ndarray) -> np.ndarray:
        r = np.hypot(xx - lon, yy - lat)
        return (100.0 * np.exp(- (r / 0.05) ** 2)).astype(np.float32)
    
    p = write_elevation_tile(cdb, lon, lat, lod, size, sampler)
    click.echo(f"Wrote {p}")


@main.command("write-points")
@click.option("--root", type=click.Path(path_type=Path), required=True)
@click.option("--lon", type=float, multiple=True, required=True, help="Longitude(s)")
@click.option("--lat", type=float, multiple=True, required=True, help="Latitude(s)")
@click.option("--lod", type=int, default=5, show_default=True)
def write_points_cmd(root: Path, lon: tuple[float], lat: tuple[float], lod: int):
    """Write some point features into vector tiles."""
    if len(lon) != len(lat):
        raise click.ClickException("lon and lat must be same length")
    
    feats = [(Point(lon[i], lat[i]), {"name": f"pt{i}"}) for i in range(len(lon))]
    cdb = CDBPaths(root)
    written = write_vector_features(cdb, lod, feats)
    
    for p in written:
        click.echo(f"Wrote {p}")


@main.command("convert-models")
@click.option("--input", "input_path", type=click.Path(exists=True, path_type=Path), required=True, 
              help="Input OBJ file or directory with OBJ files named with LOD suffix (e.g., building_lod0.obj)")
@click.option("--output", "output_path", type=click.Path(path_type=Path), required=True, 
              help="Output directory (for single file) or CDB root (for directory)")
@click.option("--dataset", type=click.Choice(["geospecific", "geotypical"], case_sensitive=False), default="geotypical", show_default=True, help="Model dataset type")
@click.option("--category", type=str, help="Feature category (e.g., A_Culture, A_Hydro). If not provided, will be auto-detected from feature code.")
@click.option("--subcategory", type=str, help="Feature subcategory (e.g., L_Misc_Feature, L_Building). If not provided, will be auto-detected from feature code.")
@click.option("--feature-code", type=str, required=True, help="Feature code (e.g., AL015 for Building, EC030 for Trees, 015 for generic Building)")
@click.option("--lat", type=float, help="Latitude for model placement (geospecific only)")
@click.option("--lon", type=float, help="Longitude for model placement (geospecific only)")
def convert_models_cmd(input_path: Path, output_path: Path, dataset: str, category: str, subcategory: str, feature_code: str, lat: float, lon: float):
    """
    Convert OBJ files to OpenFlight FLT format for CDB models.
    
    The feature code automatically determines category and subcategory if not provided.
    Examples: AL015 (Building), EC030 (Trees), BH065 (Dam), etc.
    
    Supports two modes:
    
    1. Geospecific models (individual placed models):
       cdb convert-models --input model.obj --output model.flt \\
         --dataset geospecific --category A_Culture --subcategory L_Building \\
         --feature-code 015 --lat 32.5 --lon -118.5
    
    2. Geotypical models (reusable feature types):
       cdb convert-models --input obj_dir/ --output cdb_root/ \\
         --dataset geotypical --category A_Culture --subcategory L_Building \\
         --feature-code 015
    
    Directory structure for geotypical models:
       CDB/GTModel/500_GTModelGeometry/A_Culture/L_Building/015_Feature_Name/
    
    Directory structure for geospecific models:
       CDB/Tiles/{lat_cell}/{lon_cell}/400_GeometricModels/L{lod:02d}/U{u}/
    
    For directory mode, filenames should include LOD suffix:
      building_lod0.obj  -> LOD 0
      building_lod1.obj  -> LOD 1
      building_lod5.obj  -> LOD 5
      
    The LOD is extracted from: *_lod{N}.obj
    """
    try:
        from .obj_converter import convert_obj_to_flt_single_lod
        from .datasets import GEOSPECIFIC_MODELS, GEOTYPICAL_MODELS
        
        # Auto-detect category and subcategory from feature code if not provided
        feature_info = get_feature_info(feature_code)
        if feature_info:
            if not category:
                category = feature_info["category"]
                click.echo(f"Auto-detected category: {category}")
            if not subcategory:
                subcategory = feature_info["subcategory"]
                click.echo(f"Auto-detected subcategory: {subcategory}")
        else:
            # Fallback to defaults if feature code not found
            if not category:
                category = "A_Culture"
                click.echo(f"Warning: Feature code '{feature_code}' not found in dictionary. Using default category: {category}")
            if not subcategory:
                subcategory = "L_Misc_Feature"
                click.echo(f"Warning: Feature code '{feature_code}' not found in dictionary. Using default subcategory: {subcategory}")
        
        dataset_map = {
            "geospecific": GEOSPECIFIC_MODELS,
            "geotypical": GEOTYPICAL_MODELS,
        }
        dataset_obj = dataset_map[dataset.lower()]
        
        if input_path.is_file():
            # Single file conversion
            if dataset.lower() == "geospecific" and (lat is None or lon is None):
                raise click.ClickException("--lat and --lon are required for geospecific models")
            
            click.echo(f"Converting {input_path.name} to FLT...")
            success = convert_obj_to_flt_single_lod(input_path, output_path)
            
            if success:
                click.echo(f"✓ Successfully converted to {output_path}")
                if dataset.lower() == "geospecific":
                    click.echo(f"  Location: {lat}, {lon}")
                else:
                    click.echo(f"  Feature: {category}/{subcategory}/{feature_code}")
            else:
                raise click.ClickException(f"Failed to convert {input_path}")
        
        elif input_path.is_dir():
            # Directory conversion
            click.echo(f"Converting OBJ directory to CDB models...")
            
            obj_files = list(input_path.glob("**/*.obj"))
            
            if not obj_files:
                raise click.ClickException(f"No .obj files found in {input_path}")
            
            click.echo(f"Found {len(obj_files)} OBJ files")
            
            # Parse each file for LOD
            file_info = []
            
            for obj_file in obj_files:
                # Parse LOD from filename: *_lod{N}
                lod_match = re.search(r'_lod(\d+)', obj_file.stem, re.IGNORECASE)
                lod = int(lod_match.group(1)) if lod_match else 5  # Default to LOD 5
                
                file_info.append({
                    'path': obj_file,
                    'lod': lod,
                    'name': obj_file.stem,
                })
            
            click.echo(f"\nGrouped by LOD:")
            lod_groups = {}
            for file in file_info:
                lod = file['lod']
                if lod not in lod_groups:
                    lod_groups[lod] = []
                lod_groups[lod].append(file['path'])
            
            # Create the base directory structure
            model_top = "GTModel" if dataset.lower() == "geotypical" else "GSModel"
            geom_folder = "500_GTModelGeometry" if dataset.lower() == "geotypical" else "500_GSModelGeometry"
            base_dir = Path(output_path) / model_top / geom_folder / category / subcategory / f"{feature_code}_building"
            base_dir.mkdir(parents=True, exist_ok=True)
            
            # Group files by base name (without LOD suffix) to determine instances
            instances = {}
            for file_data in file_info:
                # Extract base name by removing _lod{N} suffix
                base_name = re.sub(r'_lod\d+', '', file_data['name'], flags=re.IGNORECASE)
                if base_name not in instances:
                    instances[base_name] = []
                instances[base_name].append(file_data)
            
            click.echo(f"Found {len(instances)} unique model(s)")
            
            instance_count = 1  # Start instance count
            
            # Process each unique model instance
            for base_name, files_for_instance in instances.items():
                # Create a unique subdirectory for this building instance
                instance_dir = base_dir / f"{instance_count:06d}"  # Format as 000001, 000002, etc.
                instance_dir.mkdir(parents=True, exist_ok=True)
                
                click.echo(f"\nInstance {instance_count:06d} ({base_name}):")
                
                for file_data in files_for_instance:
                    lod = file_data['lod']
                    obj_file = file_data['path']
                    
                    # Create LOD subdirectory
                    lod_dir = instance_dir / f"LOD{lod}"
                    lod_dir.mkdir(parents=True, exist_ok=True)
                    
                    # Output filename: model.flt
                    output_filename = "model.flt"
                    output_flt = lod_dir / output_filename
                    
                    click.echo(f"  Converting: {obj_file.name}")
                    click.echo(f"    Output: {output_flt.relative_to(output_path)}")
                    
                    try:
                        success = convert_obj_to_flt_single_lod(obj_file, output_flt)
                        if success:
                            written.append(output_flt)
                            xmlp = write_model_xml(
                                out_flt=output_flt,
                                feature_code=feature_full,
                                seq=seq3 if 'seq3' in locals() else f"{seq:03d}",
                                building_id=bid8 if 'bid8' in locals() else f"{next_building_id:08d}",
                                dataset_prefix=dataset_prefix,
                                series=series,
                                tile=tile,
                                lod_str=lod_str if 'lod_str' in locals() else f"L{lod:02d}",
                                original_name=obj_file.stem
                            )
                            click.echo(f"    ✓ {out_flt.relative_to(target_root)} (xml: {xmlp.relative_to(target_root)})")
                        else:
                            click.echo(f"    ✗ conversion failed for {obj_file.name}")
                    except Exception as e:
                        click.echo(f"    ✗ Error converting {obj_file.name}: {e}")
                
                instance_count += 1  # Increment for next instance
            
            click.echo(f"\n✓ Success! Converted {len(instances)} model instance(s) to {base_dir}")
        
        else:
            raise click.ClickException(f"Input path {input_path} not found")
    
    except Exception as e:
        raise click.ClickException(f"Conversion failed: {e}")


@main.command("run-script")
@click.argument("script", type=click.Path(exists=True, path_type=Path))
@click.option("--root", type=click.Path(path_type=Path), help="Optional CDB root to override entries in the script")
@click.option("--init-root", is_flag=True, help="Create CDB root and basic CDB folders if missing")
def run_script(script: Path, root: Path | None, init_root: bool):
    """
    Execute a JSON script that lists rasters, vectors and models to convert to CDB.

    The script schema is the same as `options.json` example:
      {
        "rasters": [ 
          { "input": "D:/data/elevation.tif", "type": "elevation", "min_lod": 0, "max_lod": 3, "tile_size": 256 },
          { "input": "D:/data/ortho.tif", "type": "imagery", "min_lod": 0, "max_lod": 3, "tile_size": 256 }
        ],
        "vectors": [ 
          { "input": "D:/data/buildings.shp", "dataset": "buildings", "format": "gpkg", "min_lod": 0, "max_lod": 3 } 
        ],
        "models": [
          { "input": "D:/models/buildings", "dataset": "geotypical", "feature_code": "AL015" },
          { "input": "D:/models/trees", "dataset": "geotypical", "feature_code": "EC030", "category": "A_Vegetation", "subcategory": "L_Tree" }
        ]
      }
    
    Note: For models, feature_code is required. Category and subcategory are auto-detected from the feature code 
    dictionary but can be overridden in the JSON. Supported codes: AL015 (Building), EC030 (Trees), BH065 (Dam), etc.

    Examples:
      cdb run-script D:/path/to/options.json --root D:/CDB/test2
      python -m cdb.cli run-script D:/path/to/options.json
    """
    from .raster_tiler import tile_raster_to_cdb
    from .vector_tiler import tile_vector_to_cdb
    from .obj_converter import convert_obj_directory_to_cdb, convert_obj_to_flt_single_lod

    def log(msg: str):
        click.echo(msg)

    def _ensure_root(path: Path):
        if not path.exists():
            if not init_root:
                raise click.ClickException(f"CDB root {path} does not exist (use --init-root to create)")
            click.echo(f"Creating CDB root and basic folders at: {path}")
        
        # Always ensure the CDB directory structure exists when --init-root is set
        if init_root:
            path.mkdir(parents=True, exist_ok=True)
            #(path / "Configuration").mkdir(parents=True, exist_ok=True)
            (path / "Metadata").mkdir(parents=True, exist_ok=True)
            (path / "Tiles").mkdir(parents=True, exist_ok=True)
            # GTModel and GSModel directories are created only when models are processed
            #(path / "GTModel").mkdir(parents=True, exist_ok=True)
            #(path / "GSModel").mkdir(parents=True, exist_ok=True)
            # try to write metadata files (non-fatal)
            try:
                written = initialize_cdb_metadata(path, version="1.2")
                if written:
                    click.echo(f"Initialized metadata ({len(written)} files)")
            except Exception as e:
                click.echo(f"Warning: initialize_cdb_metadata failed: {e}")

    cfg_path = Path(script)
    with open(cfg_path, "r", encoding="utf-8") as fh:
        cfg = json.load(fh)

    cdb_root = Path(root) if root else None

    # Rasters
    for r in cfg.get("rasters", []):
        input_path = Path(r["input"])
        target_root = Path(r.get("root", cdb_root)) if (r.get("root") or cdb_root) else None
        if target_root is None:
            raise click.ClickException("CDB root must be provided either in script or --root")
        _ensure_root(target_root)
        
        # Determine dataset type from the "type" field (elevation or imagery)
        raster_type = r.get("type", "elevation").lower()
        dataset = IMAGERY if raster_type == "imagery" else ELEVATION
        tile_size = int(r.get("tile_size", 256))
        
        log(f"Tiling raster ({raster_type}): {input_path} -> {target_root}")
        log(f"  Tile size: {tile_size}x{tile_size} pixels")
        tile_raster_to_cdb(
            input_path=input_path,
            cdb=CDBPaths(target_root),
            min_lod=int(r.get("min_lod", 0)),
            max_lod=int(r.get("max_lod", 5)),
            tile_size=tile_size,
            dataset=dataset,
            progress_callback=log,
        )

    # Vectors
    for v in cfg.get("vectors", []):
        input_path = Path(v["input"])
        target_root = Path(v.get("root", cdb_root)) if (v.get("root") or cdb_root) else None
        if target_root is None:
            raise click.ClickException("CDB root must be provided either in script or --root")
        _ensure_root(target_root)
        log(f"Tiling vector: {input_path} -> {target_root}")
        
        # Use GSFeature (100) as the standard CDB vector dataset
        # Legacy dataset names like "buildings", "roads" etc. map to GSFeature
        dataset_name = v.get("dataset", "gsfeature").lower()
        if dataset_name in ("buildings", "roads", "vector", "hydrography", "boundaries"):
            vector_dataset = GSFEATURE
            log(f"  Using GSFeature (100) dataset for {dataset_name}")
        elif dataset_name == "gsfeature":
            vector_dataset = GSFEATURE
        elif dataset_name == "gtfeature":
            vector_dataset = GTFEATURE
        else:
            # Try to get from globals, fallback to GSFEATURE
            vector_dataset = globals().get(dataset_name.upper(), None) or GSFEATURE
        
        tile_vector_to_cdb(
            input_path=input_path,
            cdb=CDBPaths(target_root),
            min_lod=int(v.get("min_lod", 0)),
            max_lod=int(v.get("max_lod", 5)),
            dataset=vector_dataset,
            output_format=v.get("format", "GPKG"),
            progress_callback=log,
        )

    # Models (geotypical or geospecific) - unified handling
    for m in cfg.get("models", []):
        input_path = Path(m.get("input", "") or "")
        dataset = m.get("dataset", "geotypical").lower()
        target_root = Path(m.get("root", cdb_root)) if (m.get("root") or cdb_root) else None
        if target_root is None:
            raise click.ClickException("CDB root must be provided either in script or --root")

        _ensure_root(target_root)

        feature_code = m["feature_code"]  # e.g. AL015 or 015
        
        # Extract numeric part for feature lookup (e.g., AL015 -> 015)
        mnum = re.search(r'(\d+)', feature_code)
        feature_num = mnum.group(1) if mnum else feature_code
        
        # Ensure feature number is zero-padded to 3 digits (CDB standard)
        if feature_num.isdigit():
            feature_num = f"{int(feature_num):03d}"
        
        # Auto-detect category and subcategory from feature code if not provided in JSON
        # For GTModels, use the numeric code (015) to get simplified names like "Building"
        # instead of "Building_Superstructure"
        lookup_code = feature_num if dataset == "geotypical" else feature_code
        feature_info = get_feature_info(lookup_code)
        if feature_info:
            category = m.get("category", feature_info["category"])
            subcategory = m.get("subcategory", feature_info["subcategory"])
            feature_name = feature_info["name"]
            log(f"Feature: {feature_code} - {feature_name} ({category}/{subcategory})")
        else:
            category = m.get("category", "A_Culture")
            subcategory = m.get("subcategory", "L_Misc_Feature")
            feature_name = "Building"  # default
            log(f"Warning: Feature code '{feature_code}' not found in dictionary. Using defaults.")
        
        feature_full = feature_code

        # base: GTModel/500_GTModelGeometry/<Category>/<Subcategory>/<015_Building>/
        model_top = "GTModel" if dataset == "geotypical" else "GSModel"
        geom_folder = "500_GTModelGeometry" if dataset == "geotypical" else "500_GSModelGeometry"
        gt_base = Path(target_root) / model_top / geom_folder / category / subcategory / f"{feature_num}_{feature_name}"
        gt_base.mkdir(parents=True, exist_ok=True)
        click.echo(f"Models -> {dataset}: {input_path} -> {gt_base.resolve()}")

        # collect OBJ files
        objs: list[Path] = []
        if input_path and Path(input_path).exists():
            p = Path(input_path)
            objs = list(p.glob("**/*.obj")) if p.is_dir() else [p]
        if not objs:
            click.echo(f"  Warning: no OBJ files found for model entry: {m.get('input')}")
            continue

        # build file_info (path, name, lod) and group by base name (without _lodN)
        file_info: list[dict] = []
        for obj in objs:
            stem = obj.stem
            lod_match = re.search(r'_lod(\d+)', stem, re.IGNORECASE)
            lod = int(lod_match.group(1)) if lod_match else 0
            file_info.append({"path": obj, "name": stem, "lod": lod})

        instances: dict[str, list[dict]] = {}
        for f in file_info:
            base_name = re.sub(r'_lod\d+$', '', f['name'], flags=re.IGNORECASE)
            instances.setdefault(base_name, []).append(f)

        # determine next building id by scanning existing filenames under gt_base
        existing_ids: list[int] = []
        for p in gt_base.rglob("*.flt"):
            mid = re.search(r'_building_(\d{1,})\.flt$', p.name)
            if mid:
                try:
                    existing_ids.append(int(mid.group(1)))
                except Exception:
                    pass
        next_building_id = (max(existing_ids) + 1) if existing_ids else 1

        # choose dataset prefix: geospecific -> D500, geotypical -> D510 (adjustable)
        dataset_prefix = "D500" if dataset == "geospecific" else "D510"
        series = m.get("series", "S001")
        tile = m.get("tile", "T001")

        seq = 1
        written = []
        for base_name, files_for_instance in instances.items():
            seq3 = f"{seq:03d}"
            bid8 = f"{next_building_id:08d}"
            short_name = re.sub(r'[^a-z0-9\-_]', '_', base_name.lower())

            for file_data in files_for_instance:
                lod = file_data["lod"]
                obj_file = file_data["path"]
                lod_str = f"L{lod:02d}"        # L00, L01...
                lod_dir = gt_base / lod_str
                lod_dir.mkdir(parents=True, exist_ok=True)

                # filename: D510_S001_T001_L00_AL015_004_building_00000001.flt
                out_fname = f"{dataset_prefix}_{series}_{tile}_{lod_str}_{feature_full}_{seq3}_building_{bid8}.flt"
                out_flt = lod_dir / out_fname

                click.echo(f"  Converting {obj_file.name} -> {out_flt.resolve()}")
                try:
                    # Try OSG converter first, fallback to custom writer
                    success = convert_obj_to_flt_osg(obj_file, out_flt)
                    if not success:
                        click.echo(f"    OSG conversion failed, trying custom FLT writer...")
                        success = convert_obj_to_flt_single_lod(obj_file, out_flt)
                    
                    if success:
                        written.append(out_flt)
                        xmlp = write_model_xml(
                            out_flt=out_flt,
                            feature_code=feature_full,
                            seq=seq3 if 'seq3' in locals() else f"{seq:03d}",
                            building_id=bid8 if 'bid8' in locals() else f"{next_building_id:08d}",
                            dataset_prefix=dataset_prefix,
                            series=series,
                            tile=tile,
                            lod_str=lod_str if 'lod_str' in locals() else f"L{lod:02d}",
                            original_name=obj_file.stem
                        )
                        click.echo(f"    ✓ {out_flt.relative_to(target_root)} (xml: {xmlp.relative_to(target_root)})")
                    else:
                        click.echo(f"    ✗ conversion failed for {obj_file.name}")
                except Exception as e:
                    click.echo(f"    ✗ Error converting {obj_file.name}: {e}")

            seq += 1
            next_building_id += 1

        click.echo(f"  Converted {len(written)} model files for feature {feature_full} into {gt_base.resolve()}")
    
    log("Script complete.")

def write_model_xml(out_flt: Path,
                    feature_code: str,
                    seq: str,
                    building_id: str,
                    dataset_prefix: str,
                    series: str,
                    tile: str,
                    lod_str: str,
                    original_name: str | None = None) -> Path:
    """
    Write a simple XML metadata file next to a FLT file.
    Returns the Path to the xml file.
    """
    xml_path = out_flt.with_suffix('.xml')
    root = ET.Element("ModelMetadata")
    ET.SubElement(root, "Filename").text = out_flt.name
    ET.SubElement(root, "FeatureCode").text = feature_code
    ET.SubElement(root, "DatasetPrefix").text = dataset_prefix
    ET.SubElement(root, "Series").text = series
    ET.SubElement(root, "Tile").text = tile
    ET.SubElement(root, "LOD").text = lod_str
    ET.SubElement(root, "Sequence").text = seq
    ET.SubElement(root, "BuildingID").text = building_id
    if original_name:
        ET.SubElement(root, "OriginalName").text = original_name
    ET.SubElement(root, "CreatedBy").text = "Aphelion3D"
    ET.SubElement(root, "CreatedOn").text = datetime.utcnow().isoformat() + "Z"

    xml_path.parent.mkdir(parents=True, exist_ok=True)
    tree = ET.ElementTree(root)
    tree.write(xml_path, encoding="utf-8", xml_declaration=True)
    return xml_path