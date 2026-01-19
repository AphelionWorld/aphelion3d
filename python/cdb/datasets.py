from __future__ import annotations
from dataclasses import dataclass


@dataclass
class Dataset:
    """CDB Dataset definition"""
    code: int
    name: str
    extension: str


# CDB 1.2 Dataset Codes (OGC 15-113r6)

# Elevation Data (1-99)
ELEVATION = Dataset(code=1, name="Elevation", extension="tif")
MIN_MAX_ELEVATION = Dataset(code=2, name="MinMaxElevation", extension="tif")
TERRAIN_ELEVATION = Dataset(code=10, name="TerrainElevation", extension="tif")

# Raster Imagery (1-99)
IMAGERY = Dataset(code=4, name="Imagery", extension="tif")
IMAGERY_RGB = Dataset(code=4, name="Imagery_RGB", extension="tif")
IMAGERY_RGBA = Dataset(code=5, name="Imagery_RGBA", extension="tif")

# Vector Features (100-199)
GSFEATURE = Dataset(code=100, name="GSFeature", extension="shp")  # Geospatial features (point/line/polygon)
GTFEATURE = Dataset(code=101, name="GTFeature", extension="shp")  # Geotypical features

# Vector Data (200-299) - Legacy/alternate datasets
VECTOR = Dataset(code=200, name="VectorData", extension="gpkg")
ROADS = Dataset(code=201, name="Roads", extension="gpkg")
RAILROADS = Dataset(code=202, name="Railroads", extension="gpkg")
POWER_LINES = Dataset(code=203, name="PowerLines", extension="gpkg")
HYDROGRAPHY = Dataset(code=204, name="Hydrography", extension="gpkg")
BOUNDARIES = Dataset(code=205, name="Boundaries", extension="gpkg")
BUILDINGS = Dataset(code=206, name="Buildings", extension="gpkg")
LANDUSE = Dataset(code=207, name="LandUse", extension="gpkg")
TREES = Dataset(code=208, name="Trees", extension="gpkg")

# Geometric & Attribute Data (300-399)
POINT_FEATURES = Dataset(code=300, name="PointFeatures", extension="shp")
TEXT = Dataset(code=301, name="Text", extension="shp")
RASTER_ATTRIBUTION = Dataset(code=302, name="RasterAttribution", extension="tif")

# Model Data (400-599)
GEOSPECIFIC_MODELS = Dataset(code=400, name="GeometricModels", extension="flt")
GEOTYPICAL_MODELS = Dataset(code=500, name="GeotypicalModels", extension="flt")

# Metadata (600+)
METADATA = Dataset(code=600, name="Metadata", extension="xml")
FEATURE_DATA_DICTIONARY = Dataset(code=601, name="FeatureDataDictionary", extension="xml")

# Legacy/Compatibility
MODELS = Dataset(code=400, name="Models", extension="flt")  # Alias for GEOSPECIFIC_MODELS


# Feature Type Codes (OGC CDB 1.2 / FACC/DFDD)
# These are used in model naming and feature classification
FEATURE_CODES = {
    # Buildings & Structures (AL0XX)
    "AL013": {"name": "Building", "category": "A_Culture", "subcategory": "L_Building"},
    "AL015": {"name": "Building_Superstructure", "category": "A_Culture", "subcategory": "L_Building"},
    "AL018": {"name": "Building_Shed", "category": "A_Culture", "subcategory": "L_Building"},
    "AL020": {"name": "Built_Up_Area", "category": "A_Culture", "subcategory": "L_Settlement"},
    "AL030": {"name": "Cemetery", "category": "A_Culture", "subcategory": "L_Misc_Feature"},
    "AL036": {"name": "Castle", "category": "A_Culture", "subcategory": "L_Building"},
    "AL045": {"name": "Complex_Outline", "category": "A_Culture", "subcategory": "L_Misc_Feature"},
    "AL070": {"name": "Fence", "category": "A_Culture", "subcategory": "L_Misc_Feature"},
    "AL073": {"name": "Flagpole", "category": "A_Culture", "subcategory": "L_Misc_Feature"},
    "AL080": {"name": "Gantry", "category": "A_Culture", "subcategory": "L_Misc_Feature"},
    "AL099": {"name": "Hut", "category": "A_Culture", "subcategory": "L_Building"},
    "AL105": {"name": "Settlement", "category": "A_Culture", "subcategory": "L_Settlement"},
    "AL130": {"name": "Monument", "category": "A_Culture", "subcategory": "L_Misc_Feature"},
    "AL200": {"name": "Ruins", "category": "A_Culture", "subcategory": "L_Building"},
    "AL208": {"name": "Shanty_Town", "category": "A_Culture", "subcategory": "L_Settlement"},
    "AL241": {"name": "Tower", "category": "A_Culture", "subcategory": "L_Building"},
    "AL260": {"name": "Wall", "category": "A_Culture", "subcategory": "L_Misc_Feature"},
    
    # Transportation (AP0XX, AQ0XX)
    "AP010": {"name": "Cart_Track", "category": "A_Culture", "subcategory": "L_Road"},
    "AP020": {"name": "Interchange", "category": "A_Culture", "subcategory": "L_Road"},
    "AP030": {"name": "Road", "category": "A_Culture", "subcategory": "L_Road"},
    "AP040": {"name": "Gate", "category": "A_Culture", "subcategory": "L_Misc_Feature"},
    "AP041": {"name": "Vehicle_Barrier", "category": "A_Culture", "subcategory": "L_Misc_Feature"},
    "AP050": {"name": "Trail", "category": "A_Culture", "subcategory": "L_Road"},
    "AQ040": {"name": "Bridge", "category": "A_Culture", "subcategory": "L_Bridge"},
    "AQ045": {"name": "Bridge_Span", "category": "A_Culture", "subcategory": "L_Bridge"},
    "AQ050": {"name": "Bridge_Superstructure", "category": "A_Culture", "subcategory": "L_Bridge"},
    "AQ055": {"name": "Bridge_Tower", "category": "A_Culture", "subcategory": "L_Bridge"},
    "AQ060": {"name": "Control_Tower", "category": "A_Culture", "subcategory": "L_Building"},
    "AQ065": {"name": "Culvert", "category": "A_Culture", "subcategory": "L_Misc_Feature"},
    "AQ113": {"name": "Pipeline", "category": "A_Culture", "subcategory": "L_Misc_Feature"},
    "AQ125": {"name": "Transportation_Station", "category": "A_Culture", "subcategory": "L_Building"},
    "AQ130": {"name": "Tunnel", "category": "A_Culture", "subcategory": "L_Misc_Feature"},
    
    # Utilities (AT0XX)
    "AT005": {"name": "Cable", "category": "A_Culture", "subcategory": "L_Misc_Feature"},
    "AT010": {"name": "Dish_Aerial", "category": "A_Culture", "subcategory": "L_Misc_Feature"},
    "AT042": {"name": "Pylon", "category": "A_Culture", "subcategory": "L_Misc_Feature"},
    "AT045": {"name": "Radar_Station", "category": "A_Culture", "subcategory": "L_Building"},
    "AT080": {"name": "Communication_Tower", "category": "A_Culture", "subcategory": "L_Building"},
    
    # Hydrography (BH0XX)
    "BH020": {"name": "Canal", "category": "A_Hydro", "subcategory": "L_Water"},
    "BH040": {"name": "Water_Treatment_Facility", "category": "A_Culture", "subcategory": "L_Building"},
    "BH065": {"name": "Dam", "category": "A_Hydro", "subcategory": "L_Misc_Feature"},
    "BH070": {"name": "Ford", "category": "A_Hydro", "subcategory": "L_Water"},
    "BH082": {"name": "Inland_Water_Body", "category": "A_Hydro", "subcategory": "L_Water"},
    "BH090": {"name": "Land_Subject_to_Inundation", "category": "A_Hydro", "subcategory": "L_Water"},
    "BH110": {"name": "Moat", "category": "A_Hydro", "subcategory": "L_Water"},
    "BH140": {"name": "River", "category": "A_Hydro", "subcategory": "L_Water"},
    "BH145": {"name": "Vanishing_Point", "category": "A_Hydro", "subcategory": "L_Water"},
    "BH170": {"name": "Natural_Pool", "category": "A_Hydro", "subcategory": "L_Water"},
    
    # Vegetation (EC0XX, ED0XX)
    "EC015": {"name": "Forest", "category": "A_Vegetation", "subcategory": "L_Tree"},
    "EC030": {"name": "Trees", "category": "A_Vegetation", "subcategory": "L_Tree"},
    "EC040": {"name": "Cleared_Way", "category": "A_Vegetation", "subcategory": "L_Misc_Feature"},
    "ED010": {"name": "Marsh", "category": "A_Vegetation", "subcategory": "L_Misc_Feature"},
    "ED020": {"name": "Swamp", "category": "A_Vegetation", "subcategory": "L_Misc_Feature"},
    
    # Generic/Miscellaneous
    "015": {"name": "Building", "category": "A_Culture", "subcategory": "L_Building"},
    "100": {"name": "Tree", "category": "A_Vegetation", "subcategory": "L_Tree"},
    "999": {"name": "Generic_Feature", "category": "A_Culture", "subcategory": "L_Misc_Feature"},
}


def get_feature_info(feature_code: str) -> dict[str, str] | None:
    """
    Look up feature information by code.
    
    Args:
        feature_code: Feature code (e.g., "AL015", "015")
    
    Returns:
        Dictionary with name, category, and subcategory, or None if not found
    """
    return FEATURE_CODES.get(feature_code.upper())


def get_dataset_by_name(name: str) -> Dataset | None:
    """
    Look up a dataset by its name.
    
    Args:
        name: Dataset name (e.g., "elevation", "imagery", "buildings")
    
    Returns:
        Dataset object or None if not found
    """
    name_lower = name.lower()
    dataset_map = {
        "elevation": ELEVATION,
        "min_max_elevation": MIN_MAX_ELEVATION,
        "terrain_elevation": TERRAIN_ELEVATION,
        "imagery": IMAGERY,
        "imagery_rgb": IMAGERY_RGB,
        "imagery_rgba": IMAGERY_RGBA,
        "vector": VECTOR,
        "vectordata": VECTOR,
        "roads": ROADS,
        "railroads": RAILROADS,
        "power_lines": POWER_LINES,
        "powerlines": POWER_LINES,
        "hydrography": HYDROGRAPHY,
        "boundaries": BOUNDARIES,
        "buildings": BUILDINGS,
        "landuse": LANDUSE,
        "trees": TREES,
        "point_features": POINT_FEATURES,
        "pointfeatures": POINT_FEATURES,
        "text": TEXT,
        "raster_attribution": RASTER_ATTRIBUTION,
        "geospecific": GEOSPECIFIC_MODELS,
        "geospecific_models": GEOSPECIFIC_MODELS,
        "geotypical": GEOTYPICAL_MODELS,
        "geotypical_models": GEOTYPICAL_MODELS,
        "models": MODELS,
        "metadata": METADATA,
        "feature_data_dictionary": FEATURE_DATA_DICTIONARY,
    }
    return dataset_map.get(name_lower)