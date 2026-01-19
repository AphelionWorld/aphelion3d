#!/usr/bin/env python3
"""
Test script to demonstrate the feature code dictionary usage.
"""

from cdb.datasets import (
    get_feature_info, FEATURE_CODES,
    ELEVATION, MIN_MAX_ELEVATION, TERRAIN_ELEVATION,
    IMAGERY, IMAGERY_RGB, IMAGERY_RGBA,
    VECTOR, ROADS, RAILROADS, POWER_LINES, HYDROGRAPHY, BOUNDARIES, BUILDINGS, LANDUSE, TREES,
    POINT_FEATURES, TEXT, RASTER_ATTRIBUTION,
    GEOSPECIFIC_MODELS, GEOTYPICAL_MODELS,
    METADATA, FEATURE_DATA_DICTIONARY
)

def test_dataset_codes():
    """Display all CDB dataset codes."""
    
    print("=" * 60)
    print("CDB Dataset Codes")
    print("=" * 60)
    
    datasets = [
        ("Elevation Data", [
            ELEVATION,
            MIN_MAX_ELEVATION,
            TERRAIN_ELEVATION,
        ]),
        ("Raster Imagery", [
            IMAGERY,
            IMAGERY_RGB,
            IMAGERY_RGBA,
        ]),
        ("Vector Data", [
            VECTOR,
            ROADS,
            RAILROADS,
            POWER_LINES,
            HYDROGRAPHY,
            BOUNDARIES,
            BUILDINGS,
            LANDUSE,
            TREES,
        ]),
        ("Geometric & Attribute Data", [
            POINT_FEATURES,
            TEXT,
            RASTER_ATTRIBUTION,
        ]),
        ("Model Data", [
            GEOSPECIFIC_MODELS,
            GEOTYPICAL_MODELS,
        ]),
        ("Metadata", [
            METADATA,
            FEATURE_DATA_DICTIONARY,
        ]),
    ]
    
    for category, dataset_list in datasets:
        print(f"\n{category}:")
        for ds in dataset_list:
            print(f"  {ds.code:3d} - {ds.name:30s} (.{ds.extension})")

def test_feature_codes():
    """Test various feature codes."""
    
    print("\n" + "=" * 60)
    print("CDB Feature Code Dictionary Test")
    print("=" * 60)
    
    # Test common feature codes
    test_codes = [
        "AL015",  # Building
        "EC030",  # Trees
        "BH065",  # Dam
        "AP030",  # Road
        "AQ040",  # Bridge
        "015",    # Generic building
        "100",    # Generic tree
        "INVALID", # Should return None
    ]
    
    for code in test_codes:
        info = get_feature_info(code)
        if info:
            print(f"\n{code}:")
            print(f"  Name: {info['name']}")
            print(f"  Category: {info['category']}")
            print(f"  Subcategory: {info['subcategory']}")
        else:
            print(f"\n{code}: NOT FOUND")
    
    print("\n" + "=" * 60)
    print(f"Total feature codes in dictionary: {len(FEATURE_CODES)}")
    print("=" * 60)
    
    # List all available feature codes
    print("\nAll available feature codes:")
    for code in sorted(FEATURE_CODES.keys()):
        info = FEATURE_CODES[code]
        print(f"  {code:8} - {info['name']}")

if __name__ == "__main__":
    test_dataset_codes()
    print("\n")
    test_feature_codes()
