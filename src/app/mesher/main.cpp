#include <iostream>
#include <string>
#include <filesystem>
#include <pmp/io/io.h>
#include "algorithms.h"
#include "mesh.h"
#include "output_writer.h"
#include "3dtiler.h"

namespace fs = std::filesystem;

void printUsage(const char* programName) {
    std::cout << "Usage: " << programName << " <input_file> <output_directory> [options]\n";
    std::cout << "\n";
    std::cout << "Arguments:\n";
    std::cout << "  input_file       : Path to input vector file (e.g., .shp, .gpkg, .geojson)\n";
    std::cout << "  output_directory : Directory where output files will be written\n";
    std::cout << "\n";
    std::cout << "Options:\n";
    std::cout << "  --obj            : Generate individual OBJ files (default)\n";
    std::cout << "  --single         : Generate single combined OBJ file\n";
    std::cout << "  --groups         : Generate single OBJ with groups (one group per building)\n";
    std::cout << "  --3dtiles        : Generate 3D Tiles (B3DM format)\n";
    std::cout << "  --all            : Generate individual OBJs, grouped OBJ, and 3D Tiles\n";
    std::cout << "\n";
    std::cout << "Examples:\n";
    std::cout << "  " << programName << " buildings.shp output/\n";
    std::cout << "  " << programName << " buildings.shp output/ --groups\n";
    std::cout << "  " << programName << " buildings.shp output/ --3dtiles\n";
    std::cout << "  " << programName << " buildings.shp output/ --all\n";
}

int main(int argc, char** argv) {
    if (argc < 3) {
        printUsage(argv[0]);
        return 1;
    }

    std::string inputFile = argv[1];
    std::string outputDir = argv[2];
    
    // Parse output format options
    bool generateIndividualOBJ = true;  // Default
    bool generateSingleOBJ = false;
    bool generateGroupedOBJ = false;
    bool generate3DTiles = false;
    
    if (argc > 3) {
        std::string option = argv[3];
        if (option == "--obj") {
            generateIndividualOBJ = true;
            generateSingleOBJ = false;
            generateGroupedOBJ = false;
            generate3DTiles = false;
        } else if (option == "--single") {
            generateIndividualOBJ = false;
            generateSingleOBJ = true;
            generateGroupedOBJ = false;
            generate3DTiles = false;
        } else if (option == "--groups") {
            generateIndividualOBJ = false;
            generateSingleOBJ = false;
            generateGroupedOBJ = true;
            generate3DTiles = false;
        } else if (option == "--3dtiles") {
            generateIndividualOBJ = false;
            generateSingleOBJ = false;
            generateGroupedOBJ = false;
            generate3DTiles = true;
        } else if (option == "--all") {
            generateIndividualOBJ = true;
            generateSingleOBJ = false;
            generateGroupedOBJ = true;
            generate3DTiles = true;
        } else {
            std::cerr << "Unknown option: " << option << "\n";
            printUsage(argv[0]);
            return 1;
        }
    }

    // Check if input file exists
    if (!fs::exists(inputFile)) {
        std::cerr << "Error: Input file does not exist: " << inputFile << "\n";
        return 1;
    }

    // Create output directory if it doesn't exist
    try {
        if (!fs::exists(outputDir)) {
            fs::create_directories(outputDir);
            std::cout << "Created output directory: " << outputDir << "\n";
        }
    } catch (const fs::filesystem_error& e) {
        std::cerr << "Error creating output directory: " << e.what() << "\n";
        return 1;
    }

    std::cout << "Reading polygons from: " << inputFile << "\n";
    
    // Generate building meshes from file
    std::vector<aphelion3d::mesh::BuildingMesh> buildingMeshes;
    try {
        aphelion3d::algorithms::generateBuildingMeshesFromFile(inputFile, buildingMeshes);
    } catch (const std::exception& e) {
        std::cerr << "Error reading polygons: " << e.what() << "\n";
        return 1;
    }

    std::cout << "Generated " << buildingMeshes.size() << " building meshes\n\n";

    if (buildingMeshes.empty()) {
        std::cout << "No meshes to write.\n";
        return 0;
    }

    // Get source CRS from input file for 3D Tiles
    OGRSpatialReference* sourceSRS = nullptr;
    GDALDataset* dataset = nullptr;
    
    if (generate3DTiles) {
        GDALAllRegister();
        dataset = (GDALDataset*)GDALOpenEx(
            inputFile.c_str(), GDAL_OF_VECTOR, nullptr, nullptr, nullptr
        );
        
        if (dataset) {
            OGRLayer* layer = dataset->GetLayer(0);
            if (layer) {
                sourceSRS = layer->GetSpatialRef();
                if (sourceSRS) {
                    char* srsWKT = nullptr;
                    sourceSRS->exportToWkt(&srsWKT);
                    std::cout << "Detected CRS: " << (srsWKT ? srsWKT : "unknown") << "\n\n";
                    CPLFree(srsWKT);
                }
            }
        }
    }

    // Generate 3D Tiles if requested
    if (generate3DTiles) {
        std::cout << "=== Generating 3D Tiles ===\n";
        std::string tilesDir = (fs::path(outputDir) / "3dtiles").string();
        
        aphelion3d::tiles3d::TilesetOptions tileOpts;
        tileOpts.outputDir = tilesDir;
        tileOpts.format = aphelion3d::tiles3d::ContentFormat::B3DM;
        tileOpts.geometricError = 50.0;
        tileOpts.autoDetectCRS = true;
        tileOpts.useTileTransform = true;

        aphelion3d::tiles3d::TilesetGenerator generator(tileOpts);
        
        if (generator.generate(buildingMeshes, sourceSRS)) {
            std::cout << "✓ 3D Tiles generated in: " << tilesDir << "\n";
            std::cout << "  Tileset: " << tilesDir << "/tileset.json\n\n";
        } else {
            std::cerr << "✗ Failed to generate 3D Tiles\n\n";
        }
    }
    
    // Generate single OBJ file (no groups)
    if (generateSingleOBJ) {
        std::cout << "=== Generating Single OBJ ===\n";
        std::string singleFile = (fs::path(outputDir) / "all_buildings.obj").string();
        
        if (aphelion3d::output::writeBuildingMeshesToOBJ(buildingMeshes, singleFile, false)) {
            std::cout << "✓ Single OBJ written: " << singleFile << "\n";
            std::cout << "  Total buildings: " << buildingMeshes.size() << "\n\n";
        } else {
            std::cerr << "✗ Failed to write single OBJ file\n\n";
        }
    }
    
    // Generate grouped OBJ file
    if (generateGroupedOBJ) {
        std::cout << "=== Generating Grouped OBJ ===\n";
        std::string groupedFile = (fs::path(outputDir) / "all_buildings_grouped.obj").string();
        
        if (aphelion3d::output::writeBuildingMeshesToOBJ(buildingMeshes, groupedFile, true)) {
            std::cout << "✓ Grouped OBJ written: " << groupedFile << "\n";
            std::cout << "  Total buildings: " << buildingMeshes.size() << "\n";
            std::cout << "  Each building is in its own group (g building_N)\n\n";
        } else {
            std::cerr << "✗ Failed to write grouped OBJ file\n\n";
        }
    }

    // Generate individual OBJ files if requested
    if (generateIndividualOBJ) {
        std::cout << "=== Generating Individual OBJ Files ===\n";
        std::string objDir = (fs::path(outputDir) / "obj").string();
        
        try {
            if (!fs::exists(objDir)) {
                fs::create_directories(objDir);
            }
        } catch (const fs::filesystem_error& e) {
            std::cerr << "Error creating OBJ directory: " << e.what() << "\n";
            if (dataset) GDALClose(dataset);
            return 1;
        }

        int successCount = 0;
        int failCount = 0;
        
        for (size_t i = 0; i < buildingMeshes.size(); ++i) {
            std::string outputFile = (fs::path(objDir) / 
                ("building_" + std::to_string(i) + ".obj")).string();
            
            try {
                pmp::write(buildingMeshes[i], outputFile);
                successCount++;
            } catch (const std::exception& e) {
                std::cerr << "  [" << (i + 1) << "/" << buildingMeshes.size() 
                          << "] Failed to write " << outputFile << ": " 
                          << e.what() << "\n";
                failCount++;
            }
        }

        std::cout << "✓ Individual OBJ files written: " << successCount << "/" 
                  << buildingMeshes.size() << "\n";
        std::cout << "  Output directory: " << objDir << "\n";
        if (failCount > 0) {
            std::cout << "  Failed: " << failCount << " files\n";
        }
        std::cout << "\n";
    }

    // Clean up GDAL
    if (dataset) {
        GDALClose(dataset);
    }

    // Final summary
    std::cout << "=== Summary ===\n";
    std::cout << "Total buildings processed: " << buildingMeshes.size() << "\n";
    std::cout << "Output directory: " << outputDir << "\n";
    if (generateIndividualOBJ) {
        std::cout << "  ├─ obj/ (individual OBJ files)\n";
    }
    if (generateSingleOBJ) {
        std::cout << "  ├─ all_buildings.obj (single combined mesh)\n";
    }
    if (generateGroupedOBJ) {
        std::cout << "  ├─ all_buildings_grouped.obj (grouped by building)\n";
    }
    if (generate3DTiles) {
        std::cout << "  └─ 3dtiles/ (3D Tiles tileset)\n";
    }

    return 0;
}