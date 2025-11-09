#pragma once

#include <string>
#include <vector>
#include <array>
#include <nlohmann/json.hpp>
#include "mesh.h"

// Forward declare GDAL types to avoid exposing in header
class OGRSpatialReference;

namespace aphelion3d::tiles3d {

    enum class ContentFormat { B3DM, GLB };

    struct TilesetOptions {
        std::string outputDir = "./tileset";
        ContentFormat format = ContentFormat::B3DM;
        double geometricError = 50.0;
        
        // CRS handling
        bool autoDetectCRS = true;        // Auto-compute origin from mesh centroid
        
        // Manual origin override (used if autoDetectCRS = false)
        double originLonDeg = 0.0;
        double originLatDeg = 0.0;
        double originHeightM = 0.0;
        
        // If true, apply ENU→ECEF transform in tileset root
        bool useTileTransform = true;
    };

    class TilesetGenerator {
    public:
        explicit TilesetGenerator(const TilesetOptions& options);

        // Generate tileset from building meshes
        // sourceSRS: spatial reference of input meshes (nullptr = assume WGS84)
        bool generate(
            const std::vector<mesh::BuildingMesh>& meshes,
            OGRSpatialReference* sourceSRS = nullptr,
            const std::vector<std::array<double, 3>>& wgs84Positions = {}
        );

    private:
        TilesetOptions options_;

        // Write tileset.json
        bool writeTilesetJSON(
            const nlohmann::json& root,
            const std::string& path
        );

        // Write B3DM tile (header + feature/batch tables + GLB)
        bool writeB3DMTile(
            const std::vector<uint8_t>& glbData,
            const std::string& path
        );

        // Convert PMP meshes to GLB binary
        std::vector<uint8_t> meshesToGLB(
            const std::vector<mesh::BuildingMesh>& meshes
        );
    };

} // namespace aphelion3d::tiles3d

