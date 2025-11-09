#include "3dtiler.h"
#include <cstring>
#include <fstream>
#include <filesystem>
#include <iostream>
#include <limits>
#include <cmath>
#include <tiny_gltf.h>
#include <ogr_spatialref.h>

namespace fs = std::filesystem;
using json = nlohmann::json;

namespace aphelion3d::tiles3d {

namespace {

// WGS84 to ECEF conversion
std::array<double, 3> wgs84ToECEF(double latDeg, double lonDeg, double heightM) {
    const double a = 6378137.0;              // WGS84 semi-major axis
    const double f = 1.0 / 298.257223563;    // WGS84 flattening
    const double e2 = 2.0 * f - f * f;       // First eccentricity squared
    
    const double latRad = latDeg * M_PI / 180.0;
    const double lonRad = lonDeg * M_PI / 180.0;
    
    const double sinLat = std::sin(latRad);
    const double cosLat = std::cos(latRad);
    const double N = a / std::sqrt(1.0 - e2 * sinLat * sinLat);
    
    return {
        (N + heightM) * cosLat * std::cos(lonRad),
        (N + heightM) * cosLat * std::sin(lonRad),
        (N * (1.0 - e2) + heightM) * sinLat
    };
}

// ENU basis vectors at given lat/lon
void enuAxes(double latRad, double lonRad,
             std::array<double,3>& E,
             std::array<double,3>& N,
             std::array<double,3>& U) {
    const double sinLat = std::sin(latRad), cosLat = std::cos(latRad);
    const double sinLon = std::sin(lonRad), cosLon = std::cos(lonRad);
    E = { -sinLon,            cosLon,            0.0 };
    N = { -sinLat * cosLon,  -sinLat * sinLon,   cosLat };
    U = {  cosLat * cosLon,   cosLat * sinLon,   sinLat };
}

// Build 4x4 ENU→ECEF transform matrix (column-major for Cesium)
std::array<double,16> enuToEcefMatrix(double latDeg, double lonDeg, double heightM) {
    const double deg2rad = M_PI / 180.0;
    std::array<double,3> E, N, U;
    enuAxes(latDeg * deg2rad, lonDeg * deg2rad, E, N, U);
    
    // Use GDAL to transform origin to ECEF
    OGRSpatialReference wgs84;
    wgs84.SetWellKnownGeogCS("WGS84");
    
    OGRSpatialReference ecef;
    ecef.importFromEPSG(4978);
    
    OGRCoordinateTransformation* toEcef = OGRCreateCoordinateTransformation(&wgs84, &ecef);
    
    double ox = lonDeg, oy = latDeg, oz = heightM;
    if (toEcef && toEcef->Transform(1, &ox, &oy, &oz)) {
        OCTDestroyCoordinateTransformation(toEcef);
    } else {
        std::cerr << "Warning: GDAL ECEF transform failed, using manual calculation\n";
        if (toEcef) OCTDestroyCoordinateTransformation(toEcef);
        auto O = wgs84ToECEF(latDeg, lonDeg, heightM);  // fallback
        ox = O[0]; oy = O[1]; oz = O[2];
    }

    // Column-major: [ E N U T ]
    return {
        E[0], N[0], U[0], 0.0,
        E[1], N[1], U[1], 0.0,
        E[2], N[2], U[2], 0.0,
        ox, oy, oz, 1.0
    };
}

// Compute centroid of all meshes in WGS84
std::array<double, 3> computeMeshCentroidWGS84(
    const std::vector<mesh::BuildingMesh>& meshes,
    OGRCoordinateTransformation* toWGS84
) {
    double sumLon = 0.0, sumLat = 0.0, sumH = 0.0;
    size_t count = 0;

    for (const auto& m : meshes) {
        for (auto v : m.vertices()) {
            auto pos = m.position(v);
            double x = pos[0], y = pos[1], z = pos[2];
            
            if (toWGS84) {
                if (toWGS84->Transform(1, &x, &y, &z)) {
                    sumLat += x;  // GDAL returns lat first
                    sumLon += y;  // then lon
                    sumH += z;
                    ++count;
                }
            } else {
                // Assume already in WGS84
                sumLon += x;
                sumLat += y;
                sumH += z;
                ++count;
            }
        }
    }

    if (count == 0) return {0.0, 0.0, 0.0};
    return {sumLon / count, sumLat / count, sumH / count};  // Return as lon, lat
}

// Transform mesh vertices from source CRS to local ENU
void transformMeshToENU(
    mesh::BuildingMesh& m,
    OGRCoordinateTransformation* toWGS84,
    double originLon, double originLat, double originH
) {
    // Create WGS84 and ECEF coordinate systems
    OGRSpatialReference wgs84;
    wgs84.SetWellKnownGeogCS("WGS84");
    wgs84.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);  // Force lon, lat order
    
    OGRSpatialReference ecef;
    ecef.importFromEPSG(4978);  // EPSG:4978 is WGS84 ECEF
   // ecef.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);  // Force X, Y, Z order
    
    OGRCoordinateTransformation* wgs84ToEcef = OGRCreateCoordinateTransformation(&wgs84, &ecef);
    if (!wgs84ToEcef) {
        std::cerr << "Failed to create WGS84->ECEF transformation\n";
        return;
    }
    
    // Transform origin to ECEF
    double ox = originLon, oy = originLat, oz = originH;
    if (!wgs84ToEcef->Transform(1, &ox, &oy, &oz)) {
        std::cerr << "Failed to transform origin to ECEF\n";
        OCTDestroyCoordinateTransformation(wgs84ToEcef);
        return;
    }
    
    std::cout << "Origin ECEF: X=" << ox << " Y=" << oy << " Z=" << oz << "\n";
    
    // Compute ENU axes at origin
    const double deg2rad = M_PI / 180.0;
    std::array<double,3> E, N, U;
    enuAxes(originLat * deg2rad, originLon * deg2rad, E, N, U);

    for (auto v : m.vertices()) {
        auto pos = m.position(v);
        double x = pos[0], y = pos[1], z = pos[2];
        
        // Transform to WGS84
        if (toWGS84) {
            if (!toWGS84->Transform(1, &x, &y, &z)) {
                std::cerr << "Warning: coordinate transformation failed for vertex\n";
                continue;
            }
        }
        
        // Transform WGS84 to ECEF using GDAL
        if (!wgs84ToEcef->Transform(1, &x, &y, &z)) {
            std::cerr << "Warning: WGS84->ECEF transformation failed\n";
            continue;
        }
        
        // ECEF → ENU: local = [E N U]^T · (P - O)
        const double dx = x - ox;
        const double dy = y - oy;
        const double dz = z - oz;
        
        const pmp::Point enuPos(
            E[0]*dx + E[1]*dy + E[2]*dz,
            N[0]*dx + N[1]*dy + N[2]*dz,
            U[0]*dx + U[1]*dy + U[2]*dz
        );
        
        m.position(v) = enuPos;
    }
    
    OCTDestroyCoordinateTransformation(wgs84ToEcef);
}

// Helper: align buffer size to 4 bytes (GLTF requirement)
void alignTo4(std::vector<uint8_t>& buf) {
    while (buf.size() % 4 != 0) buf.push_back(0);
}

} // anonymous namespace

TilesetGenerator::TilesetGenerator(const TilesetOptions& options)
    : options_(options) {
    fs::create_directories(options_.outputDir);
    fs::create_directories(options_.outputDir + "/tiles");
}

bool TilesetGenerator::generate(
    const std::vector<mesh::BuildingMesh>& meshes,
    OGRSpatialReference* sourceSRS,
    const std::vector<std::array<double, 3>>& /*wgs84Positions*/
) {
    if (meshes.empty()) return false;

    // Determine origin and setup transformation
    double originLon = options_.originLonDeg;
    double originLat = options_.originLatDeg;
    double originH = options_.originHeightM;
    
    OGRCoordinateTransformation* toWGS84 = nullptr;
    
    if (sourceSRS && options_.autoDetectCRS) {
        OGRSpatialReference wgs84;
        wgs84.SetWellKnownGeogCS("WGS84");
        toWGS84 = OGRCreateCoordinateTransformation(sourceSRS, &wgs84);
        
        if (!toWGS84) {
            std::cerr << "Warning: failed to create coordinate transformation\n";
        }
    }
    
    // Copy meshes so we can transform them without modifying originals
    std::vector<mesh::BuildingMesh> transformedMeshes = meshes;
    
    if (options_.useTileTransform) {
        // Compute centroid as origin if auto-detect enabled
        if (options_.autoDetectCRS) {
            auto centroid = computeMeshCentroidWGS84(transformedMeshes, toWGS84);
            originLon = centroid[0];
            originLat = centroid[1];
            originH = centroid[2];
            
            std::cout << "Auto-detected origin: lon=" << originLon 
                      << " lat=" << originLat << " h=" << originH << "m\n";
        }
        
        // Transform all meshes to ENU
        for (auto& m : transformedMeshes) {
            transformMeshToENU(m, toWGS84, originLon, originLat, originH);
        }
    }

    std::cout << "Origin coordinates:\n";
    std::cout << "  Longitude: " << originLon << " deg\n";
    std::cout << "  Latitude: " << originLat << " deg\n";
    std::cout << "  Height: " << originH << " m\n";
    
    if (toWGS84) {
        OCTDestroyCoordinateTransformation(toWGS84);
    }

    // Build GLB from transformed meshes
    const auto glbData = meshesToGLB(transformedMeshes);
    if (glbData.empty()) return false;

    // Build tileset.json
    json tileset;
    tileset["asset"] = { {"version","1.0"}, {"generator","Aphelion3D"} };
    tileset["geometricError"] = options_.geometricError;

    json root;
    root["geometricError"] = options_.geometricError / 2.0;
    root["refine"] = "REPLACE";

    // Set ENU→ECEF transform
    if (options_.useTileTransform) {
        const auto M = enuToEcefMatrix(originLat, originLon, originH);
        root["transform"] = std::vector<double>(M.begin(), M.end());
    }

    // Simple box bounding volume (you can enhance this)
    // For now, omit or use a large default
    root["boundingVolume"]["box"] = {
        0.0, 0.0, 0.0,  // center in ENU
        1000.0, 0.0, 0.0,
        0.0, 1000.0, 0.0,
        0.0, 0.0, 100.0
    };

    std::string contentRel = (options_.format == ContentFormat::B3DM) 
        ? "tiles/tile_0.b3dm" : "tiles/tile_0.glb";
    root["content"] = { {"uri", contentRel} };
    tileset["root"] = root;

    fs::create_directories(options_.outputDir);
    if (!writeTilesetJSON(tileset, options_.outputDir + "/tileset.json"))
        return false;

    const std::string outPath = options_.outputDir + "/" + contentRel;
    if (options_.format == ContentFormat::B3DM) {
        return writeB3DMTile(glbData, outPath);
    } else {
        std::ofstream f(outPath, std::ios::binary);
        if (!f) return false;
        f.write(reinterpret_cast<const char*>(glbData.data()), glbData.size());
        return true;
    }
}

bool TilesetGenerator::writeTilesetJSON(const json& root, const std::string& path) {
    std::ofstream file(path);
    if (!file) return false;
    file << root.dump(2);
    return true;
}

bool TilesetGenerator::writeB3DMTile(const std::vector<uint8_t>& glbData, const std::string& path) {
    const uint32_t headerSize = 28;

    json ft;
    ft["BATCH_LENGTH"] = 1;
    std::string ftJSON = ft.dump();
    while (ftJSON.size() % 8 != 0) ftJSON.push_back(' ');

    json bt = json::object();
    std::string btJSON = bt.dump();
    while (btJSON.size() % 8 != 0) btJSON.push_back(' ');

    const uint32_t featureTableJSONSize = static_cast<uint32_t>(ftJSON.size());
    const uint32_t featureTableBinarySize = 0;
    const uint32_t batchTableJSONSize = static_cast<uint32_t>(btJSON.size());
    const uint32_t batchTableBinarySize = 0;

    const uint32_t totalSize =
        headerSize +
        featureTableJSONSize + featureTableBinarySize +
        batchTableJSONSize + batchTableBinarySize +
        static_cast<uint32_t>(glbData.size());

    std::ofstream out(path, std::ios::binary);
    if (!out) return false;

    out.write("b3dm", 4);
    uint32_t version = 1;
    out.write(reinterpret_cast<char*>(&version), 4);
    out.write(reinterpret_cast<const char*>(&totalSize), 4);
    out.write(reinterpret_cast<const char*>(&featureTableJSONSize), 4);
    out.write(reinterpret_cast<const char*>(&featureTableBinarySize), 4);
    out.write(reinterpret_cast<const char*>(&batchTableJSONSize), 4);
    out.write(reinterpret_cast<const char*>(&batchTableBinarySize), 4);

    out.write(ftJSON.data(), ftJSON.size());
    out.write(btJSON.data(), btJSON.size());
    out.write(reinterpret_cast<const char*>(glbData.data()), glbData.size());
    return true;
}

std::vector<uint8_t> TilesetGenerator::meshesToGLB(
    const std::vector<mesh::BuildingMesh>& meshes
) {
    tinygltf::Model model;
    model.scenes.resize(1);
    model.defaultScene = 0;

    model.buffers.resize(1);
    tinygltf::Buffer& bin = model.buffers[0];

    // Add a default material
    tinygltf::Material mat;
    mat.pbrMetallicRoughness.baseColorFactor = {0.8, 0.8, 0.8, 1.0};
    mat.pbrMetallicRoughness.metallicFactor = 0.0;
    mat.pbrMetallicRoughness.roughnessFactor = 1.0;
    mat.doubleSided = true;
    model.materials.push_back(mat);
    int defaultMat = 0;

    // Build geometry for each building mesh
    for (const auto& m : meshes) {
        std::vector<float> positions;
        std::vector<float> normals;
        std::vector<uint32_t> indices;

        // Compute face normals using PMP
        mesh::BuildingMesh& mutable_m = const_cast<mesh::BuildingMesh&>(m);
        pmp::face_normals(mutable_m);  // Compute face normals
        auto fnormals = m.get_face_property<pmp::Normal>("f:normal");

        if (!fnormals) {
            // If normals don't exist, compute them using PMP
            mesh::BuildingMesh& mutable_m = const_cast<mesh::BuildingMesh&>(m);
            pmp::vertex_normals(mutable_m);  // BuildingMesh inherits from SurfaceMesh
            fnormals = m.get_face_property<pmp::Normal>("f:normal");
        }

        // Use face normals - duplicate vertices for flat shading
        uint32_t nextIdx = 0;
        
        for (auto f : m.faces()) {
            // Get face vertices
            std::vector<pmp::Vertex> fverts;
            for (auto v : m.vertices(f)) {
                fverts.push_back(v);
            }
            
            if (fverts.size() < 3) continue;
            
            // Get the face normal
            pmp::Normal fn = fnormals[f];
            
            // Triangle fan for this face, each vertex gets the same face normal
            for (size_t k = 1; k + 1 < fverts.size(); ++k) {
                // Add 3 vertices for this triangle
                for (auto vi : {fverts[0], fverts[k], fverts[k+1]}) {
                    auto p = m.position(vi);
                    positions.push_back(static_cast<float>(p[0]));
                    positions.push_back(static_cast<float>(p[1]));
                    positions.push_back(static_cast<float>(p[2]));
                    
                    normals.push_back(static_cast<float>(fn[0]));
                    normals.push_back(static_cast<float>(fn[1]));
                    normals.push_back(static_cast<float>(fn[2]));
                    
                    indices.push_back(nextIdx++);
                }
            }
        }

        // Append positions to buffer
        const size_t posOffset = bin.data.size();
        bin.data.resize(bin.data.size() + positions.size() * sizeof(float));
        std::memcpy(bin.data.data() + posOffset, positions.data(), positions.size() * sizeof(float));
        alignTo4(bin.data);

        tinygltf::BufferView posBV;
        posBV.buffer = 0;
        posBV.byteOffset = static_cast<int>(posOffset);
        posBV.byteLength = static_cast<int>(positions.size() * sizeof(float));
        posBV.target = TINYGLTF_TARGET_ARRAY_BUFFER;
        const int posBVIndex = static_cast<int>(model.bufferViews.size());
        model.bufferViews.push_back(posBV);

        // Compute min/max for positions
        std::array<float,3> minv{std::numeric_limits<float>::max(),
                                 std::numeric_limits<float>::max(),
                                 std::numeric_limits<float>::max()};
        std::array<float,3> maxv{std::numeric_limits<float>::lowest(),
                                 std::numeric_limits<float>::lowest(),
                                 std::numeric_limits<float>::lowest()};
        for (size_t i = 0; i < positions.size(); i += 3) {
            minv[0] = std::min(minv[0], positions[i+0]);
            minv[1] = std::min(minv[1], positions[i+1]);
            minv[2] = std::min(minv[2], positions[i+2]);
            maxv[0] = std::max(maxv[0], positions[i+0]);
            maxv[1] = std::max(maxv[1], positions[i+1]);
            maxv[2] = std::max(maxv[2], positions[i+2]);
        }

        tinygltf::Accessor posAcc;
        posAcc.bufferView = posBVIndex;
        posAcc.byteOffset = 0;
        posAcc.componentType = TINYGLTF_COMPONENT_TYPE_FLOAT;
        posAcc.count = static_cast<int>(positions.size() / 3);
        posAcc.type = TINYGLTF_TYPE_VEC3;
        posAcc.minValues = { minv[0], minv[1], minv[2] };
        posAcc.maxValues = { maxv[0], maxv[1], maxv[2] };
        const int posAccIndex = static_cast<int>(model.accessors.size());
        model.accessors.push_back(posAcc);

        // Append normals to buffer
        const size_t normOffset = bin.data.size();
        bin.data.resize(bin.data.size() + normals.size() * sizeof(float));
        std::memcpy(bin.data.data() + normOffset, normals.data(), normals.size() * sizeof(float));
        alignTo4(bin.data);

        tinygltf::BufferView normBV;
        normBV.buffer = 0;
        normBV.byteOffset = static_cast<int>(normOffset);
        normBV.byteLength = static_cast<int>(normals.size() * sizeof(float));
        normBV.target = TINYGLTF_TARGET_ARRAY_BUFFER;
        const int normBVIndex = static_cast<int>(model.bufferViews.size());
        model.bufferViews.push_back(normBV);

        tinygltf::Accessor normAcc;
        normAcc.bufferView = normBVIndex;
        normAcc.byteOffset = 0;
        normAcc.componentType = TINYGLTF_COMPONENT_TYPE_FLOAT;
        normAcc.count = static_cast<int>(normals.size() / 3);
        normAcc.type = TINYGLTF_TYPE_VEC3;
        const int normAccIndex = static_cast<int>(model.accessors.size());
        model.accessors.push_back(normAcc);

        // Append indices
        const size_t idxOffset = bin.data.size();
        bin.data.resize(bin.data.size() + indices.size() * sizeof(uint32_t));
        std::memcpy(bin.data.data() + idxOffset, indices.data(), indices.size() * sizeof(uint32_t));
        alignTo4(bin.data);

        tinygltf::BufferView idxBV;
        idxBV.buffer = 0;
        idxBV.byteOffset = static_cast<int>(idxOffset);
        idxBV.byteLength = static_cast<int>(indices.size() * sizeof(uint32_t));
        idxBV.target = TINYGLTF_TARGET_ELEMENT_ARRAY_BUFFER;
        const int idxBVIndex = static_cast<int>(model.bufferViews.size());
        model.bufferViews.push_back(idxBV);

        tinygltf::Accessor idxAcc;
        idxAcc.bufferView = idxBVIndex;
        idxAcc.byteOffset = 0;
        idxAcc.componentType = TINYGLTF_COMPONENT_TYPE_UNSIGNED_INT;
        idxAcc.count = static_cast<int>(indices.size());
        idxAcc.type = TINYGLTF_TYPE_SCALAR;
        const int idxAccIndex = static_cast<int>(model.accessors.size());
        model.accessors.push_back(idxAcc);

        // Primitive
        tinygltf::Primitive prim;
        prim.indices = idxAccIndex;
        prim.mode = TINYGLTF_MODE_TRIANGLES;
        prim.material = defaultMat;
        prim.attributes["POSITION"] = posAccIndex;
        prim.attributes["NORMAL"] = normAccIndex;  // Add normals

        // Mesh
        tinygltf::Mesh gltfMesh;
        gltfMesh.primitives.push_back(prim);
        const int meshIndex = static_cast<int>(model.meshes.size());
        model.meshes.push_back(gltfMesh);

        // Node
        tinygltf::Node node;
        node.mesh = meshIndex;
        const int nodeIndex = static_cast<int>(model.nodes.size());
        model.nodes.push_back(node);

        model.scenes[0].nodes.push_back(nodeIndex);
    }

    model.buffers[0].uri.clear();

    // Write to temp file, read back
    tinygltf::TinyGLTF gltf;
    const fs::path tilesDir = fs::path(options_.outputDir) / "tiles";
    std::error_code ec;
    fs::create_directories(tilesDir, ec);

    const fs::path tmpGlb = tilesDir / "__aphelion_tmp.glb";

    const bool ok = gltf.WriteGltfSceneToFile(
        &model,
        tmpGlb.string(),
        true, true, false, true
    );

    if (!ok) {
        std::cerr << "tinygltf failed to write GLB\n";
        return {};
    }

    std::vector<uint8_t> data;
    std::ifstream in(tmpGlb, std::ios::binary | std::ios::ate);
    if (!in) return {};
    const std::streamsize size = in.tellg();
    in.seekg(0, std::ios::beg);
    data.resize(static_cast<size_t>(size));
    if (!in.read(reinterpret_cast<char*>(data.data()), size)) return {};
    in.close();

    fs::remove(tmpGlb, ec);

    return data;
}

} // namespace aphelion3d::tiles3d