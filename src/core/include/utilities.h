#pragma once

#include <pmp/surface_mesh.h>   // Public API exposes pmp::SurfaceMesh and pmp::Vertex
#include <ogrsf_frmts.h>        // Public API uses OGRPolygon/OGRLinearRing pointers
#include <vector>
#include <array>
#include <string>

namespace aphelion3d::utilities {

    // Basic geometry typedefs
    using Point2D   = std::array<double, 2>;
    using Point3D   = std::array<double, 3>;
    using Ring2D    = std::vector<Point2D>;
    using Ring3D    = std::vector<Point3D>;
    using Polygon2D = std::vector<Ring2D>;
    using Polygon3D = std::vector<Ring3D>;

    // Container for footprint + optional 3D roof information
    struct PolygonData {
        Polygon2D footprint2D;            // 2D base footprint
        Polygon3D roofProfile;            // 3D roof profile (exterior first, then holes), if available
        double    height = 0.0;           // Height (attribute for 2D; derived for 2.5D/3D)
        bool      hasRoofProfile = false; // True if geometry provided Zs for a roof profile
    };

    // Footprint conversion (2D)
    Polygon2D polygonToPointArrays(OGRPolygon* polygon);
    Ring2D    extractExteriorRing(OGRLinearRing* ring);
    std::vector<Ring2D> extractInteriorRings(OGRPolygon* polygon);

    // Wall-face creation between two rings
    void createWallFaces(
        pmp::SurfaceMesh& mesh,
        const std::vector<pmp::Vertex>& baseVertices,
        const std::vector<pmp::Vertex>& topVertices,
        bool reverse = false
    );

    // Validation and metrics
    bool   validatePolygon(OGRPolygon* polygon);
    double calculatePolygonArea(OGRPolygon* polygon);

    // Readers
    std::vector<PolygonData> readOGRPolygonsWithProfile(std::string filepath);

    // 3D variants (preserve Z)
    Ring3D extractExteriorRing3D(OGRLinearRing* ring);
    std::vector<Ring3D> extractInteriorRings3D(OGRPolygon* polygon);

    // Removes vertices in a ring that are closer than tolerance to the previously kept vertex.
    // Returns the number of removed distinct vertices (excludes the implicit closing point).
    // Does not modify the ring if the result would have fewer than 3 distinct vertices.
    std::size_t removeClosePointsFromRing(OGRLinearRing* ring, double tolerance);

    // Applies removeClosePointsFromRing to exterior and interior rings, then re-enforces winding.
    // Returns the total count of removed distinct vertices across all rings.
    std::size_t removeClosePoints(OGRPolygon* polygon, double tolerance);

    void ensureCorrectWinding(OGRPolygon* polygon);

    void transformRingToECEF(Ring2D& ring, OGRCoordinateTransformation* ct);
} // namespace aphelion3d::utilities