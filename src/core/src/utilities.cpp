#include "utilities.h"
#include <cmath>

//include earcut header
#include <poly2tri/poly2tri.h>

namespace aphelion3d::utilities {

// Small helper: squared 2D distance between points
static inline double dist2XY(const OGRPoint& a, const OGRPoint& b) {
    const double dx = a.getX() - b.getX();
    const double dy = a.getY() - b.getY();
    return dx * dx + dy * dy;
}

std::size_t removeClosePointsFromRing(OGRLinearRing* ring, double tolerance) {
    if (!ring) return 0;

    const int n = ring->getNumPoints();
    // Need at least 4 points (3 distinct + closing point) to be a valid ring.
    if (n <= 4) return 0;

    const bool hasZ = ring->getCoordinateDimension() == 3;
    const double tol2 = tolerance * tolerance;

    // Gather original points (exclude the closing point at the end)
    std::vector<OGRPoint> pts;
    pts.reserve(std::max(0, n - 1));
    {
        OGRPoint p;
        for (int i = 0; i < n - 1; ++i) {
            ring->getPoint(i, &p);
            pts.push_back(p);
        }
    }

    // Build deduplicated sequence
    std::vector<OGRPoint> out;
    out.reserve(pts.size());
    if (!pts.empty()) {
        out.push_back(pts.front());
        for (std::size_t i = 1; i < pts.size(); ++i) {
            if (dist2XY(pts[i], out.back()) >= tol2) {
                out.push_back(pts[i]);
            }
        }
        // Also ensure last kept is not too close to the first kept (before re-closing)
        if (out.size() >= 2 && dist2XY(out.back(), out.front()) < tol2) {
            out.pop_back();
        }
    }

    // If fewer than 3 distinct vertices, don't modify (would produce invalid ring)
    if (out.size() < 3) {
        return 0;
    }

    const std::size_t removed = (pts.size() > out.size()) ? (pts.size() - out.size()) : 0;

    // Rebuild the ring: clear, add kept points, and re-add the closing point (same as first)
    ring->empty();

    if (hasZ) {
        for (const auto& p : out) {
            ring->addPoint(p.getX(), p.getY(), p.getZ());
        }
        // Closing point equals first
        ring->addPoint(out.front().getX(), out.front().getY(), out.front().getZ());
    } else {
        for (const auto& p : out) {
            ring->addPoint(p.getX(), p.getY());
        }
        // Closing point equals first
        ring->addPoint(out.front().getX(), out.front().getY());
    }

    return removed;
}

std::size_t removeClosePoints(OGRPolygon* polygon, double tolerance) {
    if (!polygon) return 0;

    std::size_t totalRemoved = 0;

    // Exterior
    if (auto* exterior = polygon->getExteriorRing()) {
        totalRemoved += removeClosePointsFromRing(exterior, tolerance);
    }

    // Interiors
    const int holes = polygon->getNumInteriorRings();
    for (int i = 0; i < holes; ++i) {
        if (auto* interior = polygon->getInteriorRing(i)) {
            totalRemoved += removeClosePointsFromRing(interior, tolerance);
        }
    }

    // Ensure we still meet winding conventions after edits
    ensureCorrectWinding(polygon);

    return totalRemoved;
}

    Polygon2D polygonToPointArrays(OGRPolygon* polygon) 
    {
        Polygon2D result;
        
        if (!polygon) {
            return result;
        }

        // Add exterior ring
        OGRLinearRing* exteriorRing = polygon->getExteriorRing();
        if (exteriorRing) {
            result.push_back(extractExteriorRing(exteriorRing));
        }

        // Add interior rings
        auto holes = extractInteriorRings(polygon);
        result.insert(result.end(), holes.begin(), holes.end());

        return result;
    }

    Ring2D extractExteriorRing(OGRLinearRing* ring) 
    {
        Ring2D points;

        if (!ring) {
            return points;
        }

        int numPoints = ring->getNumPoints();
        // Skip last point as OGR rings are closed (first == last)
        for (int i = 0; i < numPoints - 1; ++i) {
            points.push_back({ring->getX(i), ring->getY(i)});
        }

        return points;
    }
    
    std::vector<Ring2D> extractInteriorRings(OGRPolygon* polygon) {
        std::vector<Ring2D> holes;

        if (!polygon) {
            return holes;
        }

        int numInteriorRings = polygon->getNumInteriorRings();
        for (int i = 0; i < numInteriorRings; ++i) {
            OGRLinearRing* interiorRing = polygon->getInteriorRing(i);
            if (interiorRing) {
                holes.push_back(extractExteriorRing(interiorRing));
            }
        }

        return holes;
    }

    void createWallFaces(
        pmp::SurfaceMesh& mesh,
        const std::vector<pmp::Vertex>& baseVertices,
        const std::vector<pmp::Vertex>& topVertices,
        bool reverse
    ) {
        int n = baseVertices.size();
        
        if (n != topVertices.size() || n < 2) {
            return;
        }

        for (int i = 0; i < n; ++i) {
            int next = (i + 1) % n;
            
            if (reverse) {
                // Reversed winding for interior holes
                //std::cout << "Creating wall faces for interior hole with " << baseVertices.size() << " vertices.\n";    
                //std::cout << "Adding wall face triangles between base vertex " << i << " and " << next << "\n";
                mesh.add_triangle(baseVertices[i], topVertices[i], baseVertices[next]);
                mesh.add_triangle(baseVertices[next], topVertices[i], topVertices[next]);
            } else {
                //std::cout << "Creating wall faces for exterior ring with " << baseVertices.size() << " vertices.\n";    
                //std::cout << "Adding wall face triangles between base vertex " << i << " and " << next << "\n";
                // Normal winding for exterior
                mesh.add_triangle(baseVertices[i], baseVertices[next], topVertices[i]);
                mesh.add_triangle(baseVertices[next], topVertices[next], topVertices[i]);
            }
        }
    }

    bool validatePolygon(OGRPolygon* polygon) {
        if (!polygon) {
            return false;
        }

        OGRLinearRing* exteriorRing = polygon->getExteriorRing();
        if (!exteriorRing) {
            return false;
        }

        if (exteriorRing->getNumPoints() < 4) { // 3 points + closing point
            return false;
        }

        return true;
    }

    double calculatePolygonArea(OGRPolygon* polygon) 
    {
        if (!polygon) {
            return 0.0;
        }

        return polygon->get_Area();
    }

    
    Ring3D extractExteriorRing3D(OGRLinearRing* ring) 
    {
        Ring3D points;

        if (!ring) {
            return points;
        }

        int numPoints = ring->getNumPoints();
        bool hasZ = ring->getCoordinateDimension() == 3;
        
        // Skip last point as OGR rings are closed (first == last)
        for (int i = 0; i < numPoints - 1; ++i) {
            if (hasZ) {
                points.push_back({ring->getX(i), ring->getY(i), ring->getZ(i)});
            } else {
                points.push_back({ring->getX(i), ring->getY(i), 0.0});
            }
        }

        return points;
    }

    std::vector<Ring3D> extractInteriorRings3D(OGRPolygon* polygon) 
    {
        std::vector<Ring3D> holes;

        if (!polygon) {
            return holes;
        }

        int numInteriorRings = polygon->getNumInteriorRings();
        for (int i = 0; i < numInteriorRings; ++i) {
            OGRLinearRing* interiorRing = polygon->getInteriorRing(i);
            if (interiorRing) {
                holes.push_back(extractExteriorRing3D(interiorRing));
            }
        }

        return holes;
    }

    void ensureCorrectWinding(OGRPolygon* polygon) 
    {
        if (!polygon) return;

        // Fix exterior ring (should be counter-clockwise)
        OGRLinearRing* exteriorRing = polygon->getExteriorRing();
        if (exteriorRing) {
            if (!exteriorRing->isClockwise()) {
                // OGR's isClockwise returns true if CW, false if CCW
                // We want CCW for exterior, so if it's already CCW (false), do nothing
                // If it returns true (CW), we need to reverse it
            } else {
                // It's clockwise, reverse it to make it CCW
                exteriorRing->reversePoints();
            }
        }

        // Fix interior rings (holes should be clockwise)
        int numInteriorRings = polygon->getNumInteriorRings();
        for (int i = 0; i < numInteriorRings; ++i) {
            OGRLinearRing* interiorRing = polygon->getInteriorRing(i);
            if (interiorRing) {
                if (interiorRing->isClockwise()) {
                    // Already clockwise, keep it
                } else {
                    // It's CCW, reverse it to make it CW
                    interiorRing->reversePoints();
                }
            }
        }
    }

    std::vector<PolygonData> readOGRPolygonsWithProfile(std::string filepath) 
    {
        std::vector<PolygonData> polygons;

        GDALAllRegister();
        GDALDataset* dataset = (GDALDataset*) GDALOpenEx(filepath.c_str(), GDAL_OF_VECTOR, nullptr, nullptr, nullptr);
        if (!dataset) {
            return polygons;
        }

        OGRLayer* layer = dataset->GetLayer(0);
        if (!layer) {
            GDALClose(dataset);
            return polygons;
        }

        layer->ResetReading();
        OGRFeature* feature;
        while ((feature = layer->GetNextFeature()) != nullptr) 
        {
            OGRGeometry* geometry = feature->GetGeometryRef();
            if (geometry != nullptr && wkbFlatten(geometry->getGeometryType()) == wkbPolygon) 
            {
                OGRPolygon* ogrPolygon = geometry->toPolygon();
                removeClosePoints(ogrPolygon, 0.005);
                ensureCorrectWinding(ogrPolygon);
                PolygonData polyData;   
                
                // Check if geometry has Z coordinates
                bool is3D = wkbHasZ(geometry->getGeometryType());
                polyData.hasRoofProfile = is3D;
                
                // Always extract 2D footprint
                polyData.footprint2D = polygonToPointArrays(ogrPolygon);
                
                if (is3D) 
                {
                    // Extract 3D roof profile
                    OGRLinearRing* exteriorRing = ogrPolygon->getExteriorRing();
                    if (exteriorRing) {
                        polyData.roofProfile.push_back(extractExteriorRing3D(exteriorRing));
                    }
                    
                    // Add interior rings in 3D
                    auto holes3D = extractInteriorRings3D(ogrPolygon);
                    polyData.roofProfile.insert(polyData.roofProfile.end(), holes3D.begin(), holes3D.end());
                    
                    // Calculate height from Z variation
                    if (exteriorRing && exteriorRing->getNumPoints() > 0) {
                        double minZ = exteriorRing->getZ(0);
                        double maxZ = exteriorRing->getZ(0);
                        
                        for (int i = 1; i < exteriorRing->getNumPoints(); ++i) {
                            double z = exteriorRing->getZ(i);
                            if (z < minZ) minZ = z;
                            if (z > maxZ) maxZ = z;
                        }
                        
                        polyData.height = maxZ - minZ;
                    }
                } 
                else 
                {
                    // Read height from attribute
                    int heightFieldIndex = feature->GetFieldIndex("HEIGHT");
                    //std::cout << "Reading height from attribute: " << heightFieldIndex << "\n";
                    if (heightFieldIndex != -1) {
                        polyData.height = feature->GetFieldAsDouble(heightFieldIndex);
                    } else {
                        heightFieldIndex = feature->GetFieldIndex("height");
                        if (heightFieldIndex != -1) {
                            polyData.height = feature->GetFieldAsDouble(heightFieldIndex);
                        } else {
                            heightFieldIndex = feature->GetFieldIndex("HAUTEUR");
                            if (heightFieldIndex != -1) {
                                polyData.height = feature->GetFieldAsDouble(heightFieldIndex);
                            } else {
                                polyData.height = 10.0; // Default height
                            }
                        }
                    }
                }

                polygons.push_back(polyData);
            }
            OGRFeature::DestroyFeature(feature);
        }

        GDALClose(dataset);
        return polygons;
    }

// Transform all ring vertices to ECEF in-place
void transformRingToECEF(Ring2D& ring, OGRCoordinateTransformation* ct) {
    for (auto& pt : ring) {
        double x = pt[0], y = pt[1], z = 0.0;
        if (ct->Transform(1, &x, &y, &z)) {
            pt[0] = x;
            pt[1] = y;
            // If you store 3D, update z as well
        }
    }
}

    
} // namespace aphelion3d::utilities