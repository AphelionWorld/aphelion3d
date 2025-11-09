#include "mesh.h"
#include "utilities.h"
#include <poly2tri/poly2tri.h>
#include <iostream> 

namespace aphelion3d::mesh {
    
    void BuildingMesh::createFromFootprintAndHeight(OGRPolygon* footprint, double height) {
        if (!utilities::validatePolygon(footprint) || height <= 0.0001) {
            return;
        }

        this->footprint = footprint;
        this->height = height;
        this->clear();

        // Get all rings (exterior + interior holes)
        auto polygon = utilities::polygonToPointArrays(footprint);
        
        // Prepare poly2tri data structures
        std::vector<p2t::Point*> exteriorPoints;
        std::vector<std::vector<p2t::Point*>> holePoints;
        
        // Add exterior ring (first ring)
        if (!polygon.empty()) {
            for (const auto& pt : polygon[0]) {
                exteriorPoints.push_back(new p2t::Point(pt[0], pt[1]));
            }
        }
        
        // Create CDT with exterior ring
        p2t::CDT cdt(exteriorPoints);
        
        // Add holes (interior rings)
        for (size_t i = 1; i < polygon.size(); ++i) {
            std::vector<p2t::Point*> hole;
            for (const auto& pt : polygon[i]) {
                hole.push_back(new p2t::Point(pt[0], pt[1]));
            }
            cdt.AddHole(hole);
            holePoints.push_back(hole);
        }
        
        // Triangulate
        cdt.Triangulate();
        
        // Get triangles
        auto triangles = cdt.GetTriangles();

        // Create all vertices (for all rings combined)
        std::vector<pmp::Vertex> baseVertices;
        std::vector<pmp::Vertex> topVertices;
        
        for (const auto& ring : polygon) {
            for (const auto& point : ring) {
                baseVertices.push_back(this->add_vertex(pmp::Point(point[0], point[1], 0.0)));
                topVertices.push_back(this->add_vertex(pmp::Point(point[0], point[1], height)));
            }
        }
        
        // Build vertex map for triangle indices
        std::map<p2t::Point*, size_t> pointToIndex;
        size_t idx = 0;
        for (const auto& ring : polygon) {
            for (size_t i = 0; i < ring.size(); ++i) {
                // Find matching point in poly2tri structures
                for (auto* pt : exteriorPoints) {
                    if (std::abs(pt->x - ring[i][0]) < 1e-9 && std::abs(pt->y - ring[i][1]) < 1e-9) {
                        pointToIndex[pt] = idx;
                        break;
                    }
                }
                for (const auto& hole : holePoints) {
                    for (auto* pt : hole) {
                        if (std::abs(pt->x - ring[i][0]) < 1e-9 && std::abs(pt->y - ring[i][1]) < 1e-9) {
                            pointToIndex[pt] = idx;
                            break;
                        }
                    }
                }
                ++idx;
            }
        }

        // Create base faces using triangulation
        for (const auto tri : triangles) {
            size_t i0 = pointToIndex[tri->GetPoint(0)];
            size_t i1 = pointToIndex[tri->GetPoint(1)];
            size_t i2 = pointToIndex[tri->GetPoint(2)];
            
            this->add_triangle(
                baseVertices[i0], 
                baseVertices[i2], 
                baseVertices[i1]
            );
        }

        // Create top faces using triangulation (reversed winding)
        for (const auto tri : triangles) {
            size_t i0 = pointToIndex[tri->GetPoint(0)];
            size_t i1 = pointToIndex[tri->GetPoint(1)];
            size_t i2 = pointToIndex[tri->GetPoint(2)];
            
            this->add_triangle(
                topVertices[i0], 
                topVertices[i1], 
                topVertices[i2]
            );
        }

        // Create wall faces for exterior ring (polygon[0])
        int offset = 0;
        std::vector<pmp::Vertex> ringBase(baseVertices.begin() + offset, 
                                          baseVertices.begin() + offset + polygon[0].size());
        std::vector<pmp::Vertex> ringTop(topVertices.begin() + offset, 
                                         topVertices.begin() + offset + polygon[0].size());

        utilities::createWallFaces(*this, ringBase, ringTop, false);
        offset += polygon[0].size();

        // Create wall faces for interior rings (holes)
        for (size_t ringIdx = 1; ringIdx < polygon.size(); ++ringIdx) 
        {
            std::vector<pmp::Vertex> holeBase(baseVertices.begin() + offset, 
                                              baseVertices.begin() + offset + polygon[ringIdx].size());
            std::vector<pmp::Vertex> holeTop(topVertices.begin() + offset, 
                                             topVertices.begin() + offset + polygon[ringIdx].size());
            utilities::createWallFaces(*this, holeBase, holeTop, true);
            offset += polygon[ringIdx].size();
        }
        
        // Clean up poly2tri memory
        for (auto* pt : exteriorPoints) delete pt;
        for (auto& hole : holePoints) {
            for (auto* pt : hole) delete pt;
        }
    }

    void BuildingMesh::createFromFootprintAndRoofProfile(
        const utilities::Polygon2D& footprint, 
        const utilities::Polygon3D& roofProfile,
        double height
    ) {
        if (footprint.empty() || roofProfile.empty() || height <= 0.0) {
            return;
        }

        this->height = height;
        this->clear();

        // Prepare poly2tri data structures
        std::vector<p2t::Point*> exteriorPoints;
        std::vector<std::vector<p2t::Point*>> holePoints;
        
        // Add exterior ring
        if (!footprint.empty()) {
            for (const auto& pt : footprint[0]) {
                exteriorPoints.push_back(new p2t::Point(pt[0], pt[1]));
            }
        }
        
        p2t::CDT cdt(exteriorPoints);
        
        // Add holes
        for (size_t i = 1; i < footprint.size(); ++i) {
            std::vector<p2t::Point*> hole;
            for (const auto& pt : footprint[i]) {
                hole.push_back(new p2t::Point(pt[0], pt[1]));
            }
            cdt.AddHole(hole);
            holePoints.push_back(hole);
        }
        
        cdt.Triangulate();
        auto triangles = cdt.GetTriangles();

        // Create base vertices from 2D footprint
        std::vector<pmp::Vertex> baseVertices;
        for (const auto& ring : footprint) {
            for (const auto& point : ring) {
                baseVertices.push_back(this->add_vertex(pmp::Point(point[0], point[1], 0.0)));
            }
        }

        // Create roof vertices from 3D profile
        std::vector<pmp::Vertex> roofVertices;
        for (const auto& ring : roofProfile) {
            for (const auto& point : ring) {
                roofVertices.push_back(this->add_vertex(pmp::Point(point[0], point[1], point[2])));
            }
        }
        
        // Build vertex map
        std::map<p2t::Point*, size_t> pointToIndex;
        size_t idx = 0;
        for (const auto& ring : footprint) {
            for (size_t i = 0; i < ring.size(); ++i) {
                for (auto* pt : exteriorPoints) {
                    if (std::abs(pt->x - ring[i][0]) < 1e-9 && std::abs(pt->y - ring[i][1]) < 1e-9) {
                        pointToIndex[pt] = idx;
                        break;
                    }
                }
                for (const auto& hole : holePoints) {
                    for (auto* pt : hole) {
                        if (std::abs(pt->x - ring[i][0]) < 1e-9 && std::abs(pt->y - ring[i][1]) < 1e-9) {
                            pointToIndex[pt] = idx;
                            break;
                        }
                    }
                }
                ++idx;
            }
        }

        // Create base faces
        for (const auto tri : triangles) {
            size_t i0 = pointToIndex[tri->GetPoint(0)];
            size_t i1 = pointToIndex[tri->GetPoint(1)];
            size_t i2 = pointToIndex[tri->GetPoint(2)];
            
            this->add_triangle(
                baseVertices[i0], 
                baseVertices[i2], 
                baseVertices[i1]
            );
        }

        // Create roof faces (with potentially different Z values)
        for (const auto tri : triangles) {
            size_t i0 = pointToIndex[tri->GetPoint(0)];
            size_t i1 = pointToIndex[tri->GetPoint(1)];
            size_t i2 = pointToIndex[tri->GetPoint(2)];
            
            this->add_triangle(
                roofVertices[i0], 
                roofVertices[i1], 
                roofVertices[i2]
            );
        }

        // Create walls between base and roof
        int offset = 0;
        for (size_t ringIdx = 0; ringIdx < footprint.size(); ++ringIdx) {
            int ringSize = footprint[ringIdx].size();
            std::vector<pmp::Vertex> ringBase(baseVertices.begin() + offset, 
                                              baseVertices.begin() + offset + ringSize);
            std::vector<pmp::Vertex> ringRoof(roofVertices.begin() + offset, 
                                              roofVertices.begin() + offset + ringSize);
            
            bool isHole = (ringIdx > 0);
            utilities::createWallFaces(*this, ringBase, ringRoof, isHole);
            offset += ringSize;
        }
        
        // Clean up
        for (auto* pt : exteriorPoints) delete pt;
        for (auto& hole : holePoints) {
            for (auto* pt : hole) delete pt;
        }
    }
}