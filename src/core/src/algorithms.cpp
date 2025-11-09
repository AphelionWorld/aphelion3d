#include "algorithms.h"


namespace aphelion3d::algorithms {

    void generateBuildingMeshesFromFile(std::string filepath, std::vector<mesh::BuildingMesh>& buildingMeshes) 
    {
        auto polygonsData = utilities::readOGRPolygonsWithProfile(filepath);
        
        //std::cout << "Read " << polygonsData.size() << " polygons from file.\n";    
        for (const auto& polyData : polygonsData) 
        {
            mesh::BuildingMesh buildingMesh;
            
            if (polyData.hasRoofProfile) {
                // Use 3D profile for slanted/custom roofs
                buildingMesh.createFromFootprintAndRoofProfile(
                    polyData.footprint2D,
                    polyData.roofProfile,
                    polyData.height
                );
            } 
            else 
            {
                //std::cout << "Polygon has no roof profile, creating flat roof.\n";
                // Create flat roof from 2D footprint
                // Convert Polygon2D to OGRPolygon for the existing method
                OGRPolygon ogrPolygon;
                OGRLinearRing exteriorRing;
                for (const auto& point : polyData.footprint2D[0]) {
                    exteriorRing.addPoint(point[0], point[1]);
                }
                exteriorRing.closeRings();
                ogrPolygon.addRing(&exteriorRing);

                // Add interior rings (holes)
                for (size_t i = 1; i < polyData.footprint2D.size(); ++i) {
                    OGRLinearRing interiorRing;
                    for (const auto& point : polyData.footprint2D[i]) {
                        interiorRing.addPoint(point[0], point[1]);
                    }
                    interiorRing.closeRings();
                    ogrPolygon.addRing(&interiorRing);
                }

                buildingMesh.createFromFootprintAndHeight(&ogrPolygon, polyData.height);
            }

            //if(buildingMesh.isValid()) 
           // {
            buildingMeshes.push_back(buildingMesh);
          //  }
        }
        std::cout << "Generated " << buildingMeshes.size() << " building meshes.\n";
    }

    void generateBuildingMesh(OGRPolygon* footprint, double height) {
        if (!utilities::validatePolygon(footprint) || height <= 0.0) {
            return;
        }

        mesh::BuildingMesh buildingMesh;
        buildingMesh.createFromFootprintAndHeight(footprint, height);
    }

} // namespace aphelion3d::algorithms