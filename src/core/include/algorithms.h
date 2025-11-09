#include "utilities.h"
#include "mesh.h"


namespace aphelion3d::algorithms {

    void generateBuildingMeshesFromFile(std::string filepath, std::vector<mesh::BuildingMesh>& buildingMeshes) ;

    void generateBuildingMesh(OGRPolygon* footprint, double height);

} // namespace aphelion3d::algorithms