#include "mesh.h"
#include "utilities.h"
#include <iostream>
#include <fstream>


namespace aphelion3d::output {
    // Output writer functionality can be added here in the future
    // Function to write objs as groups in a single file could be implemented here
     bool writeBuildingMeshesToOBJ(
        const std::vector<mesh::BuildingMesh>& meshes,
        const std::string& filename,
        bool useGroups = true
    );

    bool writeBuildingMeshesToMultiObjs(
        const std::vector<mesh::BuildingMesh>& meshes,
        const std::string& directory
    );
} // namespace aphelion3d::output