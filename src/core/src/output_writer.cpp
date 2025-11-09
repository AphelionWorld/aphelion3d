#include "output_writer.h"

namespace aphelion3d::output {
    // Output writer functionality can be added here in the future
    // Function to write objs as groups in a single file could be implemented here
     bool writeBuildingMeshesToOBJ(
        const std::vector<mesh::BuildingMesh>& meshes,
        const std::string& filename,
        bool useGroups
    )
    {
        // This function will write multiple building meshes into a single OBJ file.
        // If useGroups is true, each building mesh will be written as a separate group.
        std::ofstream file(filename);
        if (!file.is_open()) return false;

        file << std::setprecision(12);
        file << "# Aphelion3D Building Meshes\n";
        file << "# Total buildings: " << meshes.size() << "\n\n";

        int vertexOffset = 0;

        for (size_t i = 0; i < meshes.size(); ++i) {
            const auto& mesh = meshes[i];
            
            if (useGroups) {
                file << "g building_" << i << "\n";
                file << "# Vertices: " << mesh.n_vertices() 
                    << " Faces: " << mesh.n_faces() << "\n";
            }
            
            // Vertices
            for (auto v : mesh.vertices()) {
                auto pos = mesh.position(v);
                file << "v " << pos[0] << " " << pos[1] << " " << pos[2] << "\n";
            }
            
            // Faces
            for (auto f : mesh.faces()) {
                file << "f";
                for (auto v : mesh.vertices(f)) {
                    file << " " << (v.idx() + vertexOffset + 1);
                }
                file << "\n";
            }
            
            if (useGroups) file << "\n";
            vertexOffset += mesh.n_vertices();
        }

        file.close();
        return true;
    }

    bool writeBuildingMeshesToMultiObjs(
        const std::vector<mesh::BuildingMesh>& meshes,
        const std::string& directory
    ) 
    {
        //create directory if it doesn't exist

        if(!std::filesystem::exists(directory))
            std::filesystem::create_directories(directory);

        // This function will write each building mesh into a separate OBJ file
        // in the specified directory.
        for (size_t i = 0; i < meshes.size(); ++i) 
        {
            const auto& mesh = meshes[i];
            std::string filename = directory + "/building_" + std::to_string(i) + ".obj";
            std::ofstream file(filename);
            if (!file.is_open()) return false;

            file << "# Aphelion3D Building Mesh\n";
            file << "# Vertices: " << mesh.n_vertices() 
                 << " Faces: " << mesh.n_faces() << "\n\n";

            // Vertices
            for (auto v : mesh.vertices()) {
                auto pos = mesh.position(v);
                file << "v " << pos[0] << " " << pos[1] << " " << pos[2] << "\n";
            }

            // Faces
            for (auto f : mesh.faces()) {
                file << "f";
                for (auto v : mesh.vertices(f)) {
                    file << " " << (v.idx() + 1);
                }
                file << "\n";
            }

            file.close();
        }

        return true;
    }   
} // namespace aphelion3d::output