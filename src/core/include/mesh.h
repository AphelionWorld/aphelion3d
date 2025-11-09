#pragma once

#include <pmp/surface_mesh.h>
#include <pmp/algorithms/Normals.h>
#include <ogrsf_frmts.h> // GDAL for OGRPolygon
#include "utilities.h"

// namespace aphelion3d::mesh
namespace aphelion3d::mesh {
    // Mesh library functionality

// using PMP SurfaceMesh as the base mesh type
    using Mesh = pmp::SurfaceMesh;

    // We create a base mesh class that inherits from PMP SurfaceMesh
    class BaseMesh : public Mesh 
    {
    public:
        BaseMesh() = default;
        ~BaseMesh() override = default;

        // Additional mesh functionalities can be added here
        bool isValid() const 
        {
            return !this->isValid();
        }
    };

    //we create a building class that inherits from BaseMesh. The mesh will be
    // a building mesh created from a 2D footprint and extruded to a certain height
    class BuildingMesh : public BaseMesh 
    {
    public:
        BuildingMesh() = default;
        ~BuildingMesh() override = default;

        // Additional building mesh functionalities can be added here
        
        // Function to create a building mesh from a 2D footprint and height
        void createFromFootprintAndHeight(OGRPolygon* footprint, double height);

        /**
         * @brief Create a building mesh with a custom roof profile
         * @param footprint 2D footprint polygon
         * @param roofProfile 3D roof profile with Z coordinates
         * @param height Building height
         */
        void createFromFootprintAndRoofProfile(
            const utilities::Polygon2D& footprint, 
            const utilities::Polygon3D& roofProfile,
            double height
        );

    private:

        OGRPolygon* footprint; // 2D footprint of the building
        double height;        // Height of the building
    };
}
