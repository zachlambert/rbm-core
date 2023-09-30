#pragma once

#include "owl/geometry/mesh.h"
#include "owl/geometry/bvh.h"
#include "owl/geometry/mesh_conversions.h"


namespace owl {

typedef Mesh<
    make_vertex<double, 3, vertex_attributes::position>,
    unsigned short,
    3>
        CollisionMesh;

} // namespace owl
