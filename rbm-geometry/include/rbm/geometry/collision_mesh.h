#pragma once

#include "rbm/geometry/mesh.h"
#include "rbm/geometry/bvh.h"
#include "rbm/geometry/mesh_conversions.h"


namespace rbm {

typedef Mesh<
    make_vertex<double, 3, vertex_attributes::position>,
    unsigned short,
    3>
        CollisionMesh;

} // namespace rbm
