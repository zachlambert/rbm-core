#pragma once

#include "math/geometry/mesh.h"
#include "math/geometry/bvh.h"
#include "math/geometry/mesh_conversions.h"


namespace math {

typedef Mesh<
    make_vertex<double, 3, vertex_attributes::position>,
    unsigned short,
    3>
        CollisionMesh;

} // namespace math
