#pragma once

#include "mbox/geometry/mesh.h"
#include "mbox/geometry/bvh.h"
#include "mbox/geometry/mesh_conversions.h"


namespace mbox {

typedef Mesh<
    make_vertex<double, 3, vertex_attributes::position>,
    unsigned short,
    3>
        CollisionMesh;

} // namespace mbox
