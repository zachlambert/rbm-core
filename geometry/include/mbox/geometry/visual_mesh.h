#pragma once

#include "mbox/geometry/mesh.h"
#include "mbox/geometry/vertex.h"

namespace mbox {

using VisualMeshVertex = make_vertex<float, 3,
    vertex_attributes::position |
    vertex_attributes::normal |
    vertex_attributes::color
>;
using VisualMesh = Mesh<VisualMeshVertex, unsigned short, 3>;

} // namespace mbox
