#pragma once

#include "rbm/geometry/mesh.h"
#include "rbm/geometry/vertex.h"

namespace rbm {

using VisualMeshVertex = make_vertex<float, 3,
    vertex_attributes::position |
    vertex_attributes::normal |
    vertex_attributes::color
>;
using VisualMesh = Mesh<VisualMeshVertex, unsigned short, 3>;

} // namespace rbm
