#pragma once

#include "math/geometry/mesh.h"
#include "math/geometry/vertex.h"

namespace math {

using VisualMeshVertex = make_vertex<float, 3,
    vertex_attributes::position |
    vertex_attributes::normal |
    vertex_attributes::color
>;
using VisualMesh = Mesh<VisualMeshVertex, unsigned short, 3>;

} // namespace math
