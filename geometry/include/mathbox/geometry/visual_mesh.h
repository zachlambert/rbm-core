#pragma once

#include "owl/geometry/mesh.h"
#include "owl/geometry/vertex.h"

namespace owl {

using VisualMeshVertex = make_vertex<float, 3,
    vertex_attributes::position |
    vertex_attributes::normal |
    vertex_attributes::color
>;
using VisualMesh = Mesh<VisualMeshVertex, unsigned short, 3>;

} // namespace owl
