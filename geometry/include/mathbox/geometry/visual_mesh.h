#pragma once

#include "mathbox/geometry/mesh.h"
#include "mathbox/geometry/vertex.h"

namespace mbox {

using VisualMeshVertex = make_vertex<float, 3,
    vertex_attributes::position |
    vertex_attributes::normal |
    vertex_attributes::color
>;
using VisualMesh = Mesh<VisualMeshVertex, unsigned short, 3>;

} // namespace mbox
