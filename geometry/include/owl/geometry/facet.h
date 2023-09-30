#pragma once

#include "owl/types/matrix.h"
#include "owl/geometry/bounding.h"
#include <array>

namespace owl {

template <typename IndexType, int VertexCount>
struct Facet {
    std::array<IndexType, VertexCount> indices;
};


} // namespace geometry
