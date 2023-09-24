#pragma once

#include "math/types/matrix.h"
#include "math/geometry/bounding.h"
#include <array>

namespace math {

template <typename IndexType, int VertexCount>
struct Facet {
    std::array<IndexType, VertexCount> indices;
};


} // namespace geometry
