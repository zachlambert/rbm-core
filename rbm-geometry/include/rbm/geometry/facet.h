#pragma once

#include "rbm/types/matrix.h"
#include "rbm/geometry/bounding.h"
#include <array>

namespace rbm {

template <typename IndexType, int VertexCount>
struct Facet {
    std::array<IndexType, VertexCount> indices;
};


} // namespace rbm