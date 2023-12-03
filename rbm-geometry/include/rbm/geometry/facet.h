#pragma once

#include <array>
#include <rbm/types/matrix.h>
#include "rbm/geometry/bounding.h"


namespace rbm {

template <typename IndexType, int VertexCount>
struct Facet {
    std::array<IndexType, VertexCount> indices;
};


} // namespace rbm
