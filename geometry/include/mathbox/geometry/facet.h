#pragma once

#include "mbox/types/matrix.h"
#include "mbox/geometry/bounding.h"
#include <array>

namespace mbox {

template <typename IndexType, int VertexCount>
struct Facet {
    std::array<IndexType, VertexCount> indices;
};


} // namespace mbox