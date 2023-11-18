#pragma once

#include "mathbox/types/matrix.h"
#include "mathbox/geometry/bounding.h"
#include <array>

namespace mbox {

template <typename IndexType, int VertexCount>
struct Facet {
    std::array<IndexType, VertexCount> indices;
};


} // namespace mbox