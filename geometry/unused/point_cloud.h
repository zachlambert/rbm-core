#pragma once

#include <vector>
#include "math/geometry/vertex.h"

namespace math {

template <typename Vertex>
struct PointCloud {
    std::vector<Vertex> vertices;
    size_t width;
    size_t height; // Set to 1 for unstructured point cloud
};

} // namespace math
