#pragma once

#include "math/geometry/height_map.h"
#include "math/geometry/height_map_grid.h"
#include "math/geometry/visual_mesh.h"


namespace math {

VisualMesh height_map_to_visual_mesh(
    const HeightMapInterface& height_map,
    const math::Vector2d& origin,
    const math::Vector2d& size,
    const math::Vector2i& divisions,
    const ColorRGBd& color);

HeightMapGrid height_map_to_height_map_grid(
    const HeightMapInterface& height_map,
    const math::Vector2d& origin,
    const math::Vector2d& size,
    const math::Vector2i& divisions);

} // namespace math
