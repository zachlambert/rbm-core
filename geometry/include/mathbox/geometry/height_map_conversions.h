#pragma once

#include "mbox/geometry/height_map.h"
#include "mbox/geometry/height_map_grid.h"
#include "mbox/geometry/visual_mesh.h"


namespace mbox {

VisualMesh height_map_to_visual_mesh(
    const HeightMapInterface& height_map,
    const Vector2d& origin,
    const Vector2d& size,
    const Vector2i& divisions,
    const ColorRGBd& color);

HeightMapGrid height_map_to_height_map_grid(
    const HeightMapInterface& height_map,
    const Vector2d& origin,
    const Vector2d& size,
    const Vector2i& divisions);

} // namespace mbox
