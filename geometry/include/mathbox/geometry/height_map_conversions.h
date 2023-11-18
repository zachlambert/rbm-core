#pragma once

#include "mathbox/geometry/height_map.h"
#include "mathbox/geometry/height_map_grid.h"
#include "mathbox/geometry/visual_mesh.h"


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
