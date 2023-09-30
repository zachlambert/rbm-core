#pragma once

#include "owl/geometry/height_map.h"
#include "owl/geometry/height_map_grid.h"
#include "owl/geometry/visual_mesh.h"


namespace owl {

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

} // namespace owl
