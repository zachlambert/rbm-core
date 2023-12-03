
#include "rbm/geometry/height_map_conversions.h"
#include <iostream>


namespace rbm {

VisualMesh height_map_to_visual_mesh(
    const HeightMapInterface& height_map,
    const Vector2d& origin,
    const Vector2d& size,
    const Vector2i& divisions,
    const ColorRGBd& color) {

    VisualMesh mesh;
    VisualMesh::Vertex vertex;

    // Order of vertices in cell (cw = cell width, ch = cell height)
    // (relative to cell centre)
    // (-cw/2, -ch/2)
    // (+cw/2, -ch/2)
    // (-cw/2, +ch/2)
    // (+cw/2, +ch/2)

    Matrixd<2, 3> projection = Matrixd<2, 3>::Zero();
    projection.block<2, 2>(0, 0).setIdentity();

    Vector2d cell_size = size.array() / divisions.cast<double>().array();

    for (Vector2i index = Vector2i::Zero(); index != divisions; step_grid_index(index, divisions)) {
        Vector2d cell_centre;
        cell_centre = origin.array() + (index.cast<double>().array() + 0.5).array() * cell_size.array();

        Vector3d normal = height_map.normal(cell_centre);

        std::size_t vertices_offset = mesh.vertices.size();
        for (int k = 0; k < 4; k++) {
            Vector2i corner_index = index;
            if (k % 2 == 1) {
                corner_index(0)++;
            }
            if (k >= 2) {
                corner_index(1)++;
            }
            Vector2d corner_position =
                origin + Vector2d(corner_index.cast<double>().array() * cell_size.array());
            double height = height_map.height(corner_position);

            Vector3d position;
            position.head<2>() = corner_position;
            position.z() = height;

            vertex.position = position.cast<float>();
            vertex.normal = normal.cast<float>();
            vertex.color = color.cast<float>();
            mesh.vertices.push_back(vertex);
        }

        VisualMesh::Facet cell_triangle_1;
        cell_triangle_1.indices = { 0, 1, 2 };
        for (auto& index: cell_triangle_1.indices) index += vertices_offset;
        VisualMesh::Facet cell_triangle_2;
        cell_triangle_2.indices = { 3, 2, 1 };
        for (auto& index: cell_triangle_2.indices) index += vertices_offset;

        mesh.facets.push_back(cell_triangle_1);
        mesh.facets.push_back(cell_triangle_2);
    }
    return mesh;
}

HeightMapGrid height_map_to_height_map_grid(
    const HeightMapInterface& height_map,
    const Vector2d& origin,
    const Vector2d& size,
    const Vector2i& divisions)
{
    Vector2i divisions_plus_one = divisions + Vector2i::Ones();

    HeightMapGrid grid;
    grid.resize(divisions);
    grid.origin() = origin;
    grid.size() = size;

    Vector2d cell_size = size.array() / divisions.cast<double>().array();

    for (
        Vector2i index = Vector2i::Zero();
        index != divisions_plus_one;
        step_grid_index(index, divisions_plus_one))
    {
        Vector2d point = origin + Vector2d(index.cast<double>().array() * cell_size.array());
        grid.get_height(index) = height_map.height(point);
        std::cout << "Point: " << point.transpose() << std::endl;
        std::cout << "Height: " << grid.get_height(index) << std::endl;
    }

    Vector2i i00(0, 0);
    Vector2i i01(0, 1);
    Vector2i i10(1, 0);
    Vector2i i11(1, 1);

    for (
        Vector2i index = Vector2i::Zero();
        index != divisions;
        step_grid_index(index, divisions))
    {
        Vector2d p00 = origin + Vector2d((index + i00).cast<double>().array() * cell_size.array());
        Vector2d p01 = origin + Vector2d((index + i01).cast<double>().array() * cell_size.array());
        Vector2d p10 = origin + Vector2d((index + i10).cast<double>().array() * cell_size.array());
        Vector2d p11 = origin + Vector2d((index + i11).cast<double>().array() * cell_size.array());

        // Displacement from bottom left to top right
        Vector3d a;
        a.head<2>() = p11 - p00;
        a.z() = height_map.height(p11) - height_map.height(p00);

        // Displacement from top left to bottom right
        Vector3d b;
        b.head<2>() = p10 - p01;
        b.z() = height_map.height(p10) - height_map.height(p01);

        grid.get_normal(index) = (b.cross(a)).normalized();
    }

    return grid;
}

} // namespace rbm
