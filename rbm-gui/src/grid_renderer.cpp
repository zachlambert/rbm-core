#include "rbm/gui/grid_renderer.h"
#include "rbm/gui/shader.h"
#include "rbm/gui/shader/grid_fs.h"
#include "rbm/gui/shader/grid_vs.h"
#include <imgui.h>


namespace rbm {

GridRenderer::GridRenderer()
{
    program_id = load_shader(
        shader::grid_vs,
        shader::grid_fs
    );
    mvp_loc = glGetUniformLocation(program_id, "MVP");
}

int GridRenderer::create_grid()
{
    grids.emplace_back();
    return grids.size() - 1;
}

void GridRenderer::configure_grid(int grid_index, double width, double line_width, int num_major_divisions, int num_minor_divisions)
{
    if (grid_index < 0 || grid_index >= grids.size()) {
        return;
    }

    auto& grid = grids[grid_index];

    if (!grid.initialised) {
        glGenVertexArrays(1, &grid.static_VAO);
        glGenBuffers(1, &grid.static_VBO);
        grid.initialised = true;
    }

    // Generate vertices

    grid.line_width = line_width;
    grid.vertices.clear();

    const Vector4f line_color_major = Vector4f(0.8, 0.8, 0.8, 1);
    const Vector4f line_color_minor = Vector4f(0.5, 0.5, 0.5, 1);
    const Vector4f line_color_x = Vector4f(0.8, 0, 0, 1);
    const Vector4f line_color_y = Vector4f(0, 0.8, 0, 1);
    const Vector4f line_color_z = Vector4f(0, 0, 0.8, 1);

    auto create_line = [&](const Vector3f& a, const Vector3f& b, const Vector4f color) {
        grid.vertices.push_back(Vertex(a, color));
        grid.vertices.push_back(Vertex(b, color));
    };

    // Vertical line at origin
    create_line(Vector3f(0, 0, 0), Vector3f(0, 0, width), line_color_z);

    // Other lines
    int num_divisions = num_minor_divisions * num_major_divisions;
    num_divisions += (num_divisions % 2); // Ensure it is even

    for (int i = -num_divisions/2; i <= num_divisions/2; i++) {
        bool is_major = i % num_minor_divisions == 0;

        Vector4f line_color_i_x, line_color_i_y;
        if (i == 0) {
            line_color_i_x = line_color_x;
            line_color_i_y = line_color_y;
        } else if (is_major) {
            line_color_i_x = line_color_major;
            line_color_i_y = line_color_major;
        } else {
            line_color_i_x = line_color_minor;
            line_color_i_y = line_color_minor;
        }

        float pos = i * width / num_divisions;
        create_line(Vector3f(-width/2, pos, 0), Vector3f(width/2, pos, 0), line_color_i_x);
        create_line(Vector3f(pos, -width/2, 0), Vector3f(pos, width/2, 0), line_color_i_y);
    }

    // Load buffer data

    glBindVertexArray(grid.static_VAO);

    glBindBuffer(GL_ARRAY_BUFFER, grid.static_VBO);
    glBufferData(
        GL_ARRAY_BUFFER,
        grid.vertices.size() * sizeof(Vertex),
        &grid.vertices[0],
        GL_STATIC_DRAW
    );

    // Position
    glVertexAttribPointer(
        0, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex),
        (void*)offsetof(Vertex, position)
    );
    glEnableVertexAttribArray(0);

    // Color
    glVertexAttribPointer(
        1, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex),
        (void*)offsetof(Vertex, color)
    );
    glEnableVertexAttribArray(1);
}

void GridRenderer::queue_grid(int grid_index)
{
    if (grid_index < 0 || grid_index >= grids.size()) {
        return;
    }
    grids[grid_index].queued = true;
}

void GridRenderer::render(
    const Matrix4f& view,
    const Matrix4f& projection)
{
    glUseProgram(program_id);
    Matrix4f mvp = projection * view;
    glUniformMatrix4fv(mvp_loc, 1, GL_FALSE, mvp.data());

    for (const auto& grid: grids) {
        if (!grid.initialised) {
            continue;
        }
        if (!grid.queued) {
            continue;
        }
        grid.queued = false;

        glBindVertexArray(grid.static_VAO);
        glLineWidth(grid.line_width);
        glDrawArrays(
            GL_LINES,
            0,
            grid.vertices.size()
        );
    }
}

} // namespace rbm
