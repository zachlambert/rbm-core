#pragma once

#include "rbm/gui/renderer.h"


namespace rbm {

class GridRenderer: public Renderer {
public:
    GridRenderer();

    void render(
        const Matrix4f& view,
        const Matrix4f& projection) override;

    int create_grid();
    void configure_grid(int grid, double width, double line_width, int num_major_divisions, int num_minor_divisions);
    void queue_grid(int grid);

private:
    struct Vertex {
        Vector3f position;
        Vector4f color;
        Vertex() {}
        Vertex(const Vector3f& position, const Vector4f& color):
            position(position),
            color(color)
        {}
    };

    unsigned int program_id;
    unsigned int mvp_loc;

    struct Grid {
        bool initialised = false;
        mutable bool queued = false;
        unsigned int static_VAO;
        unsigned int static_VBO;
        float line_width;
        std::vector<Vertex> vertices;
    };
    std::vector<Grid> grids;
};

} // namespace rbm
