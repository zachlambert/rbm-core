#pragma once

#include <rbm/types/color.h>
#include <rbm/geometry/visual_mesh.h>
#include "rbm/gui/renderer.h"


namespace rbm {

class MeshRenderer: public Renderer {
public:
    MeshRenderer();

    int load_mesh(const VisualMesh& mesh);
    bool queue_mesh(
        int mesh,
        const Transform3d& pose,
        bool wireframe);

    void render(
        const Matrix4f& view,
        const Matrix4f& projection) override;

private:
    struct Command {
        int mesh_index;
        Matrix4f model;
        bool wireframe;
    };
    std::vector<Command> commands;

    unsigned int program_id;
    unsigned int mvp_loc;
    unsigned int m_loc;
    unsigned int static_VAO, static_VBO, static_EBO;

    VisualMesh mesh_data;
    std::vector<MeshRange> mesh_ranges;
};

} // namespace rbm
