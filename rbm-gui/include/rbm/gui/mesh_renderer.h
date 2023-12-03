#pragma once

#include "sviz/render/renderer.h"
#include "mbox/types/color.h"
#include "mbox/geometry/visual_mesh.h"


namespace sviz {

class MeshRenderer: public Renderer {
public:
    MeshRenderer();

    int load_mesh(const mbox::VisualMesh& mesh);
    bool queue_mesh(
        int mesh,
        const mbox::Transform3d& pose,
        bool wireframe);

    void render(
        const mbox::Matrix4f& view,
        const mbox::Matrix4f& projection) override;

private:
    struct Command {
        int mesh_index;
        mbox::Matrix4f model;
        bool wireframe;
    };
    std::vector<Command> commands;

    unsigned int program_id;
    unsigned int mvp_loc;
    unsigned int m_loc;
    unsigned int static_VAO, static_VBO, static_EBO;

    mbox::VisualMesh mesh_data;
    std::vector<mbox::MeshRange> mesh_ranges;
};

} // namespace sviz
