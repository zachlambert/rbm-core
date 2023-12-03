#pragma once

#include "sviz/render/entity.h"
#include "sviz/render/mesh_renderer.h"
#include <mbox/geometry/visual_mesh.h>

namespace sviz {

class Mesh: public Entity {
public:
    Mesh(const std::shared_ptr<const mbox::VisualMesh>& mesh):
        mesh(mesh),
        pose(mbox::Transform3d::Identity()),
        mesh_index(-1)
    {}

    bool init(Renderers& renderers)
    {
        mesh_index = renderers.get<MeshRenderer>()->load_mesh(*mesh);
        return true;
    }
    void render(Renderers& renderers)
    {
        renderers.get<MeshRenderer>()->queue_mesh(mesh_index, pose, false);
    }

    void set_pose(const mbox::Transform3d& pose)
    {
        this->pose = pose;
    }

private:
    std::shared_ptr<const mbox::VisualMesh> mesh;
    int mesh_index;
    mbox::Transform3d pose;
};

} // namespace sviz
