#pragma once

#include <rbm/geometry/visual_mesh.h>
#include "rbm/gui/entity.h"
#include "rbm/gui/mesh_renderer.h"


namespace rbm {

class MeshEntity: public Entity {
public:
    MeshEntity(const std::shared_ptr<const VisualMesh>& mesh):
        mesh(mesh),
        pose(Transform3d::Identity()),
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

    void set_pose(const Transform3d& pose)
    {
        this->pose = pose;
    }

private:
    std::shared_ptr<const VisualMesh> mesh;
    int mesh_index;
    Transform3d pose;
};

} // namespace rbm
