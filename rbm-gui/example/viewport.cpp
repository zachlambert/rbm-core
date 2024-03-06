#include <rbm/geometry/primitive.h>
#include <rbm/geometry/primitive_conversions.h>
#include <rbm/geometry/visual_mesh.h>
#include <rbm/gui/window.h>
#include <rbm/gui/viewport.h>
#include <rbm/gui/grid.h>
#include <rbm/gui/mesh_entity.h>

int main()
{
    rbm::Sphere3f sphere;
    sphere.radius = 1;
    auto mesh = std::make_shared<rbm::VisualMesh>();
    rbm::primitive_to_mesh(sphere, 0.5, *mesh);
    for (auto& vertex: mesh->vertices) {
        vertex.color = rbm::ColorRGBf::Red();
    }

    rbm::Window window("render");

    auto viewport = std::make_shared<rbm::Viewport>(1200, 900);
    viewport->add_entity("grid", std::make_shared<rbm::GridEntity>());
    viewport->add_entity("mesh", std::make_shared<rbm::MeshEntity>(mesh));

    window.add_widget("viewport", viewport);
    window.run();
    return 0;
}
