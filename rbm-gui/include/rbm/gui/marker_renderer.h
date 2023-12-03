#pragma once

#include <thread>
#include <rbm/geometry/primitive.h>
#include <rbm/geometry/mesh.h>
#include <rbm/types/matrix.h>
#include <rbm/types/color.h>
#include <rbm/types/interval.h>
#include "rbm/gui/renderer.h"


namespace rbm {

class MarkerRenderer: public Renderer {
public:
    MarkerRenderer();
    ~MarkerRenderer();

    void queue_primitive(const InstancedPrimitive3d& primitive, const ColorRGBd& color);

    // Helper methods
    void queue_line(const Vector3d& begin, const Vector3d& end, double width, const ColorRGBd& color);
    void queue_arrow(const Vector3d& begin, const Vector3d& end, double width, double head_length, double head_width, const ColorRGBd& color);
    void queue_frame(const Transform3d& pose, double axes_length, double axes_width);
    void queue_bounding_box(const Transform3d& pose, const BoundingBox3d& bounding_box, const ColorRGBd& color);

    void render(
        const Eigen::Matrix4f& view,
        const Eigen::Matrix4f& projection) override;

private:
    void generate_meshes();
    bool init_opengl();

    std::jthread generate_meshes_thread;
    std::atomic<bool> running;
    std::atomic<bool> generate_done;
    bool init_done;

    float resolution;

    using Mesh = ::rbm::Mesh<
        make_vertex<float, 3,
            vertex_attributes::position |
            vertex_attributes::normal>,
        unsigned short,
        3
    >;

    struct Command {
        ColorRGBf color;
        Eigen::Matrix4f model;
        std::size_t mesh_datas_index;
    };
    mutable std::vector<Command> commands;

    struct MeshData {
        std::size_t primitive;
        MeshRange mesh_range;
    };
    std::vector<MeshData> mesh_datas;
    Mesh mesh_data;

    struct {
        // OpenGL shader parameter locations
        unsigned int program_id;
        unsigned int mvp_loc;
        unsigned int mv_loc;
        unsigned int color_loc;
        // OpenGL array/buffer objects
        unsigned int static_VAO, static_VBO, static_EBO;
    } glData;
};

} // namespace rbm
