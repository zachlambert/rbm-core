#pragma once

#include "sviz/render/renderer.h"
#include "mbox/geometry/primitive.h"
#include "mbox/geometry/mesh.h"
#include "mbox/types/matrix.h"
#include "mbox/types/color.h"
#include "mbox/types/interval.h"
#include <thread>


namespace sviz {

class MarkerRenderer: public Renderer {
public:
    MarkerRenderer();
    ~MarkerRenderer();

    void queue_primitive(const mbox::InstancedPrimitive3d& primitive, const mbox::ColorRGBd& color);

    // Helper methods
    void queue_line(const mbox::Vector3d& begin, const mbox::Vector3d& end, double width, const mbox::ColorRGBd& color);
    void queue_arrow(const mbox::Vector3d& begin, const mbox::Vector3d& end, double width, double head_length, double head_width, const mbox::ColorRGBd& color);
    void queue_frame(const mbox::Transform3d& pose, double axes_length, double axes_width);
    void queue_bounding_box(const mbox::Transform3d& pose, const mbox::BoundingBox3d& bounding_box, const mbox::ColorRGBd& color);

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

    using Mesh = mbox::Mesh<
        mbox::make_vertex<float, 3,
            mbox::vertex_attributes::position |
            mbox::vertex_attributes::normal>,
        unsigned short,
        3
    >;

    struct Command {
        mbox::ColorRGBf color;
        Eigen::Matrix4f model;
        std::size_t mesh_datas_index;
    };
    mutable std::vector<Command> commands;

    struct MeshData {
        std::size_t primitive;
        mbox::MeshRange mesh_range;
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

} // namespace sviz
