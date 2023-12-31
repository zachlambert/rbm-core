#include "rbm/gui/marker_renderer.h"

#include <GL/glew.h>
#include <rbm/geometry/primitive_conversions.h>
#include <rbm/transform/miscellaneous.h>
#include "rbm/gui/shader.h"
#include "rbm/gui/shader/marker_fs.h"
#include "rbm/gui/shader/marker_vs.h"


namespace rbm {

MarkerRenderer::MarkerRenderer():
    resolution(0.2),
    init_done(false),
    generate_done(false),
    running(true)
{
    generate_meshes_thread = std::jthread(std::bind(&MarkerRenderer::generate_meshes, this));
}

MarkerRenderer::~MarkerRenderer()
{
    running = false;
}

void MarkerRenderer::generate_meshes() {
    Mesh primitive_mesh;
    for (std::size_t primitive_index = 0; primitive_index < rbm::variant_count<InstancedPrimitive3f>; primitive_index++) {
        if (!running) return;

        InstancedPrimitive3f primitive = std::visit([](const auto& value) -> InstancedPrimitive3f{
            typedef std::decay_t<decltype(value)> T;
            return T::Instance();

        }, rbm::variant_construct<InstancedPrimitive3f>(primitive_index));

        primitive_to_mesh(primitive, resolution, primitive_mesh);

        MeshData data;
        data.primitive = primitive_index;
        data.mesh_range = insert_mesh(primitive_mesh, mesh_data);
        mesh_datas.push_back(data);
    }
    generate_done = true;
}

bool MarkerRenderer::init_opengl() {
    if (!generate_done) return false;

    // Get a program_id and uniform locations
    glData.program_id = load_shader(
        shader::marker_vs,
        shader::marker_fs
    );
    glData.mvp_loc = glGetUniformLocation(glData.program_id, "MVP");
    glData.mv_loc = glGetUniformLocation(glData.program_id, "MV");
    glData.color_loc = glGetUniformLocation(glData.program_id, "color");

    // Generate ids
    glGenVertexArrays(1, &glData.static_VAO);
    glGenBuffers(1, &glData.static_VBO);
    glGenBuffers(1, &glData.static_EBO);

    // Bind vertex array
    glBindVertexArray(glData.static_VAO);

    // Bind and configure buffer for vertex attributes
    glBindBuffer(GL_ARRAY_BUFFER, glData.static_VBO);
    glBufferData(
        GL_ARRAY_BUFFER,
        mesh_data.vertices.size() * sizeof(typename Mesh::Vertex),
        &mesh_data.vertices[0],
        GL_STATIC_DRAW
    );

    // Bind and configure buffer for indices
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, glData.static_EBO);
    glBufferData(
        GL_ELEMENT_ARRAY_BUFFER,
        mesh_data.facets.size() * 3 * sizeof(Mesh::IndexType),
        &mesh_data.facets[0].indices[0],
        GL_STATIC_DRAW
    );
    static_assert(Mesh::VertexCount == 3);
    static_assert(sizeof(Mesh::Facet) == 3 * sizeof(Mesh::IndexType));

    // Define vertex buffer attributes

    // Position
    glVertexAttribPointer(
        0, 3, GL_FLOAT, GL_FALSE, sizeof(Mesh::Vertex),
        (void*)offsetof(Mesh::Vertex, position)
    );
    glEnableVertexAttribArray(0);
    // Normal
    glVertexAttribPointer(
        1, 3, GL_FLOAT, GL_FALSE, sizeof(Mesh::Vertex),
        (void*)offsetof(Mesh::Vertex, normal)
    );
    glEnableVertexAttribArray(1);

    return true;
}

void MarkerRenderer::queue_primitive(const InstancedPrimitive3d& primitive, const ColorRGBd& color) {
    if (!init_done) return;

    Command command;
    command.color = color.cast<float>();
    command.model = std::visit([](const auto& value) -> Matrix4f {
        return value.instance_transform().template cast<float>();
    }, primitive);

    std::size_t primitive_index = rbm::variant_index(primitive);
    for (std::size_t i = 0; i < mesh_datas.size(); i++) {
        const MeshData& mesh_data = mesh_datas[i];
        if (mesh_data.primitive == primitive_index) {
            command.mesh_datas_index = i;
            break;
        }
    }
    commands.push_back(command);
}

void MarkerRenderer::queue_arrow(
        const Vector3d& begin,
        const Vector3d& end,
        double width,
        double head_length,
        double head_width,
        const ColorRGBd& color)
{
    Vector3d direction = (end - begin).normalized();
    double length = (end - begin).norm();

    head_length = std::min(head_length, length);
    Vector3d line_end = end - direction * head_length;
    double line_length = length - head_length;

    Cylinder3d cylinder;
    cylinder.pose.translation() = begin + direction * line_length / 2;
    cylinder.pose.rotation() = rotation_from_normal(direction);
    cylinder.length = line_length;
    cylinder.radius = width / 2;
    queue_primitive(cylinder, color);

    Sphere3d back_cap;
    back_cap.position = begin;
    back_cap.radius = width / 2;
    queue_primitive(back_cap, color);

    Cone3d head;
    head.pose.translation() = begin + direction * line_length;
    head.pose.rotation() = rotation_from_normal(direction);
    head.radius = head_width / 2;
    head.length = head_length;
    queue_primitive(head, color);
}

void MarkerRenderer::queue_frame(
        const Transform3d& pose,
        double axes_length,
        double axes_width)
{
    Vector3d origin = pose.translation();
    Vector3d end_x = origin + axes_length * pose.rotation().block<3, 1>(0, 0);
    Vector3d end_y = origin + axes_length * pose.rotation().block<3, 1>(0, 1);
    Vector3d end_z = origin + axes_length * pose.rotation().block<3, 1>(0, 2);
    queue_line(origin, end_x, axes_width, ColorRGBd::Red());
    queue_line(origin, end_y, axes_width, ColorRGBd::Green());
    queue_line(origin, end_z, axes_width, ColorRGBd::Blue());
}

void MarkerRenderer::queue_bounding_box(
        const Transform3d& pose,
        const BoundingBox3d& bounding_box,
        const ColorRGBd& color)
{
    Box3d box;
    box.pose.translation() = pose * bounding_box.centre();
    box.pose.rotation() = pose.rotation();
    box.size = bounding_box.size();
    for (auto& scalar: box.size) {
        scalar = std::max(scalar, 1e-6);
    }
    queue_primitive(box, color);
}

void MarkerRenderer::queue_line(const Vector3d& begin, const Vector3d& end, double width, const ColorRGBd& color) {
    Cylinder3d cylinder;
    cylinder.pose.translation() = (begin + end) / 2;
    cylinder.pose.rotation() = rotation_from_normal((end - begin).normalized());
    cylinder.length = (end - begin).norm();
    cylinder.radius = width / 2;
    queue_primitive(cylinder, color);

    Sphere3d cap;
    cap.position = begin;
    cap.radius = width / 2;
    queue_primitive(cap, color);
    cap.position = end;
    cap.radius = width / 2;
    queue_primitive(cap, color);
}

void MarkerRenderer::render(
    const Matrix4f& view,
    const Matrix4f& projection)
{
    if (!init_done && !init_opengl()) {
        return;
    }
    init_done = true;

    glUseProgram(glData.program_id);
    glBindVertexArray(glData.static_VAO);

    // TODO: Use instanced rendering instead

    Matrix4f mvp;
    Matrix4f mv;
    for (const auto& command: commands) {
        mvp = projection * view * command.model;
        mv = view * command.model;
        glUniformMatrix4fv(glData.mvp_loc, 1, GL_FALSE, &mvp(0, 0));
        glUniformMatrix4fv(glData.mv_loc, 1, GL_FALSE, &mv(0, 0));
        glUniform4f(glData.color_loc, command.color.r, command.color.g, command.color.b, 1);

        const MeshData& mesh_data = mesh_datas[command.mesh_datas_index];

        glDrawElements(
            GL_TRIANGLES,
            mesh_data.mesh_range.facet_count * 3,
            GL_UNSIGNED_SHORT,
            (void*)(sizeof(unsigned short) * mesh_data.mesh_range.facets_offset * 3)
        );
    }
    commands.clear();
}

} // namespace rbm
