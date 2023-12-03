#include "sviz/render/mesh_renderer.h"

#include "shader/mesh_vs.h"
#include "shader/mesh_fs.h"
#include "sviz/render/shader.h"


namespace sviz {

MeshRenderer::MeshRenderer()
{
    program_id = load_shader(
        shader::mesh_vs,
        shader::mesh_fs
    );
    mvp_loc = glGetUniformLocation(program_id, "MVP");
    m_loc = glGetUniformLocation(program_id, "M");

    glGenVertexArrays(1, &static_VAO);
    glGenBuffers(1, &static_VBO);
    glGenBuffers(1, &static_EBO);

    // Bind vertex array and vertex buffer
    glBindVertexArray(static_VAO);
    glBindBuffer(GL_ARRAY_BUFFER, static_VBO);

    glVertexAttribPointer(
        0, 3, GL_FLOAT, GL_FALSE, sizeof(mbox::VisualMesh::Vertex),
        (void*)offsetof(mbox::VisualMesh::Vertex, position)
    );
    glEnableVertexAttribArray(0);

    glVertexAttribPointer(
        1, 3, GL_FLOAT, GL_FALSE, sizeof(mbox::VisualMesh::Vertex),
        (void*)offsetof(mbox::VisualMesh::Vertex, normal)
    );
    glEnableVertexAttribArray(1);

    glVertexAttribPointer(
        2, 3, GL_FLOAT, GL_FALSE, sizeof(mbox::VisualMesh::Vertex),
        (void*)offsetof(mbox::VisualMesh::Vertex, color)
    );
    glEnableVertexAttribArray(2);
    static_assert(sizeof(mbox::ColorRGBf) == 3 * sizeof(float));
}

int MeshRenderer::load_mesh(const mbox::VisualMesh& mesh) {
    mbox::MeshRange mesh_range = insert_mesh(mesh, mesh_data);
    mesh_ranges.emplace_back(mesh_range);

    // Rebind vertices and indices

    glBindVertexArray(static_VAO);

    // Bind and configure buffer for vertex attributes
    glBindBuffer(GL_ARRAY_BUFFER, static_VBO);
    glBufferData(
        GL_ARRAY_BUFFER,
        mesh_data.vertices.size() * sizeof(mbox::VisualMesh::Vertex),
        &mesh_data.vertices[0],
        GL_STATIC_DRAW
    );

    // Bind and configure buffer for indices
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, static_EBO);
    glBufferData(
        GL_ELEMENT_ARRAY_BUFFER,
        mesh_data.facets.size() * mbox::VisualMesh::VertexCount * sizeof(mbox::VisualMesh::IndexType),
        &mesh_data.facets[0].indices[0],
        GL_STATIC_DRAW
    );
    static_assert(mbox::VisualMesh::VertexCount == 3);
    static_assert(sizeof(mbox::VisualMesh::Facet) == mbox::VisualMesh::VertexCount * sizeof(mbox::VisualMesh::IndexType));

    return mesh_ranges.size() - 1;
}

bool MeshRenderer::queue_mesh(
    int mesh,
    const mbox::Transform3d& pose,
    bool wireframe)
{
    if (mesh < 0 || mesh >= mesh_ranges.size()) return false;
    Command command;
    command.mesh_index = mesh;
    command.model = pose.cast<float>().matrix();
    command.wireframe = wireframe;
    commands.push_back(command);
    return true;
}

void MeshRenderer::render(
    const mbox::Matrix4f& view,
    const mbox::Matrix4f& projection)
{
    mbox::Matrix4f mvp;
    mbox::Matrix4f m;

    // Visual meshes

    glUseProgram(program_id);
    glBindVertexArray(static_VAO);

    for (const auto& command: commands) {
        mvp = projection * view * command.model;
        m = command.model;
        glUniformMatrix4fv(mvp_loc, 1, GL_FALSE, &mvp(0, 0));
        glUniformMatrix4fv(m_loc, 1, GL_FALSE, &m(0, 0));

        const mbox::MeshRange& mesh_range = mesh_ranges[command.mesh_index];

        if (command.wireframe) {
            glPolygonMode(GL_FRONT, GL_LINE);
            glPolygonMode(GL_BACK, GL_LINE);
        }
        glDrawElements(
            GL_TRIANGLES,
            mesh_range.facet_count * mbox::VisualMesh::VertexCount,
            GL_UNSIGNED_SHORT,
            (void*)(sizeof(mbox::VisualMesh::IndexType) * mesh_range.facets_offset * mbox::VisualMesh::VertexCount)
        );
        static_assert(std::is_same_v<unsigned short, mbox::VisualMesh::IndexType>);

        if (command.wireframe) {
            glPolygonMode(GL_FRONT, GL_FILL);
            glPolygonMode(GL_BACK, GL_FILL);
        }
    }

    commands.clear();
}

} // namespace sviz
