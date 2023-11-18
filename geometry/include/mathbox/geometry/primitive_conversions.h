#pragma once

#include "mbox/geometry/primitive.h"
#include "mbox/geometry/mesh.h"
#include "mbox/geometry/vertex.h"


namespace mbox {

template <typename Scalar, int Dim, typename Mesh>
requires has_position<Scalar, Dim, typename Mesh::Vertex>
void primitive_to_mesh(const Box<Scalar, Dim>& box, double resolution, Mesh& mesh) {
    static_assert(Dim == 3);
    static_assert(Mesh::VertexCount == 3);
    typedef typename Mesh::Vertex Vertex;
    typedef typename Mesh::Facet Facet;

    mesh.vertices.clear();
    mesh.facets.clear();

    Facet facet1;
    Facet facet2;

    Vertex vertex;
    for (int i = 0; i < 6; i++) {
        int dim = i % 3;
        int dim2 = (dim + 1) % Dim;
        int dim3 = (dim + 2) % Dim;

        Scalar dir = i/3 == 0 ? -1 : 1;

        if constexpr (has_normal<Scalar, Dim, Vertex>) {
            vertex.normal(dim) = dir;
            vertex.normal(dim2) = 0;
            vertex.normal(dim3) = 0;
        }

        vertex.position(dim) = dir * box.size(dim) / 2;
        for (int j = 0; j < 4; j++) {
            int j_ = i/3 == 0 ? (3 - j) : j;
            Scalar dir2 = j_%2 == 0 ? -1 : 1;
            Scalar dir3 = (j_/2 == 0 ? -1 : 1) * dir;
            vertex.position(dim2) = dir2 * box.size(dim2) / 2;
            vertex.position(dim3) = dir3 * box.size(dim3) / 2;

            mesh.vertices.push_back(vertex);
        }

        for (int ind = 0; ind < 3; ind++) {
            facet1.indices[ind] = i * 4 + ind;
            facet2.indices[ind] = (i + 1) * 4 - 1 - ind;
        }
        mesh.facets.push_back(facet1);
        mesh.facets.push_back(facet2);
    }

    for (auto& vertex: mesh.vertices) {
        vertex.position = box.pose * vertex.position;
        if constexpr (has_normal<Scalar, Dim, Vertex>) {
            vertex.normal = box.pose.rotation() * vertex.normal;
        }
    }
}

template <typename Scalar, int Dim, typename Mesh>
requires has_position<Scalar, Dim, typename Mesh::Vertex>
void primitive_to_mesh(const Sphere<Scalar, Dim>& sphere, double resolution, Mesh& mesh) {
    static_assert(Dim == 3);
    static_assert(Mesh::VertexCount == 3);
    typedef typename Mesh::Vertex Vertex;
    typedef typename Mesh::Facet Facet;
    typedef Vector<Scalar, Dim> Vector;

    mesh.vertices.clear();
    mesh.facets.clear();

    const size_t N = std::max((size_t)std::ceil(2 * M_PI * sphere.radius / resolution), 2lu);
    Scalar phi, theta;

    Vector direction;
    Vertex vertex;
    vertex.position = sphere.radius * Vector(0, 0, 1);
    if constexpr (has_normal<Scalar, Dim, Vertex>) {
        vertex.normal = Vector(0, 0, 1);
    }
    mesh.vertices.push_back(vertex);
    for (std::size_t i = 1; i < N; i++) {
        phi = (M_PI*static_cast<Scalar>(i))/N;
        direction.z() = std::cos(phi);
        for (std::size_t j = 0; j < 2*N; j++) {
            theta = (2*M_PI*static_cast<Scalar>(j))/(2*N);
            direction.x() = std::sin(phi) * std::cos(theta);
            direction.y() = std::sin(phi) * std::sin(theta);
            vertex.position = sphere.radius * direction;
            if constexpr (has_normal<Scalar, Dim, Vertex>) {
                vertex.normal = direction;
            }
            mesh.vertices.push_back(vertex);
        }
    }
    vertex.position = sphere.radius * Vector(0, 0, -1);
    if constexpr (has_normal<Scalar, Dim, Vertex>) {
        vertex.normal = Vector(0, 0, -1);
    }
    mesh.vertices.push_back(vertex);

    // Append to facets

    Facet facet;

    // Top strip, common vertex = top
    for (std::size_t j = 0; j < 2*N; j++) {
        facet.indices[0] = 0;
        facet.indices[1] = 1 + j;
        if (j < 2*N-1) {
            facet.indices[2] = 2+j;
        } else {
            facet.indices[2] = 1;
        }
        mesh.facets.push_back(facet);
    }
    // Bottom strip, common vertex = bottom
    std::size_t bot = mesh.vertices.size() - 1;
    for (std::size_t j = 0; j < 2*N; j++) {
        facet.indices[0] = bot;
        facet.indices[1] = bot-1-j;
        if (j < 2*N-1) {
            facet.indices[2] = bot-2-j;
        } else {
            facet.indices[2] = bot-1;
        }
        mesh.facets.push_back(facet);
    }
    // Remaining strips, made up of rectangles between
    std::size_t strip_1_start, strip_2_start;
    for (std::size_t i = 0; i < N-2; i++) {
        strip_1_start = 1 + 2*N*i;
        strip_2_start = 1 + 2*N*(i+1);
        for (std::size_t j = 0; j < 2*N; j++) {
            // First triangle of rectangle
            facet.indices[0] = strip_1_start + j;
            facet.indices[1] = strip_2_start + j;
            if (j < 2*N-1) {
                facet.indices[2] = strip_2_start + j + 1;
            } else {
                facet.indices[2] = strip_2_start;
            }
            mesh.facets.push_back(facet);
            // Second triangle of rectangle
            if (j < 2*N-1) {
                facet.indices[0] = strip_2_start + j + 1;
                facet.indices[1] = strip_1_start + j + 1;
            } else {
                facet.indices[0] = strip_2_start;
                facet.indices[1] = strip_1_start;
            }
            facet.indices[2] = strip_1_start + j;
            mesh.facets.push_back(facet);
        }
    }

    for (auto& vertex: mesh.vertices) {
        vertex.position = sphere.position + vertex.position;
    }
}

template <typename Scalar, int Dim, typename Mesh>
requires has_position<Scalar, Dim, typename Mesh::Vertex>
void primitive_to_mesh(const Cylinder<Scalar, Dim>& cylinder, double resolution, Mesh& mesh) {
    static_assert(Dim == 3);
    static_assert(Mesh::VertexCount == 3);
    typedef typename Mesh::Vertex Vertex;
    typedef typename Mesh::Facet Facet;
    typedef Vector<Scalar, Dim> Vector;

    mesh.vertices.clear();
    mesh.facets.clear();

    const std::size_t N = std::max((size_t)std::ceil(2 * M_PI * cylinder.radius / resolution), 2lu);

    Vertex vertex;
    Facet facet;

    const Vector top_centre = Vector::UnitZ() * cylinder.length / 2;
    const Vector bot_centre = -Vector::UnitZ() * cylinder.length / 2;

    // Top face

    std::size_t start = 0;
    vertex.position = top_centre;
    if constexpr (has_normal<Scalar, Dim, Vertex>) {
        vertex.normal = Vector(0, 0, 1);
    }
    mesh.vertices.push_back(vertex);
    for (std::size_t i = 0; i < N; i++) {
        Scalar theta = i*2*M_PI / N;
        vertex.position = top_centre + cylinder.radius * Vector(std::cos(theta), std::sin(theta), 0);
        mesh.vertices.push_back(vertex);
    }
    for (std::size_t i = 0; i < N; i++) {
        facet.indices[0] = start;
        facet.indices[1] = start + 1 + i;
        facet.indices[2] = start + 1 + (i+1)%N;
        mesh.facets.push_back(facet);
    }

    // Curved surface

    start = mesh.vertices.size();
    for (std::size_t i = 0; i < N; i++) {
        Scalar theta = i*2*M_PI / N;
        if constexpr (has_normal<Scalar, Dim, Vertex>) {
            vertex.normal = Vector(std::cos(theta), std::sin(theta), 0);
        }
        vertex.position = top_centre + cylinder.radius * Vector(std::cos(theta), std::sin(theta), 0);
        mesh.vertices.push_back(vertex);
        vertex.position = bot_centre + cylinder.radius * Vector(std::cos(theta), std::sin(theta), 0);
        mesh.vertices.push_back(vertex);
    }
    for (std::size_t i = 0; i < 2*N; i+=2) {
        facet.indices[0] = start + i;
        facet.indices[1] = start + i+1;
        facet.indices[2] = start + (i+2)%(2*N);
        mesh.facets.push_back(facet);
        facet.indices[0] = start + i+1;
        facet.indices[1] = start + (i+3)%(2*N);
        facet.indices[2] = start + (i+2)%(2*N);
        mesh.facets.push_back(facet);
    }

    // Bottom face

    start = mesh.vertices.size();
    vertex.position = bot_centre;
    if constexpr (has_normal<Scalar, Dim, Vertex>) {
        vertex.normal = Vector(0, 0, -1);
    }
    mesh.vertices.push_back(vertex);
    for (std::size_t i = 0; i < N; i++) {
        Scalar theta = i*2*M_PI / N;
        vertex.position = bot_centre + cylinder.radius * Vector(std::cos(-theta), std::sin(-theta), 0);
        mesh.vertices.push_back(vertex);
    }
    for (std::size_t i = 0; i < N; i++) {
        facet.indices[0] = start;
        facet.indices[1] = start + 1 + i;
        facet.indices[2] = start + 1 + (i+1)%N;
        mesh.facets.push_back(facet);
    }

    for (auto& vertex: mesh.vertices) {
        vertex.position = cylinder.pose * vertex.position;
        if constexpr (has_normal<Scalar, Dim, Vertex>) {
            vertex.normal = cylinder.pose.rotation() * vertex.normal;
        }
    }
}

template <typename Scalar, int Dim, typename Mesh>
requires has_position<Scalar, Dim, typename Mesh::Vertex>
void primitive_to_mesh(const Cone<Scalar, Dim>& cone, double resolution, Mesh& mesh) {
    static_assert(Dim == 3);
    static_assert(Mesh::VertexCount == 3);
    typedef typename Mesh::Vertex Vertex;
    typedef typename Mesh::Facet Facet;
    typedef Vector<Scalar, Dim> Vector;

    mesh.vertices.clear();
    mesh.facets.clear();

    const std::size_t N = std::max((size_t)std::ceil(2 * M_PI * cone.radius / resolution), 2lu);

    Vertex vertex;
    Facet facet;

    // Bottom face

    std::size_t start = 0;
    vertex.position.setZero();
    if constexpr (has_normal<Scalar, Dim, Vertex>) {
        vertex.normal = Vector(0, 0, -1);
    }
    mesh.vertices.push_back(vertex);
    for (std::size_t i = 0; i < N; i++) {
        Scalar theta = 2 * M_PI * static_cast<Scalar>(i) / N;
        vertex.position = cone.radius * Vector(std::cos(-theta), std::sin(-theta), 0);
        mesh.vertices.push_back(vertex);
    }
    for (size_t i = 0; i < N; i++) {
        facet.indices[0] = start;
        facet.indices[1] = start + 1 + i;
        facet.indices[2] = start + 1 + (i+1)%N;
        mesh.facets.push_back(facet);
    }

    // Curved surface

    Scalar pitch = std::atan(cone.radius / cone.length);
    Vector normal(cos(pitch), sin(pitch), 0);

    start = mesh.vertices.size();
    for (std::size_t i = 0; i < N; i++) {
        Scalar theta = 2*M_PI*static_cast<Scalar>(i) / N;
        Scalar theta_plus_half = 2*M_PI*(static_cast<Scalar>(i) + 0.5) / N;
        vertex.position = cone.radius * Vector(std::cos(theta), std::sin(theta), 0);
        if constexpr (has_normal<Scalar, Dim, Vertex>) {
            vertex.normal = Eigen::AngleAxis<Scalar>(theta, Vector::UnitZ()) * normal;
        }
        mesh.vertices.push_back(vertex);
        vertex.position = cone.length * Vector::UnitZ();
        if constexpr (has_normal<Scalar, Dim, Vertex>) {
            vertex.normal = Eigen::AngleAxis<Scalar>(theta_plus_half, Vector::UnitZ()) * normal;
        }
        mesh.vertices.push_back(vertex);
    }
    for (std::size_t i = 0; i < (2*N); i+=2) {
        facet.indices[0] = start + i;
        facet.indices[1] = start + (i+2)%(2*N);
        facet.indices[2] = start + i+1;
        mesh.facets.push_back(facet);
    }

    for (auto& vertex: mesh.vertices) {
        vertex.position = cone.pose * vertex.position;
        if constexpr (has_normal<Scalar, Dim, Vertex>) {
            vertex.normal = cone.pose.rotation() * vertex.normal;
        }
    }
}

template <typename Scalar, int Dim, typename Mesh>
requires has_position<Scalar, Dim, typename Mesh::Vertex>
void primitive_to_mesh(const InstancedPrimitive<Scalar, Dim>& primitive, double resolution, Mesh& mesh) {
    std::visit([&mesh, &resolution](auto& value) {
        primitive_to_mesh(value, resolution, mesh);
    }, primitive);
}

} // namespace mbox
