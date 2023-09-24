#pragma once

#include <vector>
#include <algorithm>
#include <array>
#include <unordered_map>
#include "math/geometry/vertex.h"
#include "math/geometry/facet.h"
#include "math/geometry/primitive.h"


namespace math {

template <typename Vertex_, typename IndexType_, int VertexCount_>
struct Mesh {
    using Vertex = Vertex_;
    using IndexType = IndexType_;
    static constexpr int VertexCount = VertexCount_;
    using Facet = ::math::Facet<IndexType, VertexCount>;

    std::vector<Vertex> vertices;
    std::vector<Facet> facets;

    void clear() {
        vertices.clear();
        facets.clear();
    }
};

struct MeshRange {
    std::size_t facet_count;
    std::size_t facets_offset;
};

template <typename Mesh>
MeshRange insert_mesh(const Mesh& mesh, Mesh& data) {
    MeshRange range;
    range.facets_offset = data.facets.size();
    range.facet_count = mesh.facets.size();

    std::size_t vertices_offset = data.vertices.size();
    for (const auto& vertex: mesh.vertices) {
        data.vertices.push_back(vertex);
    }
    for (const auto& facet: mesh.facets) {
        typename Mesh::Facet new_facet = facet;
        for (auto& index: new_facet.indices) {
            index += vertices_offset;
        }
        data.facets.push_back(new_facet);
    }
    return range;
}

#if 0
template <typename Scalar, int Dim, typename Mesh>
void transform(const Transform<Scalar, Dim>& transform, Mesh& mesh) {
    static_assert(has_position<Scalar, Dim, Mesh::Vertex>)
    for (auto& vertex: mesh.vertices) {
        vertex.position = transform * vertex.position;
        if constexpr (has_normal<Scalar, Dim, Mesh::Vertex>) {
            vertex.normal = transform.rotation() * vertex.normal;
        }
    }
}

template <typename Scalar, int Dim, typename Mesh>
void scale(Scalar scaling, Mesh& mesh) {
    static_assert(has_position<Sacalar, Dim, Mesh::Vertex>)
    for (auto& vertex: mesh.vertices) {
        vertex.position.array() *= scaling;
    }
}
#endif

} // namespace math
