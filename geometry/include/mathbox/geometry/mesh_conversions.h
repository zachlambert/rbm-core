#pragma once

#include "mathbox/geometry/mesh.h"
#include "mathbox/geometry/bvh.h"

namespace mbox {

template <typename Scalar, int Dim, typename Mesh>
using MeshPrimitive = FacetShape<Scalar, Dim, Mesh::VertexCount>;

template <typename Scalar, int Dim, typename Mesh>
requires has_position<Scalar, Dim, typename Mesh::Vertex>
void mesh_to_primitives(
    const Mesh& mesh,
    std::vector<MeshPrimitive<Scalar, Dim, Mesh>>& primitives)
{
    primitives.resize(mesh.facets.size());
    for (std::size_t facet_i = 0; facet_i < mesh.facets.size(); facet_i++) {
        const auto& facet = mesh.facets[facet_i];
        auto& primitive = primitives[facet_i];
        for (std::size_t i = 0; i < Mesh::VertexCount; i++) {
            primitive.vertices[i] = mesh.vertices[facet.indices[i]].position;
        }
    }
}

template <typename Scalar, int Dim, typename Mesh>
requires has_position<Scalar, Dim, typename Mesh::Vertex>
void mesh_to_bvh(const Mesh& mesh, Bvh<MeshPrimitive<Scalar, Dim, Mesh>>& bvh)
{
    std::vector<MeshPrimitive<Scalar, Dim, Mesh>> primitives;
    mesh_to_primitives(mesh, primitives);
    bvh.construct(primitives);
}

} // namespace mbox
