#pragma once

// TODO: Do other primitives too

#include "owl/geometry/collision.h"
#include "owl/geometry/primitive.h"
#include <iostream>

namespace owl {

template <typename Scalar>
struct collision_details<Edge<Scalar, 3>, Triangle<Scalar, 3>> {
    typedef Vector<Scalar, 3> intersection_type;

    static bool intersects(const Edge<Scalar, 3>& edge, const Triangle<Scalar, 3>& triangle) {
        return intersection(edge, triangle).has_value();
    }

    static std::optional<Vector<Scalar, 3>> intersection(const Edge<Scalar, 3>& edge, const Triangle<Scalar, 3>& triangle) {
        Vector<Scalar, 3> normal = (triangle.vertices[1] - triangle.vertices[0]).cross(triangle.vertices[2] - triangle.vertices[0]);
        normal.normalize();

        // Intersects if:
        // 1. Edge passes through plane of the triangle
        // 2. Point on the plane is to the "left" of every face edge (inside the face)

        // If these have the same sign (such that product is positive), then both
        // endpoints of the edge lie on the same side of the plane
        Scalar perp_0 = (edge.vertices[0] - triangle.vertices[0]).dot(normal);
        Scalar perp_1 = (edge.vertices[1] - triangle.vertices[0]).dot(normal);
        if (perp_0 * perp_1 >= 0) return std::nullopt;

        Scalar fraction = std::abs(perp_0) / (std::abs(perp_0) + std::abs(perp_1));
        Vector<Scalar, 3> point = edge.vertices[0] * (1 - fraction) + edge.vertices[1] * fraction;

        for (std::size_t i = 0; i < 3; i++) {
            Vector<Scalar, 3> inward_vector = normal.cross(triangle.vertices[(i + 1) % 3] - triangle.vertices[i]);
            inward_vector.normalize();
            Scalar inward_component = (point - triangle.vertices[i]).dot(inward_vector);
            if (inward_component < 0) return std::nullopt;
        }
        return point;
    }
};

} // namespace owl

template <typename Scalar>
struct collision_details<owl::Triangle<Scalar, 3>, owl::Triangle<Scalar, 3>> {
    typedef owl::Edge<Scalar, 3> intersection_type;

    static bool intersects(const owl::Triangle<Scalar, 3>& a, const owl::Triangle<Scalar, 3>& b) {
        int intersection_count = 0;
        // Query edges of a with triangle b
        for (std::size_t i = 0; i < 3; i++) {
            if (intersection_count == 2) break;
            owl::Edge<Scalar, 3> edge;
            edge.vertices[0] = a.vertices[i];
            edge.vertices[1] = a.vertices[(i + 1) % 3];
            intersection_count += intersects(edge, b);
        }
        // Query edges of b with triangle a
        for (std::size_t i = 0; i < 3; i++) {
            if (intersection_count == 2) break;
            owl::Edge<Scalar, 3> edge;
            edge.vertices[0] = b.vertices[i];
            edge.vertices[1] = b.vertices[(i + 1) % 3];
            intersection_count += intersects(edge, a);
        }
        return intersection_count == 2;
    }

    static std::optional<owl::Edge<Scalar, 3>> intersection(
            const owl::Triangle<Scalar, 3>& a,
            const owl::Triangle<Scalar, 3>& b)
    {
        int intersection_count = 0;
        owl::Edge<Scalar, 3> result;
        // Query edges of a with triangle b
        for (std::size_t i = 0; i < 3; i++) {
            if (intersection_count == 2) break;
            owl::Edge<Scalar, 3> edge;
            edge.vertices[0] = a.vertices[i];
            edge.vertices[1] = a.vertices[(i + 1) % 3];
            auto intersection = ::math::intersection(edge, b);
            if (intersection.has_value()) {
                result.vertices[intersection_count] = intersection.value();
                intersection_count++;
            }
        }
        // Query edges of b with triangle a
        for (std::size_t i = 0; i < 3; i++) {
            if (intersection_count == 2) break;
            owl::Edge<Scalar, 3> edge;
            edge.vertices[0] = b.vertices[i];
            edge.vertices[1] = b.vertices[(i + 1) % 3];
            auto intersection = ::math::intersection(edge, a);
            if (intersection.has_value()) {
                result.vertices[intersection_count] = intersection.value();
                intersection_count++;
            }
        }
        if (intersection_count == 2) {
            return result;
        }
        return std::nullopt;
    }
};

// namespace math
