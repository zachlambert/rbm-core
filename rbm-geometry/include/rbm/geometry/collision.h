#pragma once

#include <concepts>
#include <optional>
#include <rbm/types/matrix.h>

namespace rbm {

template <typename A, typename B>
struct collision_details {};

template <typename A, typename B>
concept has_collision_intersects = requires(const A& a, const B& b) {
    { collision_details<A, B>::intersects(a, b) } -> std::convertible_to<bool>;
};

template <typename A, typename B>
concept has_collision_intersection = requires(const A& a, const B& b) {
    // { collision_details<A, B>::intersection_type };
    { collision_details<A, B>::intersection(a, b) } -> std::convertible_to<std::optional<typename collision_details<A, B>::intersection_type>>;
};

template <typename A, typename B>
requires has_collision_intersects<A, B>
bool intersects(const A& a, const B& b) {
    return collision_details<A, B>::intersects(a, b);
}

template <typename A, typename B>
using intersection_type = typename collision_details<A, B>::intersection_type;

template <typename A, typename B>
requires has_collision_intersection<A, B>
std::optional<intersection_type<A, B>> intersection(const A& a, const B& b) {
    return collision_details<A, B>::intersection(a, b);
}

} // namespace rbm
