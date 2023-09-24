#pragma once

#include <Eigen/Core>
#include <array>
#include <concepts>
#include <string>
#include <cstring>
#include "math/types/matrix.h"
#include "math/types/color.h"
#include "math/transform/transform.h"


namespace math  {

template <typename Scalar, int Dim, typename Vertex>
concept has_position = requires(const Vertex& vertex) {
    { vertex.position } -> std::convertible_to<Vector<Scalar, Dim>>;
};

template <typename Scalar, int Dim, typename Vertex>
concept has_normal = requires(const Vertex& vertex) {
    { vertex.normal } -> std::convertible_to<Vector<Scalar, Dim>>;
};

namespace vertex_attributes {
    static constexpr int position = 1<<0;
    static constexpr int normal = 1<<1;
    static constexpr int color = 1<<2;
    static constexpr int texture = 1<<3;
};

template <typename Scalar, int Dim, int Flags>
struct make_vertex {};

template <typename Scalar, int Dim>
struct make_vertex<Scalar, Dim, vertex_attributes::position> {
    Vector<Scalar, Dim> position;
};

template <typename Scalar, int Dim>
struct make_vertex<Scalar, Dim, vertex_attributes::position | vertex_attributes::normal> {
    Vector<Scalar, Dim> position;
    Vector<Scalar, Dim> normal;
};

template <typename Scalar, int Dim>
struct make_vertex<Scalar, Dim, vertex_attributes::position | vertex_attributes::normal | vertex_attributes::color> {
    Vector<Scalar, Dim> position;
    Vector<Scalar, Dim> normal;
    ColorRGB<Scalar> color;
};

template <typename Scalar, int Dim>
struct make_vertex<Scalar, Dim, vertex_attributes::position | vertex_attributes::normal | vertex_attributes::texture> {
    Vector<Scalar, Dim> position;
    Vector<Scalar, Dim> normal;
    std::array<Scalar, 2> uv;
};

} // namespace math
