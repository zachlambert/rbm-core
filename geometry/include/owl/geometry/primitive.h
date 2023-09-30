#pragma once

#include <variant>
#include <array>
#include <concepts>
#include "cpp_utils/annotated_variant.h"
#include "owl/types/matrix.h"
#include "owl/types/color.h"
#include "owl/transform/transform.h"
#include "owl/transform/homogeneous.h"
#include "owl/geometry/bounding.h"

namespace owl {

template <typename T>
concept is_primitive = requires(T t, Transform<typename T::Scalar, T::Dim> transform, typename T::Scalar scaling) {
    // { typename T::Scalar };
    { T::Dim } -> std::convertible_to<int>;
    { t.bounding_box() } -> std::convertible_to<BoundingBox<typename T::Scalar, T::Dim>>;
    { t.transform(transform) };
};

template <typename A, typename B>
concept primitives_compatible = requires(A a, B b) {
    { is_primitive<A> };
    { is_primitive<B> };
    { std::is_same_v<typename A::Scalar, typename B::Scalar> };
    { A::Dim == B::Dim };
};

template <typename T>
concept is_instanced_primitive = requires(const T& t) {
    { is_primitive<T> };
    { T::Instance() } -> std::convertible_to<T>;
    { t.instance_transform() } -> std::convertible_to<Matrix<typename T::Scalar, T::Dim + 1, T::Dim + 1>>;
};

template <typename Scalar_, int Dim_, int VertexCount>
struct FacetShape {
    using Scalar = Scalar_;
    static constexpr int Dim = Dim_;
    std::array<Vector<Scalar, Dim>, VertexCount> vertices;

    Vector<Scalar, Dim> centroid() const {
        Vector<Scalar, Dim> sum;
        sum.setZero();
        for (const auto& vertex: vertices) {
            sum += vertex;
        }
        return sum / vertices.size();
    }

    BoundingBox<Scalar, Dim> bounding_box() const {
        BoundingBox<Scalar, Dim> result;
        for (const auto& vertex: vertices) {
            result.merge(vertex);
        }
        return result;
    }

    void transform(const Transform<Scalar, Dim>& transform) {
        for (auto& vertex: vertices) {
            vertex = transform * vertex;
        }
    }

    template <typename T_>
    FacetShape<T_, Dim, VertexCount> cast() const {
        FacetShape<T_, Dim, VertexCount> result;
        for (std::size_t i = 0; i < vertices.size(); i++) {
            result.vertices[i] = vertices[i].template cast<T_>();
        }
        return result;
    }
};

template <typename Scalar, int Dim>
using Edge = FacetShape<Scalar, Dim, 2>;
typedef Edge<float, 2> Edge2f;
typedef Edge<double, 2> Edge2d;
typedef Edge<float, 3> Edge3f;
typedef Edge<double, 3> Edge3d;

template <typename Scalar, int Dim>
using Triangle = FacetShape<Scalar, Dim, 3>;
typedef Triangle<float, 2> Triangle2f;
typedef Triangle<double, 2> Triangle2d;
typedef Triangle<float, 3> Triangle3f;
typedef Triangle<double, 3> Triangle3d;

template <typename Scalar, int Dim>
using Tetrahedron = FacetShape<Scalar, Dim, 4>;
typedef Tetrahedron<float, 2> Tetrahedron2f;
typedef Tetrahedron<double, 2> Tetrahedron2d;
typedef Tetrahedron<float, 3> Tetrahedron3f;
typedef Tetrahedron<double, 3> Tetrahedron3d;

template <typename Scalar_, int Dim_>
struct Box {
    using Scalar = Scalar_;
    static constexpr int Dim = Dim_;

    Transform<Scalar, Dim> pose;
    Vector<Scalar, Dim> size;

    BoundingBox<Scalar, Dim> bounding_box() const {
        BoundingBox<Scalar, Dim> result;
        Vector<Scalar, Dim> point;
        for (int i = 0; i < (2<<Dim); i++) {
            for (int dim = 0; dim < Dim; dim++) {
                point(dim) = (i % (1<<dim) == 1 ? Scalar(1) : Scalar(-1)) * size(dim);
            }
            result.merge(pose * point);
        }
        return result;
    }

    void transform(const Transform<Scalar, 3>& transform) {
        pose = transform * pose;
    }

    static Box Instance() {
        Box result;
        result.pose.setIdentity();
        result.size.setOnes();
        return result;
    }

    Matrix<Scalar, Dim+1, Dim+1> instance_transform() const {
        return pose.matrix() * make_homogeneous_scaling<Scalar, Dim>(size);
    }

    template <typename T_>
    Box<T_, Dim> cast() const {
        Box<T_, Dim> result;
        result.pose = pose.template cast<T_>();
        result.size = size.template cast<T_>();
        return result;
    }
};
typedef Box<float, 2> Box2f;
typedef Box<float, 3> Box3f;
typedef Box<double, 2> Box2d;
typedef Box<double, 3> Box3d;

template <typename Scalar_, int Dim_>
struct Sphere {
    using Scalar = Scalar_;
    static constexpr int Dim = Dim_;

    Vector<Scalar, Dim> position;
    Scalar radius;

    BoundingBox<Scalar, Dim> bounding_box() const {
        BoundingBox<Scalar, Dim> result;
        result.lower = position - Eigen::Vector<Scalar, Dim>::Ones() * radius;
        result.upper = position + Eigen::Vector<Scalar, Dim>::Ones() * radius;
        return result;
    }

    void transform(const Transform<Scalar, 3>& transform) {
        position = transform * position;
    }

    static Sphere Instance() {
        Sphere result;
        result.position.setZero();
        result.radius = 0.5;
        return result;
    }

    Matrix<Scalar, Dim+1, Dim+1> instance_transform() const {
        return make_homogeneous_translation(position) * make_homogeneous_scaling<Scalar, Dim>(radius * 2);
    }

    template <typename T_>
    Sphere<T_, Dim> cast() const {
        Sphere<T_, Dim> result;
        result.position = position.template cast<T_>();
        result.radius = static_cast<T_>(radius);
        return result;
    }
};
typedef Sphere<float, 2> Sphere2f;
typedef Sphere<float, 3> Sphere3f;
typedef Sphere<double, 2> Sphere2d;
typedef Sphere<double, 3> Sphere3d;

template <typename Scalar_, int Dim_>
struct Cylinder {
    using Scalar = Scalar_;
    static constexpr int Dim = Dim_;

    Transform<Scalar, Dim> pose;
    Scalar radius; // Radius for dimensions 0 - (dim-2)
    Scalar length; // Length for dimension (dim-1)

    BoundingBox<Scalar, Dim> bounding_box() const {
        BoundingBox<Scalar, Dim> result;

        Vector<Scalar, Dim> top;
        top.setZero();
        top(Dim-1) = length/2;
        top = pose * top;
        Vector<Scalar, Dim> bottom;
        bottom.setZero();
        top(Dim-1) = -length/2;
        bottom = pose * bottom;

        Vector<Scalar, Dim> axis = (top - bottom).normalized();
        for (int dim = 0; dim < Dim; dim++) {
            // A cyinder projected onto a plane is a rectangle * ellipses at either end
            Scalar axis_dim_2 = std::pow(axis(dim), 2);
            Scalar ellipse_dist = 0;
            if (axis_dim_2 > 1e-12) { // arbitrary small number
                Scalar minor_radius_scale_2 = 1 - axis_dim_2;
                Scalar tan_theta_2 = (axis.norm2() - axis_dim_2) / axis_dim_2;
                ellipse_dist = std::sqrt(minor_radius_scale_2) * radius * std::sqrt(
                    (1 + tan_theta_2) / (1 + minor_radius_scale_2 * tan_theta_2)
                );
            }
            result.lower(dim) = bottom(dim) - ellipse_dist;
            result.upper(dim) = top(dim) + ellipse_dist;
        }
        return result;
    }

    void transform(const Transform<Scalar, 3>& transform) {
        pose = transform * pose;
    }

    static Cylinder Instance() {
        Cylinder result;
        result.pose.setIdentity();
        result.radius = 0.5;
        result.length = 1;
        return result;
    }

    Matrix<Scalar, Dim+1, Dim+1> instance_transform() const {
        Vector<Scalar, Dim> scaling;
        scaling.template head<Dim-1>().array() = radius * 2;
        scaling(Dim-1) = length;
        return pose.matrix() * make_homogeneous_scaling<Scalar, Dim>(scaling);
    }

    template <typename T_>
    Cylinder<T_, Dim> cast() const {
        Cylinder<T_, Dim> result;
        result.pose = pose.template cast<T_>();
        result.radius = static_cast<T_>(radius);
        result.length = static_cast<T_>(length);
        return result;
    }
};
typedef Cylinder<float, 3> Cylinder3f;
typedef Cylinder<double, 3> Cylinder3d;

template <typename Scalar_, int Dim_>
struct Cone {
    using Scalar = Scalar_;
    static constexpr int Dim = Dim_;

    Transform<Scalar, Dim> pose;
    Scalar radius;
    Scalar length;

    BoundingBox<Scalar, Dim> bounding_box() const {
        BoundingBox<Scalar, Dim> result;

        // Same as cylinder, but with top radius = 0

        Vector<Scalar, Dim> top;
        top.setZero();
        top(Dim-1) = length/2;
        top = pose * top;
        Vector<Scalar, Dim> bottom;
        bottom.setZero();
        top(Dim-1) = -length/2;
        bottom = pose * bottom;

        Vector<Scalar, Dim> axis = (top - bottom).normalized();
        for (int dim = 0; dim < Dim; dim++) {
            // A cyinder projected onto a plane is a rectangle * ellipses at either end
            Scalar axis_dim_2 = std::pow(axis(dim), 2);
            Scalar ellipse_dist = 0;
            if (axis_dim_2 > 1e-12) { // arbitrary small number
                Scalar minor_radius_scale_2 = 1 - axis_dim_2;
                Scalar tan_theta_2 = (axis.norm2() - axis_dim_2) / axis_dim_2;
                ellipse_dist = std::sqrt(minor_radius_scale_2) * radius * std::sqrt(
                    (1 + tan_theta_2) / (1 + minor_radius_scale_2 * tan_theta_2)
                );
            }
            result.lower(dim) = bottom(dim) - ellipse_dist;
            result.upper(dim) = std::max(top(dim), bottom(dim) + ellipse_dist);
        }
        return result;
    }

    void transform(const Transform<Scalar, 3>& transform) {
        pose = transform * pose;
    }

    static Cone Instance() {
        Cone result;
        result.pose.setIdentity();
        result.radius = 0.5;
        result.length = 1;
        return result;
    }

    Matrix<Scalar, Dim+1, Dim+1> instance_transform() const {
        Vector<Scalar, Dim> scaling;
        scaling.template head<Dim-1>().array() = radius * 2;
        scaling(Dim-1) = length;
        return pose.matrix() * make_homogeneous_scaling<Scalar, Dim>(scaling);
    }

    template <typename T_>
    Cone<T_, Dim> cast() const {
        Cone<T_, Dim> result;
        result.pose = pose.template cast<T_>();
        result.radius = static_cast<T_>(radius);
        result.length = static_cast<T_>(length);
        return result;
    }
};
typedef Cone<float, 2> Cone2f;
typedef Cone<double, 2> Cone2d;
typedef Cone<float, 3> Cone3f;
typedef Cone<double, 3> Cone3d;

template <typename Scalar, int Dim>
using InstancedPrimitive = std::variant<
    Box<Scalar, Dim>,
    Sphere<Scalar, Dim>,
    Cylinder<Scalar, Dim>,
    Cone<Scalar, Dim>>;
typedef InstancedPrimitive<float, 2> InstancedPrimitive2f;
typedef InstancedPrimitive<float, 3> InstancedPrimitive3f;
typedef InstancedPrimitive<double, 2> InstancedPrimitive2d;
typedef InstancedPrimitive<double, 3> InstancedPrimitive3d;

template <typename Scalar, int Dim, typename OtherScalar>
InstancedPrimitive<Scalar, Dim> cast_primitive(const InstancedPrimitive<OtherScalar, Dim>& primitive) {
    return std::visit([](const auto& primitive) -> InstancedPrimitive<Scalar, Dim> {
        return primitive.template cast<Scalar>();
    }, primitive);
}

template <typename Scalar, int Dim>
struct ColoredPrimitive {
    InstancedPrimitive<Scalar, Dim> primitive;
    ColorRGBd color;
    ColoredPrimitive() {}
    ColoredPrimitive(
            const InstancedPrimitive<Scalar, Dim>& primitive,
            const ColorRGBd& color):
        primitive(primitive), color(color)
    {}
};
typedef ColoredPrimitive<float, 2> ColoredPrimitive2f;
typedef ColoredPrimitive<float, 3> ColoredPrimitive3f;
typedef ColoredPrimitive<double, 2> ColoredPrimitive2d;
typedef ColoredPrimitive<double, 3> ColoredPrimitive3d;

} // namespace owl

template <typename Scalar, int Dim>
struct cpp_utils::variant_details<owl::InstancedPrimitive<Scalar, Dim>> {
    typedef owl::InstancedPrimitive<Scalar, Dim> T;
    static constexpr std::size_t count = 4;
    static T construct(std::size_t i) {
        switch(i) {
        case 0:
            return owl::Box<Scalar, Dim>();
        case 1:
            return owl::Sphere<Scalar, Dim>();
        case 2:
            return owl::Cylinder<Scalar, Dim>();
        case 3:
            return owl::Cone<Scalar, Dim>();
        default:
            assert(false);
            return owl::Box<Scalar, Dim>();
        }
    }
    static std::size_t index(const T& variant) {
        if (std::get_if<owl::Box<Scalar, Dim>>(&variant)) {
            return 0;
        }
        else if (std::get_if<owl::Sphere<Scalar, Dim>>(&variant)) {
            return 1;
        }
        else if (std::get_if<owl::Cylinder<Scalar, Dim>>(&variant)) {
            return 2;
        }
        else if (std::get_if<owl::Cone<Scalar, Dim>>(&variant)) {
            return 3;
        }
        else {
            assert(false);
            return 4;
        }
    }
    static std::string_view label(std::size_t i) {
        static const char* labels[count] = {
            "box", "sphere", "cylinder", "cone"
        };
        assert(i < count);
        return labels[i];
    }
};
