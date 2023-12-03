#pragma once

#include "rbm/transform/transform.h"
#include "rbm/geometry/primitive.h"
#include "rbm/transform/adjoint.h"
#include <Eigen/Jacobi>


namespace rbm {

template <typename Scalar, int Dim>
struct SpatialInertia {
    typedef Matrix<Scalar, dim_se<Dim>(), dim_se<Dim>()> matrix_t;

    double mass;
    matrix_t gyration;

    SpatialInertia()
    {
        setZero();
    }


    matrix_t matrix() const {
        return gyration * mass;
    }

    void operator+=(const SpatialInertia& rhs) {
        double new_mass = mass + rhs.mass;
        gyration = gyration * (mass/new_mass) + rhs.gyration * (rhs.mass/new_mass);
        mass = new_mass;
    }

    friend SpatialInertia operator+(SpatialInertia lhs, const SpatialInertia& rhs) {
        lhs += rhs;
        return lhs;
    }

    template <typename Scalar_>
    SpatialInertia<Scalar_, Dim> cast() const {
        SpatialInertia<Scalar_, Dim> result;
        result.mass = static_cast<Scalar_>(mass);
        result.gyration = gyration.template cast<Scalar_>();
        return result;
    }

    void setZero() {
        mass = 0;
        gyration.setZero();
        gyration.template block<3, 3>(3, 3).setIdentity();
    }

    static SpatialInertia Zero() {
        SpatialInertia result;
        result.setZero();
        return result;
    }

    void transform(const Transform<Scalar, Dim>& transform)
    {
        gyration =
            AdjointTransform(transform).matrix().transpose()
            * gyration
            * AdjointTransform(transform).matrix();
    }

    Vector3d com() const {
        return -Vector3d(
            -gyration(4, 2),
            gyration(3, 2),
            -gyration(3, 1)
        );
    }

    Transform3d principle_frame() const {
        SpatialInertia temp = *this;
        temp.transform(Transform3d(Rotation3d::Identity(), temp.com()));

        Eigen::JacobiSVD<Matrix3d> svd;
        svd.compute(temp.gyration.template block<3, 3>(0, 0), Eigen::ComputeFullU | Eigen::ComputeFullV);
        Matrix3d principle_axes = svd.matrixV();
        // Symmetric matrix, so singular vectors define an orthonormal basis
        // satisfies |det(V)| = 1
        // To be a valid rotation, require det(V) = +1
        // Therefore, mirror one of the columns if det(V) = -1
        double det = principle_axes.determinant();
        if (det < 0) {
            principle_axes.block<3, 1>(0, 0) *= -1;
        }

        return Transform3d(principle_axes, com());
    }

};
typedef SpatialInertia<float, 2> SpatialInertia2f;
typedef SpatialInertia<double, 2> SpatialInertia2d;
typedef SpatialInertia<float, 3> SpatialInertia3f;
typedef SpatialInertia<double, 3> SpatialInertia3d;

template <typename Scalar, int Dim>
SpatialInertia<Scalar, Dim> primitive_to_inertia(const Box<Scalar, Dim>& box, double density) {
    static_assert(Dim == 2 || Dim == 3);
    SpatialInertia<Scalar, Dim> result;

    double volume = 1;
    for (int dim = 0; dim < Dim; dim++) volume *= box.size(dim);
    result.mass = volume * density;

    if (Dim == 2) {
        result.gyration(0, 0) = std::pow(box.size(0), 2) + std::pow(box.size(1), 2);
    }
    else {
        for (int dim = 0; dim < Dim; dim++) {
            double moment = 0;
            for (int other = 0; other < Dim; other++) {
                if (dim == other) continue;
                moment += std::pow(box.size(other), 2);
            }
            result.gyration(dim, dim) = moment;
        }
    }

    result.transform(box.pose.inverse());

    return result;
}

template <typename Scalar, int Dim>
SpatialInertia<Scalar, Dim> primitive_to_inertia(const Sphere<Scalar, Dim>& sphere, double density) {
    static_assert(Dim == 2 || Dim == 3);
    SpatialInertia<Scalar, Dim> result;

    double volume;
    if (Dim == 2) {
        volume = M_PI * std::pow(sphere.radius, 2);
    }
    else {
        volume = (4.0 / 3) * M_PI * std::pow(sphere.radius, 3);
    }
    result.mass = volume * density;

    if (Dim == 2) {
        result.gyration(0, 0) = 0.5 * std::pow(sphere.radius, 2);
    }
    else {
        result.gyration.template block<3, 3>(0, 0) = Matrix<Scalar, dim_so<Dim>(), dim_so<Dim>()>::Identity() * 2.0 / 5 * std::pow(sphere.radius, 2);
    }

    Transform<Scalar, Dim> offset = Transform<Scalar, Dim>::Identity();
    offset.translation() = sphere.position;
    result.transform(offset.inverse());

    return result;
}

template <typename Scalar, int Dim>
SpatialInertia<Scalar, Dim> primitive_to_inertia(const Cylinder<Scalar, Dim>& cylinder, double density) {
    static_assert(Dim == 3);
    SpatialInertia<Scalar, Dim> result;

    double volume = cylinder.length * M_PI * std::pow(cylinder.radius, 2);
    result.mass = volume * density;

    for (int dim = 0; dim < Dim; dim++) {
        if (dim == Dim-1) {
            result.gyration(dim, dim) = std::pow(cylinder.radius, 2) / 2;
        }
        else {
            result.gyration(dim, dim) = std::pow(cylinder.radius, 2) / 4 + std::pow(cylinder.length, 2) / 12;
        }
    }

    result.transform(cylinder.pose.inverse());

    return result;
}

template <typename Scalar, int Dim>
SpatialInertia<Scalar, Dim> primitive_to_inertia(const Cone<Scalar, Dim>& cone, double density) {
    static_assert(Dim == 3);
    SpatialInertia<Scalar, Dim> result;

    double volume = (1.0 / 3) * cone.length * M_PI * std::pow(cone.radius, 2);
    result.mass = volume * density;

    for (int dim = 0; dim < Dim; dim++) {
        if (dim == Dim-1) {
            result.gyration(dim, dim) = 0.3 * std::pow(cone.radius, 2);
        }
        else {
            result.gyration(dim, dim) = 3 * (4 * std::pow(cone.radius, 2) + std::pow(cone.length, 2)) / 80;
        }
    }

    Transform3d com_offset = Transform3d::Identity();
    com_offset.translation() = Vector3d::UnitZ() * cone.length / 4;
    result.transform((cone.pose * com_offset).inverse());

    return result;
}

template <typename Scalar, int Dim>
SpatialInertia<Scalar, Dim> primitive_to_inertia(const InstancedPrimitive<Scalar, Dim>& primitive, double density) {
    return std::visit([density](const auto& value) { return primitive_to_inertia(value, density); }, primitive);
}

} // namespace rbm
