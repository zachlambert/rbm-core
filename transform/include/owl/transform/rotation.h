#pragma once

#include <Eigen/Geometry>
#include "owl/transform/group_dimensions.h"
#include "owl/transform/angle.h"
#include "owl/transform/cross_product_matrix.h"
#include "owl/transform/screw_transform.h"
#include <cmath>


namespace owl {

template <typename Scalar, int Dim>
class Rotation: public Eigen::Matrix<Scalar, Dim, Dim> {
public:
    Rotation() {
        (*this).setIdentity();
    }

    // Copied from:
    // https://eigen.tuxfamily.org/dox/TopicCustomizing_InheritingMatrix.html

    template<typename OtherDerived>
    Rotation(const Eigen::MatrixBase<OtherDerived>& other):
        Eigen::Matrix<Scalar, Dim, Dim>(other)
    {}

    template<typename OtherDerived>
    Rotation& operator=(const Eigen::MatrixBase<OtherDerived>& other) {
        this->Eigen::Matrix<Scalar, Dim, Dim>::operator=(other);
        return *this;
    }

    static Rotation Identity() {
        return Rotation();
    }
};
typedef Rotation<float, 2> Rotation2f;
typedef Rotation<float, 3> Rotation3f;
typedef Rotation<double, 2> Rotation2d;
typedef Rotation<double, 3> Rotation3d;
template <int Dim>
using Rotationf = Rotation<float, Dim>;
template <int Dim>
using Rotationd = Rotation<double, Dim>;


template <typename Scalar, int Dim>
using CompactRotation = std::conditional_t<Dim==2, Angle<Scalar>, Eigen::Quaternion<Scalar>>;
typedef CompactRotation<float, 2> CompactRotation2f;
typedef CompactRotation<double, 2> CompactRotation2d;
typedef CompactRotation<float, 3> CompactRotation3f;
typedef CompactRotation<double, 3> CompactRotation3d;


template <typename Scalar, int Dim>
class LogRotation {
public:
    static_assert(Dim == 2 || Dim == 3);
    typedef Vectord<dim_so<Dim>()> coords_t;
    typedef Eigen::Matrix<Scalar, Dim, Dim> matrix_t;

    LogRotation(const coords_t& coords):
        coords_(coords)
    {}

    LogRotation(const Rotation<Scalar, Dim>& rotation) {
        if constexpr (Dim == 2) {
            coords_ = std::atan2(rotation(1, 0), rotation(0,0));
        }
        if constexpr (Dim == 3) {
            // Note: With an identity rotation, 0.5*(trace - 1) = 1, but can be very slightly above 1
            // due to floating point errors, so need to clamp it
            Scalar angle = std::acos(std::clamp<double>(0.5*(rotation.trace() - 1), -1, 1));
            coords_ = 0.5 * sinc(angle) * Eigen::Matrix<Scalar, 3, 1>(
                rotation(2, 1) - rotation(1, 2),
                rotation(0, 2) - rotation(2, 0),
                rotation(1, 0) - rotation(0, 1)
            );
        }
    }

    LogRotation(const CompactRotation<Scalar, Dim>& rotation) {
        if constexpr (Dim == 2) {
            coords_ = rotation;
        }
        if constexpr (Dim == 3) {
            Eigen::AngleAxis<Scalar> angle_axis;
            angle_axis = rotation;
            coords_ = angle_axis.angle() * angle_axis.axis();
        }
    }

    coords_t& coords() { return coords_; }
    const coords_t& coords()const { return coords_; }
    matrix_t matrix() { return cross_product_matrix(coords_); }

    Rotation<Scalar, Dim> exp() const {
        if constexpr (Dim == 2) {
            return Angle<Scalar>(coords_.value()).toRotationMatrix();
        }
        if constexpr (Dim == 3) {
            Eigen::Vector<Scalar, 3> axis = coords_.normalized();
            Scalar angle = axis.dot(coords_);
            return Eigen::AngleAxis<Scalar>(angle, axis).toRotationMatrix();
        }
    }
    CompactRotation<Scalar, Dim> exp_compact() const {
        if constexpr (Dim == 2) {
            return coords_;
        }
        if constexpr (Dim == 3) {
            Eigen::Vector<Scalar, 3> axis = coords_.normalized();
            Scalar angle = axis.dot(coords_);
            Eigen::Quaternion<Scalar> quat;
            quat = Eigen::AngleAxis<Scalar>(angle, axis);
            return quat;
        }
    }
private:
    coords_t coords_;
};
typedef LogRotation<float, 2> LogRotation2f;
typedef LogRotation<double, 2> LogRotation2d;
typedef LogRotation<float, 3> LogRotation3f;
typedef LogRotation<double, 3> LogRotation3d;


template <typename Scalar>
Eigen::Quaternion<Scalar> rotation_matrix_to_quaternion(Eigen::Matrix<Scalar, 3, 3>& rotation)
{
    const auto& R = rotation; // Shorter name
    Eigen::Quaternion<Scalar> q;

    q.w() = 0.5*std::sqrt(1 + R.matrix()(0, 0) + R.matrix()(1, 1) + R.matrix()(2, 2));
    if (q.w() == 0) {
        q.x() = std::sqrt(0.5*(1+R.matrix()(0, 0)));
        q.y() = std::sqrt(0.5*(1+R.matrix()(1, 1)));
        q.z() = std::sqrt(0.5*(1+R.matrix()(2, 2)));
        return q;
    }
    q.x() = 0.25 * (R.matrix()(2, 1) - R.matrix()(1, 2)) / q.w();
    q.y() = 0.25 * (R.matrix()(0, 2) - R.matrix()(2, 0)) / q.w();
    q.z() = 0.25 * (R.matrix()(1, 0) - R.matrix()(0, 1)) / q.w();
    return q;
}

} // namespace owl
