#pragma once

#include "rbm/types/matrix.h"
#include "rbm/transform/transform.h"

namespace rbm {

template <typename Scalar, int From, int To>
Matrix<Scalar, To, From> velocity_mapping()
{
    static_assert(From == 2 || From == 3 || From == 6);
    static_assert(To == 2 || To == 3 || To == 6);

    if constexpr (From == To) {
        return Matrix<Scalar, To, From>::Identity();
    }
    if constexpr (From == 2 && To == 3) {
        Matrix<Scalar, 3, 2> result;
        result <<
            0, 1,
            1, 0,
            0, 0;
        return result;
    }
    if constexpr (From == 2 && To == 6) {
        Matrix<Scalar, 6, 2> result;
        result <<
            0, 0,
            0, 0,
            0, 1,
            1, 0,
            0, 0,
            0, 0;
        return result;
    }
    if constexpr (From == 3 && To == 2) {
        Matrix<Scalar, 2, 3> result;
        result <<
            0, 1, 0,
            1, 0, 0;
        return result;
    }
    if constexpr (From == 3 && To == 6) {
        Matrix<Scalar, 6, 3> result;
        result <<
            0, 0, 0,
            0, 0, 0,
            1, 0, 0,
            0, 1, 0,
            0, 0, 1,
            0, 0, 0;
        return result;
    }
    if constexpr (From == 6 && To == 2) {
        Matrix<Scalar, 2, 6> result;
        result <<
            0, 0, 0, 1, 0, 0,
            0, 0, 1, 0, 0, 0;
        return result;
    }
    if constexpr (From == 6 && To == 3) {
        Matrix<Scalar, 3, 6> result;
        result <<
            0, 0, 1, 0, 0, 0,
            0, 0, 0, 1, 0, 0,
            0, 0, 0, 0, 1, 0;
        return result;
    }
}

template <typename Scalar>
Matrix<Scalar, 3, 2> spatial_2d_nonholonomic_to_holonomic()
{
    Matrix<Scalar, 3, 2> result;
    result.setZero();
    result(0, 1) = 1;
    result(1, 0) = 1;
    return result;
}

template <typename Scalar>
Matrix<Scalar, 2, 3> spatial_2d_holonomic_to_nonholonomic()
{
    Matrix<Scalar, 2, 3> result;
    result.setZero();
    result(0, 1) = 1;
    result(1, 0) = 1;
    return result;
}

template <typename Scalar>
Matrix<Scalar, 6, 3> spatial_2d_to_3d()
{
    Matrix<Scalar, 6, 3> result;
    result.setZero();
    result.template block<3, 3>(2, 0).setIdentity();
    return result;
}

template <typename Scalar>
Matrix<Scalar, 3, 6> spatial_3d_to_2d()
{
    Matrix<Scalar, 3, 6> result;
    result.setZero();
    result.template block<3, 3>(0, 2).setIdentity();
    return result;
}

template <typename Scalar>
Matrix3<Scalar> normalise_orientation(Matrix3<Scalar> rotation, const Vector3<Scalar>& target_normal)
{
    typedef Vector3<Scalar> Vector;
    const Vector normal = rotation.template block<3, 1>(0, 2);
    const Vector cp = normal.cross(target_normal);
    Scalar sin_angle = cp.norm();
    Scalar angle = std::asin(sin_angle);
    return Eigen::AngleAxis<Scalar>(angle, cp.normalized()) * rotation;
}

} // namespace rbm
