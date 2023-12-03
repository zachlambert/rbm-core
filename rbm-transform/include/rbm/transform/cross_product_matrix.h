#pragma once

#include "rbm/types/matrix.h"


namespace rbm {

template <typename Scalar>
inline Eigen::Matrix<Scalar, 2, 2> cross_product_matrix(Scalar scalar)
{
    Eigen::Matrix<Scalar, 2, 2> result;
    result << static_cast<Scalar>(0), -scalar, scalar, static_cast<Scalar>(0);
    return result;
}

template <typename Scalar>
inline Eigen::Matrix<Scalar, 2, 2> cross_product_matrix(Eigen::Vector<Scalar, 1> scalar)
{
    Eigen::Matrix<Scalar, 2, 2> result;
    result << static_cast<Scalar>(0), -scalar, scalar, static_cast<Scalar>(0);
    return result;
}

template <typename Scalar>
inline Eigen::Matrix<Scalar, 2, 1> cross_product_matrix(const Eigen::Vector<Scalar, 2>& vector)
{
    return Eigen::Matrix<Scalar, 2, 1>(-vector.y(), vector.x());
}

template <typename Scalar>
inline Eigen::Matrix<Scalar, 3, 3> cross_product_matrix(const Eigen::Vector<Scalar, 3>& vector)
{
    Eigen::Matrix<Scalar, 3, 3> result;
    result << 0, -vector.z(), vector.y(),
              vector.z(), 0, -vector.x(),
              -vector.y(), vector.x(), 0;
    return result;
}

} // namespace rbm
