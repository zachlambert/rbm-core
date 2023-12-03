#pragma once

#include "rbm/types/matrix.h"


namespace rbm {

template <typename Scalar>
Scalar angle_between(const Vector<Scalar, 3>& a, const Vector<Scalar, 3>& b)
{
    typedef ::rbm::Vector<Scalar, 3> Vector;

    if (a.norm() == 0 || b.norm() == 0) return 0;

    Vector normal = a.cross(b);
    if (normal.norm() == 0) {
        if (a.dot(b) > 0) {
            return 0;
        } else {
            return M_PI;
        }
    }

    Vector u1 = a.normalized();
    Vector u2 = normal.cross(a).normalized();
    return atan2(u2.dot(b), u1.dot(b));
}

// Returns a rotation where z axis is the normal given, and the other
// axes are arbitrary perpendicular vectors on the plane defined by the normal.
template <typename Scalar>
Matrix<Scalar, 3, 3> rotation_from_direction(const Vector<Scalar, 3>& direction, int coordinate) {
    // n x n = 0
    // S(n)n = 0
    // [ s1 // s2 // s3 ]n = 0
    // -> Rows of cross product matrix give candidate normal vectors
    // one or two of these may be zero, so pick the one with the
    // largest norm
    Matrix<Scalar, 3, 3> S = cross_product_matrix(direction);
    Scalar max = -1;
    std::size_t max_i = 0;
    for (std::size_t i = 0; i < 3; i++) {
        Scalar norm = S.template block<1, 3>(i, 0).norm();
        if (norm > max) {
            max_i = i;
            max = norm;
        }
    }
    Matrix<Scalar, 3, 3> R;
    Vector<Scalar, 3> perp1 = S.template block<1, 3>(max_i, 0).transpose().normalized();
    Vector<Scalar, 3> perp2 = direction.cross(perp1);

    R.template block<3, 1>(0, coordinate) = direction;
    R.template block<3, 1>(0, (coordinate+1)%3) = perp1;
    R.template block<3, 1>(0, (coordinate+2)%3) = perp2;
    return R;
}

template <typename Scalar>
Matrix<Scalar, 3, 3> rotation_from_normal(const Vector<Scalar, 3>& normal) {
    return rotation_from_direction(normal, 2);
}

} // namespace rbm
