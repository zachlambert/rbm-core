#pragma once

#include <cmath>
#include "math/types/matrix.h"


namespace math {

template <typename T>
T sinc(T x)
{
    static constexpr T epsilon = 1e-6;
    if (std::fabs(x) < epsilon) {
        // 2nd order taylor expansion about x=0
        return 1 - (2.0/3) * std::pow(x, 2);
    } else {
        return std::sin(x) / x;
    }
}

template <typename Scalar, int Dim>
Matrix<Scalar, Dim, Dim> screw_displacement_matrix(
    const Eigen::Matrix<Scalar, dim_so<Dim>(), 1>& rot_coords)
{
    Scalar angle = rot_coords.norm();
    auto S = cross_product_matrix(rot_coords);

    // = (1 - cos(theta)) / theta^2
    Scalar cos_comp;
    if (std::fabs(angle) < 1e-4){
        cos_comp = 0.5 - (1.0/12) * std::pow(angle, 2);
    } else {
        cos_comp = (1 - std::cos(angle)) / std::pow(angle, 2);
    }
    // = (1 - sinc(theta)) / theta^2
    Scalar sin_comp;
    if (std::fabs(angle) < 1e-4) {
        sin_comp =  (1.0/6) - (1.0/120) * std::pow(angle, 2);
    } else {
        sin_comp = (1 - sinc(angle)) / std::pow(angle, 2);
    }

    return decltype(S)::Identity() + cos_comp * S + sin_comp * S * S;
}

} // namespace math
