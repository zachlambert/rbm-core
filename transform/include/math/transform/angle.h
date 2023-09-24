#pragma once

#include <cmath>
#include "math/types/matrix.h"


namespace math {

template <typename T>
void normalize_angle(T& angle) {
    static constexpr T two_pi = 2 * M_PI;
    static constexpr T pi = M_PI;
    angle -= two_pi * std::floor((angle + pi) / two_pi);
}

template <typename T>
T normalized_angle(T angle) {
    normalize_angle(angle);
    return angle;
}

template <typename T>
class Angle {
public:
    Angle(): angle(0) {}
    Angle(T angle): angle(angle) {
        normalize();
    }
    operator const T&()const {
        return angle;
    }

    Angle<T>& operator+=(const Angle<T>& rhs) {
        angle += rhs.angle;
        normalize();
        return *this;
    }
    Angle<T>& operator-=(const Angle<T>& rhs) {
        angle -= rhs.angle;
        normalize();
        return *this;
    }
    Angle<T>& operator*=(T rhs) {
        angle *= rhs;
        normalize();
        return *this;
    }
    Angle<T>& operator/=(T rhs) {
        angle /= rhs.angle;
        normalize();
        return *this;
    }

    Angle<T> inverse() const {
        return Angle<T>(-angle);
    }

    static Angle<T> Zero() {
        return Angle<T>(0);
    }

    template <typename T_>
    Angle<T_> cast() const {
        return Angle<T_>(static_cast<T_>(angle));
    }

    Angle<T>& operator=(const Eigen::Matrix<T, 2, 2>& rotation) {
        angle = std::atan2(rotation(1, 0), rotation(0, 0));
        return angle;
    }

    Matrix<T, 2, 2> toRotationMatrix() const {
        Matrix<T, 2, 2> R;
        R << std::cos(angle), -std::sin(angle), std::sin(angle), std::cos(angle);
        return R;
    }

private:
    void normalize() {
        normalize_angle(angle);
    }
    T angle;
};
typedef Angle<float> Anglef;
typedef Angle<double> Angled;

template <typename T>
Angle<T> operator+(Angle<T> lhs, const Angle<T>& rhs)
{
    lhs += rhs;
    return lhs;
}
template <typename T>
Angle<T> operator-(Angle<T> lhs, const Angle<T>& rhs)
{
    lhs -= rhs;
    return lhs;
}
template <typename T>
Angle<T> operator*(Angle<T> lhs, T rhs)
{
    lhs *= rhs;
    return lhs;
}
template <typename T>
Angle<T> operator/(Angle<T> lhs, T rhs)
{
    lhs /= rhs;
    return lhs;
}
template <typename T>
Angle<T> operator-(Angle<T> lhs)
{
    lhs = -static_cast<T>(lhs);
    return lhs;
}

template <typename T>
Eigen::Vector<T, 2> operator*(const Angle<T>& lhs, const Eigen::Vector<T, 2>& rhs) {
    return Eigen::Rotation2D<T>(lhs) * rhs;
}

} // namespace math
