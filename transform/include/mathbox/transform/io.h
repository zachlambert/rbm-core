#pragma once

#include "mathbox/transform/euler.h"
#include <iostream>


template <typename Scalar, int Dim>
std::ostream& operator<<(std::ostream& os, const mbox::EulerRotation<Scalar, Dim>& rotation)
{
    std::cout << "[ ";
    if constexpr (Dim == 3) {
        std::cout << "roll: " << rotation.roll() << ", ";
        std::cout << "pitch: " << rotation.pitch() << ", ";
    }
    std::cout << "yaw: " << rotation.yaw() << " ]";
    return os;
}

template <typename Scalar, int Dim>
std::ostream& operator<<(std::ostream& os, const mbox::EulerTransform<Scalar, Dim>& transform)
{
    std::cout << "rot: " << transform.rotation() << ", ";
    std::cout << "trans: " << transform.translation().transpose();
    return os;
}
