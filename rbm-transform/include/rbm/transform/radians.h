#pragma once

#include <cmath>


namespace rbm {

template <typename T>
T to_radians(T degrees) {
    return degrees * M_PI / 180;
}

template <typename T>
T to_degrees(T radians) {
    return radians * 180 / M_PI;
}

} // namespace rbm
