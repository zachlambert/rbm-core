
#include "math/geometry/height_map_functional.h"

// Note: Must include this eigen header, otherwise there are linker errorV
// with the corss product
#include <Eigen/Geometry>


namespace math {

double HeightMapGaussianMixture::height(const math::Vector2d& position) const {
    double z = 0;
    for (const auto& component: components_) {
        z += component.z(position);
    }
    return z;
}
math::Vector3d HeightMapGaussianMixture::normal(const math::Vector2d& position) const {
    math::Vector2d gradient = math::Vector2d::Zero();
    for (const auto& component: components_) {
        gradient += component.gradient(position);
    }

    math::Vector3d forward;
    forward.head<2>() = gradient.normalized() * std::sqrt(1 - gradient.squaredNorm());
    forward(2) = gradient.norm();

    math::Vector3d normal = forward.cross((math::Vector3d::UnitZ().cross(forward)).normalized());
    return normal;
}

} // namespace math
