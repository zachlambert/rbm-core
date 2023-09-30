
#include "owl/geometry/height_map_functional.h"

// Note: Must include this eigen header, otherwise there are linker error
// with the cross product
#include <Eigen/Geometry>


namespace owl {

double HeightMapGaussianMixture::height(const Vector2d& position) const {
    double z = 0;
    for (const auto& component: components_) {
        z += component.z(position);
    }
    return z;
}
Vector3d HeightMapGaussianMixture::normal(const Vector2d& position) const {
    Vector2d gradient = Vector2d::Zero();
    for (const auto& component: components_) {
        gradient += component.gradient(position);
    }

    Vector3d forward;
    forward.head<2>() = gradient.normalized() * std::sqrt(1 - gradient.squaredNorm());
    forward(2) = gradient.norm();

    Vector3d normal = forward.cross((Vector3d::UnitZ().cross(forward)).normalized());
    return normal;
}

} // namespace owl
