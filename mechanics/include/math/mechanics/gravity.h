#pragma once

#include "math/types/matrix.h"
#include "math/transform/transform.h"
#include "math/mechanics/inertia.h"


namespace math {

static constexpr double gravity = 9.81;

inline math::Vector3d gravity_3d()
{
    return math::Vector3d(0, 0, -gravity);
}

inline math::Vector6d gravity_6d()
{
    math::Vector6d g = math::Vector6d::Zero();
    g(5) = -gravity;
    return g;
}

inline math::Vector6d frame_gravity_force(const math::Rotation3d& rotation, const math::SpatialInertia3d& inertia)
{
    // T = link -> gravity_aligned = Translate(com) * (Link Rotation)^T
    // Force mapped with:
    // F_link = AdointTransform(T)^{-1}^T * F_gravity_aligned

    math::Transform3d transform(rotation.transpose(), inertia.com());
    return AdjointTransform3d(transform).inverse().matrix().transpose() * gravity_6d() * inertia.mass;
}


} // namespace math
