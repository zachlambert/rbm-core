#pragma once

#include "owl/types/matrix.h"
#include "owl/transform/transform.h"
#include "owl/mechanics/inertia.h"


namespace owl {

static constexpr double gravity = 9.81;

inline Vector3d gravity_3d()
{
    return Vector3d(0, 0, -gravity);
}

inline Vector6d gravity_6d()
{
    Vector6d g = Vector6d::Zero();
    g(5) = -gravity;
    return g;
}

inline Vector6d frame_gravity_force(const Rotation3d& rotation, const SpatialInertia3d& inertia)
{
    // T = link -> gravity_aligned = Translate(com) * (Link Rotation)^T
    // Force mapped with:
    // F_link = AdointTransform(T)^{-1}^T * F_gravity_aligned

    Transform3d transform(rotation.transpose(), inertia.com());
    return AdjointTransform3d(transform).inverse().matrix().transpose() * gravity_6d() * inertia.mass;
}


} // namespace owl
