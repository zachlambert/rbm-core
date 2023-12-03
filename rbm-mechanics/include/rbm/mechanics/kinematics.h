#pragma once

#include "rbm/transform/transform.h"
#include "rbm/transform/adjoint.h"

namespace rbm {

template <typename Scalar, int Dim>
Vector<Scalar, dim_se_adj<Dim>()> frame_velocity(
    const Vector<Scalar, dim_se_adj<Dim>()>& parent_velocity,
    const Transform<Scalar, Dim>& relative_transform,
    const Vector<Scalar, dim_se_adj<Dim>()>& relative_velocity)
{
    return
        AdjointTransform<Scalar, Dim>(relative_transform).inverse().matrix() * parent_velocity
        + relative_velocity;
}

template <typename Scalar, int Dim>
Vector<Scalar, dim_se_adj<Dim>()> frame_acceleration(
    const Vector<Scalar, dim_se_adj<Dim>()>& frame_velocity,
    const Vector<Scalar, dim_se_adj<Dim>()>& parent_acceleration,
    const Transform<Scalar, Dim>& relative_transform,
    const Vector<Scalar, dim_se_adj<Dim>()>& relative_velocity,
    const Vector<Scalar, dim_se_adj<Dim>()>& relative_acceleration)
{
    return
        AdjointTransform<Scalar, Dim>(relative_transform).inverse().matrix() * parent_acceleration
        + AdjointLieBracket<Scalar, Dim>(frame_velocity).matrix() * relative_velocity
        + relative_acceleration;
}

} // namespace rbm
