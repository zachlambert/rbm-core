#pragma once

#include "mbox/probability//gaussian.h"

namespace mbox {

// Weiner process = Integral of white noise (continuous)
// Value at time t is zero mean with variance proportional to t

template <int Dim>
Vectord<Dim> white_noise_integral(const Matrixd<Dim, Dim>& dfdw, double dt, std::default_random_engine& rng)
{
    GaussianDistribution<Vectord<Dim>> distribution;
    if constexpr(Dim == -1) {
        distribution.set_mean(Vectord<Dim>::Zero(dfdw.rows()));
    }
    distribution.set_covariance(dfdw * dfdw.transpose() * dt);
    return distribution.sample(rng);
}


} // namespace mbox
