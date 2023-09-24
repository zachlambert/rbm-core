#pragma once

#include <random>
#include "core/std/cloneable.h"


namespace math {

template <typename T>
class DistributionImpl {
public:
    virtual T sample(std::default_random_engine& rng) const = 0;
    virtual double evaluate(const T& x) const = 0;
    virtual std::unique_ptr<DistributionImpl<T>> clone() const = 0;
};


template <typename T>
class Distribution: public core::Cloneable<DistributionImpl<T>> {
public:
    template <typename Impl>
    Distribution(const Impl& impl):
        core::Cloneable<DistributionImpl<T>>(impl)
    {}
    T sample(std::default_random_engine& rng) const {
        return this->interface->sample(rng);
    }
    double evaluate(const T& x) const {
        return this->interface->evaluate(x);
    }
};

} // namespace math
