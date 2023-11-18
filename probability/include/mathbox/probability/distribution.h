#pragma once

#include <random>
#include "cpp_utils/cloneable.h"


namespace owl {

template <typename T>
class DistributionImpl {
public:
    virtual T sample(std::default_random_engine& rng) const = 0;
    virtual double evaluate(const T& x) const = 0;
    virtual std::unique_ptr<DistributionImpl<T>> clone() const = 0;
};


template <typename T>
class Distribution: public cpp_utils::Cloneable<DistributionImpl<T>> {
public:
    template <typename Impl>
    Distribution(const Impl& impl):
        cpp_utils::Cloneable<DistributionImpl<T>>(impl)
    {}
    T sample(std::default_random_engine& rng) const {
        return this->interface->sample(rng);
    }
    double evaluate(const T& x) const {
        return this->interface->evaluate(x);
    }
};

} // namespace owl
