#pragma once

#include "rbm/types/matrix.h"
#include "rbm/types/interval.h"
#include <numeric>

namespace rbm {

template <typename Scalar, int Dim>
struct BoundingBox {
    Vector<Scalar, Dim> lower;
    Vector<Scalar, Dim> upper;

    BoundingBox():
        lower(Vector<Scalar, Dim>::Ones() * std::numeric_limits<Scalar>::max()),
        upper(-Vector<Scalar, Dim>::Ones() * std::numeric_limits<Scalar>::max())
    {}
    Vector<Scalar, Dim> centre() const {
        return (lower + upper) / 2;
    }
    Vector<Scalar, Dim> size() const {
        return upper - lower;
    }
    Interval<Scalar> interval(int dim) const {
        return Interval<Scalar>(lower(dim), upper(dim));
    }
    void merge(const Vector<Scalar, Dim>& point) {
        for (int dim = 0; dim < Dim; ++dim) {
            lower(dim) = std::min(point(dim), lower(dim));
            upper(dim) = std::max(point(dim), upper(dim));
        }
    }
    void merge(const BoundingBox& box) {
        for (int dim = 0; dim < Dim; ++dim) {
            lower(dim) = std::min(box.lower(dim), lower(dim));
            upper(dim) = std::max(box.upper(dim), upper(dim));
        }
    }
    void split(int dim, Scalar partition, BoundingBox& lower, BoundingBox& upper) {
        lower = *this;
        upper = *this;
        lower.upper(dim) = partition;
        upper.lower(dim) = partition;
    }

    bool intersects(const BoundingBox<Scalar, Dim>& other) const {
        for (int dim = 0; dim < Dim; dim++) {
            if (lower(dim) > other.upper(dim)) return false;
            if (upper(dim) < other.lower(dim)) return false;
        }
        return true;
    }

    template <typename T_>
    BoundingBox<T_, Dim> cast() const {
        BoundingBox<T_, Dim> result;
        result.lower = lower.template cast<T_>();
        result.upper = upper.template cast<T_>();
        return result;
    }
};
typedef BoundingBox<float, 2> BoundingBox2f;
typedef BoundingBox<double, 2> BoundingBox2d;
typedef BoundingBox<float, 3> BoundingBox3f;
typedef BoundingBox<double, 3> BoundingBox3d;

} // namespace rbm
