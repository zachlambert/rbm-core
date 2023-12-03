#pragma once

#include <optional>


namespace rbm {

// Represents the interval [lower, upper]
// (ie: includes endpoints)
template <typename T>
struct Interval {
    std::optional<T> lower;
    std::optional<T> upper;
    void clamp(T& value) const {
        if (lower.has_value() && value < lower.value()) {
            value = lower.value();
        }
        if (upper.has_value() && value > upper.value()) {
            value = upper.value();
        }
    }
    bool contains(const T& value) const {
        return
            (!lower.has_value() || lower.value() <= value)
            && (!upper.has_value() || value <= upper.value());
    }
    std::optional<T> center() const {
        if (lower.has_value() && upper.has_value()) {
            return 0.5 * (lower.value() + upper.value());
        }
        return std::nullopt;
    }
    T error(const T& value) const {
        if (lower.has_value() && value < lower.value()) {
            return value - lower.value();
        }
        if (upper.has_value() && value > upper.value()) {
            return value - upper.value();
        }
        return 0;
    }
    T error_derivative(const T& value, const T& velocity) const {
        if (contains(value)) {
            return 0;
        }
        return velocity;
    }

    bool closed() const {
        return lower.has_value() && upper.has_value();
    }
    Interval closed_or(const Interval& other) const {
        if (closed()) return *this;
        return other;
    }

    static Interval Open() {
        return Interval();
    }
    static Interval Closed(const T& lower, const T& upper) {
        Interval interval;
        interval.lower = lower;
        interval.upper = upper;
        return interval;
    }
    static Interval LowerBound(const T& lower) {
        Interval interval;
        interval.lower = lower;
        return interval;
    }
    static Interval UpperBound(const T& upper) {
        Interval interval;
        interval.upper = upper;
        return interval;
    }
};
typedef Interval<float> Intervalf;
typedef Interval<double> Intervald;

} // namespace rbm
