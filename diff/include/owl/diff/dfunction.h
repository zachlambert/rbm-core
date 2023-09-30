#pragma once

#include "owl/diff/manifold.h"
#include "owl/types/matrix.h"
#include "owl/types/mtensor.h"


namespace owl {

// DFunction = Differentiable function
// DFunction1 = Once differentiable
// DFunction2 = Twice differentiable
// Define derivatives using numerator-layout convention. This means:
// For mapping X -> Y, with dimensions n and m:
// - df/dx = mxn matrix
// - d^2f/dx^2 = mxnxn tensor, represented by an array/vector of mxn matrices

// ScalarDFunction[1/2] = Scalar output
// Uses denominator-layout convention since this is more convenient for
// the scalar case.
// Instead of referrring to as df/dx and d^2f/dx^2, this provides the gradient
// and hessian.
// For a mapping X -> R, with X of dimension n:
// - gradient = n vector
// - hessian = nxn matrix

// MultiDFunction[1/2] = Differentiable function with multiple inputs.
// Provides df/dx_1, df/dx_2, ..., d^2f/dx_1^2 d^2f/dx_1dx_2, d^2/dx_2^2, ...
// Same format as standard DFunction

// Each type of function is purely design to provide a wrapper around
// a group of related functions, to pass to another part of code that
// requires them.
// The specific functions f, gradient_f, etc, should be provided as
// simple functions or methods, and only converted to a specific function
// object when required for a given operation.


template <typename X>
class ScalarDFunction1 {
public:
    using f_t = std::function<double(const X& x)>;
    using gradient_t = manifold_delta<X>;
    using gradient_f_t = std::function<gradient_t(const X&)>;

    ScalarDFunction1(const f_t& f, const gradient_f_t& gradient_f):
        f_(f), gradient_f_(gradient_f)
    {}
    double operator()(const X& x) const { return f_(x); }
    gradient_t gradient(const X& x) const { return gradient_f_(x); }

    const f_t& f() const { return f_; }
    const gradient_f_t& gradient_f() const { return gradient_f_; }

private:
    f_t f_;
    gradient_f_t gradient_f_;
};

template <typename X>
class ScalarDFunction2 {
    using f_t = std::function<double(const X& x)>;
    using gradient_t = manifold_delta<X>;
    using gradient_f_t = std::function<gradient_t(const X&)>;
    using hessian_t = manifold_hessian<X>;
    using hessian_f_t = std::function<hessian_t(const X&)>;

public:
    ScalarDFunction2(const f_t& f, const gradient_f_t& gradient_f, const hessian_f_t& hessian_f):
        f_(f), gradient_f_(gradient_f), hessian_f_(hessian_f)
    {}
    double operator()(const X& x) const { return f_(x); }
    gradient_t gradient(const X& x) const { return gradient_f_(x); }
    hessian_t hessian(const X& x) const { return hessian_f_(x); }

    const f_t& f() const { return f_; }
    const gradient_f_t& gradient_f() const { return gradient_f_; }
    const hessian_f_t& hessian_f() const { return hessian_f_; }

private:
    f_t f_;
    gradient_f_t gradient_f_;
    hessian_f_t hessian_f_;
};

template <typename X, typename Y>
class DFunction1 {
public:
    static constexpr int XDim = manifold_dim<X>;
    static constexpr int YDim = manifold_dim<Y>;
    using f_t = std::function<Y(const X& x)>;
    using dfdx_t = Matrixd<YDim, XDim>;
    using dfdx_f_t = std::function<dfdx_t(const X& x)>;

    DFunction1(const f_t& f, const dfdx_f_t& dfdx_f):
        f_(f), dfdx_f_(dfdx_f)
    {}
    Y operator()(const X& x) const { return f_(x); }
    dfdx_t dfdx(const X& x) const { return dfdx_f_(x); }

    const f_t& f() const { return f_; }
    const dfdx_f_t& dfdx_f() const { return dfdx_f_; }

private:
    f_t f_;
    dfdx_f_t dfdx_f_;
};

template <typename X, typename Y>
class DFunction2 {
public:
    static constexpr int XDim = manifold_dim<X>;
    static constexpr int YDim = manifold_dim<Y>;
    using f_t = std::function<Y(const X& x)>;
    using dfdx_t = Matrixd<YDim, XDim>;
    using dfdx_f_t = std::function<dfdx_t(const X& x)>;
    using d2fdx2_t = MTensord<YDim, XDim, XDim>;
    using d2fdx2_f_t = std::function<d2fdx2_t(const X&)>;

    DFunction2(const f_t& f, const dfdx_f_t& dfdx_f, d2fdx2_f_t d2fdx2_f):
        f_(f), dfdx_f_(dfdx_f), d2fdx2_f_(d2fdx2_f)
    {}
    Y operator()(const X& x) const { return f_(x); }
    dfdx_t dfdx(const X& x) const { return dfdx_f_(x); }
    d2fdx2_t d2fdx2(const X& x) const { return d2fdx2_f_(x); }

    const f_t& f() const { return f_; }
    const dfdx_f_t& dfdx_f() const { return dfdx_f_; }
    const d2fdx2_f_t& d2fdx2_f() const { return d2fdx2_f_; }

private:
    f_t f_;
    dfdx_f_t dfdx_f_;
    d2fdx2_f_t d2fdx2_f_;
};

template <typename Y, typename... Xs>
class MultiDFunction1 {
public:
    static constexpr int YDim = manifold_dim<Y>;
    using f_t = std::function<Y(const Xs&...)>;
    using dfdx_fs_t = std::tuple<
        std::function<
            Matrixd<YDim, manifold_dim<Xs>>
            (const Xs&...)>
        ...>;

    MultiDFunction1(const f_t& f, const dfdx_fs_t& dfdx_fs):
        f_(f), dfdx_fs_(dfdx_fs)
    {}
    Y operator()(const Xs&... xs) const {
        return f(xs...);
    }

    // TODO

    // template <int Index>
    // Matrixd<YDim, manifold_dim<std::tuple_element_t<std::tuple<Xs...>, Index>>> dfdx(const Xs&... xs) const {
    //     return std::get<Index>(dfdx_fs_)(xs...);
    // }

private:
    f_t f_;
    dfdx_fs_t dfdx_fs_;
};

// MultiDFunction2: TODO


} // namespace math
