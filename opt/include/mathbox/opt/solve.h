#pragma once

#include <concepts>
#include <vector>


// In this library, all optimisation is the following:
// x^* = argmin_x[ f(x) ]
// OR
// x^* = argmin_x[ 0.5 * |W * e(x)|^2 ]
// where the weight W is optional
// In the vast majority of cases, optimisation will be expressed as
// a "least squares" problem, with f(x) = 0.5 * |W * e(x)|^2

// If constrained, constraints are expressed as error vectors where:
// h(x) = 0
// g(x) <= 0

// Methods can be designed for working directly with f(x), or e(x)
// Methods using e(x) will be similar to those with f(x), but will differ
// in how the hessian of f(x) = 0.5 * e(x)^T * e(x) is simplified or adjusted.

// 1. Gradient descent:
// dx = -alpha * grad_f
// Change in f(x) is guaranteed to be non-positive, ie: negative unless |grad_f| = 0, in which case f(x) is unchanged.
// This method is robust, but converges slowly
// It also only requires the jacobian of f(x), not the hessian

// 2. Newton-based methods:
// a) Using f(x)
//    => Solve: H * dx = -grad_f
// b) Using e(x) with hessian of e(x) available
//    grad_f = J^T * e
//    H = J^TJ + e^T * d2e/dx2
//    => Solve: (J^TJ + e^T * d2e/dx2) * dx = - J^T * e
// b) Using e(x) without hessian of e(x)
//    H approx= J^TJ
//    => Solve: (J^T J) * dx = -J^T * e
//    => dx = -J^+ e
// c) Using e(x) with a damped hessian
//    H = J^TJ + lambda I
//    => Solve: (J^T J + lambda I) * dx = -J^T * e
//    (optional: replace lambda I -> lambda * diag(J^T J))
//    With lambda dynamically adjusted to reflect accuracy of quadratic approximation of surface
//    i) lambda=0 : Same as (b)
//    ii) lambda->inf : Gradient descent
// d) Using f(x) with a damped hessian
//    => Solve: (H + lambda I) * dx = -grad_f
//    (optional: replace lambda I -> lambda * diag(H))
//
// For the above methods, the following names will be used:
// a) Newtons-method: Essentially root-finding for grad_f
// b) Gauss-newton
// c) Levenberg-Marquardt

// 4. Lagrangian multipliers
// Optimisation under equality constraints only
// Minise over z = [x, lambda] the lagrangian:
// L(x) = f(x) - lambda^T * h(x)
// Unconstrained optimisation of L(x) may use any available method
// that works with f(x)

// 5. Interior point method
// Optimisation over equality and inequality constraints
// Minimise over z = [x, lambda], the lagrangian:
// L(x) = f(x) - lambda^T * h(x) - alpha * sum_i log(g_i(x))
// While satisfying:
// mu_i = alpha / g_i(x)
// Unconstrained optimisation of L(x) may use any available method
// that works with f(x) and has an iterative update step
// (since the update step of this algorithm must simultaneously execute
//  the unconstrained optimisation step, and enforce constraints over
//  the step)


namespace owl {

template <typename X, typename Solver>
concept is_solver = requires(const Solver& solver, typename Solver::State& state, const X& initial_x) {
    { solver.start(initial_x) } -> std::convertible_to<typename Solver::State>;
    { solver.step(state) } -> std::convertible_to<bool>;
    { state.x } -> std::convertible_to<X>;
    { state.converged } -> std::convertible_to<bool>;
};

template <typename X, typename Solver>
requires is_solver<X, Solver>
bool solve(const X& initial_x, const Solver& solver, X& result) {
    auto state = solver.start(initial_x);
    while (solver.step(state)) {}
    result = state.x;
    return state.converged;
}

template <typename X, typename Solver>
requires is_solver<X, Solver>
typename Solver::State solve(const X& initial_x, const Solver& solver) {
    auto state = solver.start(initial_x);
    while (solver.step(state)) {}
    return state;
}

template <typename X, typename Solver>
requires is_solver<X, Solver>
typename Solver::State solve(const X& initial_x, const Solver& solver, std::vector<X>& trajectory) {
    auto state = solver.start(initial_x);
    trajectory.push_back(state.x);
    while (solver.step(state)) {
        trajectory.push_back(state.x);
    }
    return state;
}

} // namespace owl
