#pragma once

#include <memory>
#include <iostream>
#include <rbm/diff/dfunction.h>


// Least squares solver work with the function for the error directly
// They will use numerator-layout convention for derivatives (ie: jacobian layout)
// since hessian aren't needed.

namespace rbm::least_squares {

// Gradient descent with f(x) = 0.5 * e(x)^T * e(x)
// grad_fx = [e^T * de/dx]^T = (de/dx)^T * e
// Therefore:
// delta_x = -alpha * (de/dx)^T * e
// Instead of an "optimal" alpha, use an empirically chosen constant for now
// Can perhaps use a policy for decreasing this over time

template <typename X>
class GradientDescent {
    static constexpr int XDim = manifold_dim<X>;
public:
    struct Config {
        std::size_t max_iterations;
        double convergence_delta_norm;
        double initial_alpha;

        Config():
            max_iterations(100),
            convergence_delta_norm(1e-8),
            initial_alpha(0.1)
        {}
    };

    struct State {
        std::size_t iteration;
        X x;
        double alpha;
        bool converged;
        State(const X& initial_x):
            iteration(0),
            x(initial_x),
            alpha(0),
            converged(false)
        {}

        friend std::ostream& operator<<(std::ostream& os, const State& state) {
            os << "Iteration: " << state.iteration << "\n";
            os << "X: " << state.x << "\n";
            os << "alpha: " << state.alpha << "\n";
            os << "converged: " << (state.converged ? "true" : "false");
            return os;
        }
    };

    GradientDescent(const DFunction1<X, VectorXd>& e):
        e(e)
    {}

    State start(const X& initial_x) const {
        State state(initial_x);
        state.alpha = config_.initial_alpha;
        return state;
    }

    bool step(State& state) const {
        Vectord<XDim> delta_x = calculate_delta_x(state.x, state.alpha);
        manifold_add(state.x, delta_x);
        state.converged = (delta_x.norm() <= config_.convergence_delta_norm);
        if (!state.converged) {
            state.iteration++;
        }
        return !state.converged && state.iteration < config_.max_iterations;
    }

    Config& config() { return config_; }
    const Config& config() const { return config_; }

private:
    Vectord<XDim> calculate_delta_x(const X& x, double& alpha) const {
        auto e_x = e(x);
        auto de_x = e.dfdx(x);
        // alpha unchanged
        return -alpha * de_x.transpose() * e_x;
    }

    const DFunction1<X, VectorXd>& e;
    Config config_;
};

// Gauss-newton for f(x) = 0.5 * e(x)^T * e(x)
// using the approximate hessian (de/dx)^T(de/dx)
// Such that:
// delta_x = -[(de/dx)^T(de/dx)]^{-1}(de/dx)^T * e
// delta_x = -(de/de)^{+} * e
// Equivalent to multivariate newton's method for trying to solving e(x) = 0
// Also equivalent to exact least-squares for linear Ax = b, if we work with
// the first-order approximation of e(x): (de/dx) * delta_x = -e

template <typename X>
class NewtonsMethod {
    static constexpr int XDim = manifold_dim<X>;
public:
    struct Config {
        std::size_t max_iterations;
        double convergence_delta_norm;
        Config():
            max_iterations(100),
            convergence_delta_norm(1e-8)
        {}
    };

    struct State {
        std::size_t iteration;
        X x;
        bool converged;
        State(const X& initial_x):
            iteration(0),
            x(initial_x),
            converged(false)
        {}

        friend std::ostream& operator<<(std::ostream& os, const State& state) {
            os << "Iteration: " << state.iteration << "\n";
            os << "X: " << state.x << "\n";
            os << "converged: " << (state.converged ? "true" : "false");
            return os;
        }
    };

    NewtonsMethod(const DFunction1<X, VectorXd>& e):
        e(e)
    {}

    State start(const X& initial_x) const {
        return State(initial_x);
    }

    bool step(State& state) const {
        Vectord<XDim> delta_x = calculate_delta_x(state.x);
        manifold_add(state.x, delta_x);
        state.converged = (delta_x.norm() <= config_.convergence_delta_norm);
        if (!state.converged) {
            state.iteration++;
        }
        return !state.converged && state.iteration < config_.max_iterations;
    }

    Config& config() { return config_; }
    const Config& config() const { return config_; }

private:
    Vectord<XDim> calculate_delta_x(const X& x) const {
        // Least squares solution of: (de/dx) * delta_x = -e
        auto e_x = e(x);
        auto de_x = e.dfdx(x);
        Eigen::JacobiSVD<Matrix<double, -1, XDim>> svd_solver(de_x, Eigen::ComputeFullU | Eigen::ComputeFullV);
        auto delta_x = -svd_solver.solve(e_x);
        return delta_x;
    }

    const DFunction1<X, VectorXd>& e;
    Config config_;
};


template <typename X>
class LevenbergMarquardt {
    static constexpr int XDim = manifold_dim<X>;
public:
    struct Config {
        std::size_t max_iterations;
        double convergence_delta_norm;
        bool verify_minima;
        double maxima_perturbation_amount;
        double initial_lambda;
        double lambda_multiplier;

        Config():
            max_iterations(100),
            convergence_delta_norm(1e-8),
            verify_minima(true),
            maxima_perturbation_amount(1e-3),
            initial_lambda(0.1),
            lambda_multiplier(1.5)
        {}
    };

    struct State {
        std::size_t iteration;
        X x;
        double lambda;
        bool converged;
        State(const X& initial_x):
            iteration(0),
            x(initial_x),
            lambda(0),
            converged(false)
        {}
    };

    LevenbergMarquardt(const DFunction1<X, VectorXd>& e):
        e(e)
    {}

    State start(const X& initial_x) const {
        State state(initial_x);
        state.lambda = config_.initial_lambda;
        return state;
    }

    bool step(State& state) const {
        Vectord<XDim> delta_x = calculate_delta_x(state.x, state.lambda);
        manifold_add(state.x, delta_x);
        state.converged = (delta_x.norm() <= config_.convergence_delta_norm);

        if (!state.converged) {
            state.iteration++;
        }
        return !state.converged && state.iteration < config_.max_iterations;
    }

    Config& config() { return config_; }
    const Config& config() const { return config_; }

private:
    Vectord<XDim> calculate_delta_x(const X& x, double& lambda) const {
        VectorXd e_x = e(x);
        Matrixd<Eigen::Dynamic, XDim> J_x = e.dfdx(x);
        std::size_t x_dim = J_x.cols();

        Matrixd<XDim, XDim> damped_H_f_x = J_x.transpose() * J_x + lambda * Matrixd<XDim, XDim>::Identity(x_dim, x_dim);
        Vectord<XDim> grad_f_x = J_x.transpose() * e_x;
        Eigen::JacobiSVD<Matrixd<XDim, XDim>> svd_solver(damped_H_f_x, Eigen::ComputeFullU | Eigen::ComputeFullV);
        Vectord<XDim> delta_x = -svd_solver.solve(grad_f_x);

        X trial_x = manifold_sum(x, delta_x);
        VectorXd e_trial_x = e(trial_x);
        double f_x = e_x.transpose() * e_x;
        double f_trial_x = e_trial_x.transpose() * e_trial_x;
        double delta_f_actual = f_trial_x - f_x;
        double delta_f_pred = (e_x.transpose() * J_x * delta_x).value();
        double ratio = delta_f_actual / delta_f_pred;
        if (ratio < 0.25) {
            // Step size small, can reduce lambda
            lambda *= (1.0 / config_.lambda_multiplier);
        }
        else if (ratio > 0.5) {
            // Step size large, increase lambda
            lambda *= config_.lambda_multiplier;
        }

        return delta_x;
    }

    const DFunction1<X, VectorXd>& e;
    mutable X initial_x;
    Config config_;
};


} // namespace rbm::least_squares
