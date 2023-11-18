#pragma once

#include <optional>
#include <variant>
#include <Eigen/Eigenvalues>
#include "mathbox/diff/dfunction.h"
#include "mathbox/diff/manifold.h"


namespace mbox::unconstrained {


template <typename Scalar, int Dim>
bool verify_minima(const Matrix<Scalar, Dim, Dim>& hessian, Vector<Scalar, Dim>& perturbation, double perturbation_amount) {
    Eigen::SelfAdjointEigenSolver<Matrix<Scalar, Dim, Dim>> solver(hessian);
    bool minima = true;
    if constexpr (Dim == Eigen::Dynamic) {
        perturbation.resize(hessian.rows(), 1);
    }
    perturbation.setZero();
    for (std::size_t i = 0; i < solver.eigenvalues().size(); i++) {
        if (solver.eigenvalues()(i) < 0) {
            minima = false;
            if constexpr (Dim != Eigen::Dynamic) {
                perturbation += solver.eigenvectors().template block<Dim, 1>(0, i) * perturbation_amount;
            }
            if constexpr (Dim == Eigen::Dynamic) {
                perturbation += solver.eigenvectors().block(0, i, Dim, 1) * perturbation_amount;
            }
        }
    }
    return minima;
}


template <typename X>
class GradientDescent {
public:
    struct Config {
        std::size_t max_iterations;
        double convergence_delta_norm;
        double learning_rate;
        Config():
            max_iterations(100),
            convergence_delta_norm(1e-8),
            learning_rate(learning_rate)
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

    GradientDescent(const ScalarDFunction1<X>& f):
        f(f)
    {}

    State start(const X& initial_x) const {
        return State(initial_x);
    }

    bool step(State& state) const {
        auto delta_x = calculate_delta_x(state.x);
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
    manifold_delta<X> calculate_delta_x(const X& x) const {
        return -config_.learning_rate * f.gradient(x);
    }

    const ScalarDFunction1<X>& f;
    Config config_;
};

template <typename X>
class GaussNewton {
public:
    struct Config {
        std::size_t max_iterations;
        double convergence_delta_norm;
        bool verify_minima;
        double maxima_perturbation_amount;

        Config():
            max_iterations(100),
            convergence_delta_norm(1e-8),
            verify_minima(true),
            maxima_perturbation_amount(1e-3)
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

    GaussNewton(const ScalarDFunction2<X>& f):
        f(f)
    {}

    State start(const X& initial_x) const {
        return State(initial_x);
    }

    bool step(State& state) const {
        auto delta_x = calculate_delta_x(state.x);
        manifold_add(state.x, delta_x);
        state.converged = (delta_x.norm() <= config_.convergence_delta_norm);

        if (state.converged && config_.verify_minima) {
            auto hessian = f.hessian(state.x);
            manifold_delta<X> perturbation;
            if (!verify_minima(hessian, perturbation, config_.maxima_perturbation_amount)) {
                state.converged = false;
                manifold_add<X>(state.x, perturbation);
            }
        }
        if (!state.converged) {
            state.iteration++;
        }
        return !state.converged && state.iteration < config_.max_iterations;
    }

    Config& config() { return config_; }
    const Config& config() const { return config_; }

private:
    manifold_delta<X> calculate_delta_x(const X& x) const {
        auto grad_fx = f.gradient(x);
        auto hessian_fx = f.hessian(x);
        Eigen::JacobiSVD<manifold_hessian<X>> svd_solver(hessian_fx, Eigen::ComputeFullU | Eigen::ComputeFullV);
        manifold_delta<X> delta_x = -svd_solver.solve(grad_fx);
        return delta_x;
    }

    const ScalarDFunction2<X>& f;
    Config config_;
};

} // namespace mbox::unconstrained
