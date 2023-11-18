#pragma once

#include "mathbox/diff/dfunction.h"
#include <iostream>


namespace mbox {

template <typename X>
class InteriorPointSolver {
    static constexpr int XDim = manifold_dim<X>;

public:
    struct Config {
        std::size_t max_iterations;
        double initial_mu;
        double convergence_delta_norm;
        bool verify_constraints;
        Config():
            max_iterations(100),
            initial_mu(1e-4),
            convergence_delta_norm(1e-14),
            verify_constraints(true)
        {}
    };

    struct State {
        std::size_t iteration;
        X x;
        VectorXd lambda;
        double mu;
        bool converged;
        State(const X& initial_x):
            iteration(0),
            x(initial_x),
            mu(0),
            converged(false)
        {}
    };

    InteriorPointSolver(
        const ScalarDFunction2<X>& f,
        const DFunction2<X, VectorXd>& g
    ):
        f(f),
        g(g)
    {}

    std::function<double(const X&)> make_barrier_function() {
        double mu = config_.initial_mu;
        return [this, mu](const X& x) -> double {
            // Only using this for visualisation, so apply reasonable thresholding
            auto e = g(x);
            bool invalid = false;
            for (std::size_t i = 0; i < e.size(); i++) {
                if (e(i) <= 0) {
                    invalid = true;
                    break;
                }
            }
            if (invalid) {
                return 0.0;
            }
            else {
                return f(x) - mu * e.array().log().sum();
            }
        };
    };

    State solve() const {
        State state = start();
        while (step(f, state)) {}
        return state;
    }

    State start(const X& initial_x) const {
        State state(initial_x);
        state.mu = config_.initial_mu;
        auto g_initial = g(state.x);
        state.lambda = g_initial.array().inverse() * state.mu;
        return state;
    }

    bool step(State& state) const {
        auto g_x = g(state.x);
        if ((g_x.array() < 0).any()) {
            std::cout << "Failed, negative g(x): " << g_x.transpose() << std::endl;
            return false;
        }

        auto f_grad_x = f.gradient(state.x);
        auto f_hessian_x = f.hessian(state.x);
        auto dg_x = g.dfdx(state.x);
        auto d2g_x = g.d2fdx2(state.x);
        auto d2g_x_rev = d2g_x.reverse();

        const std::size_t x_dim = f_grad_x.size();
        const std::size_t g_dim = g_x.size();
        const std::size_t z_dim = x_dim + g_dim;

        state.mu = (state.lambda.array() * g_x.array()).sum() / g_x.size();
        state.lambda = g_x.array().inverse() * state.mu;
#if 0
        for (std::size_t i = 0; i < state.lambda.size(); i++) {
            if (std::fabs(state.lambda(i) * g_x(i) - state.mu) > 1e-6) {
                std::cout << "Failed, inconsistent mu: " << (state.lambda.array() * g_x.array()).transpose() << std::endl;
                return false;
            }
        }
#endif
        if ((state.lambda.array() < 0).any()) {
            std::cout << "Negative lambda: " << state.lambda.transpose() << std::endl;
            return false;
        }

        // Solve Az = b
        MatrixXd A(z_dim, z_dim);
        VectorXd b(z_dim, 1);
        A.setZero();
        b.setZero();

        A.topLeftCorner(x_dim, x_dim) = f_hessian_x - d2g_x_rev * state.lambda;
        A.topRightCorner(x_dim, g_dim) = -dg_x.transpose();
        A.bottomLeftCorner(g_dim, x_dim) = state.lambda.asDiagonal() * dg_x;
        A.bottomRightCorner(g_dim, g_dim) = g_x.asDiagonal();

        b.head(x_dim) = -f_grad_x + dg_x.transpose() * state.lambda;
        b.tail(g_dim).setConstant(state.mu);
        b.tail(g_dim) = b.tail(g_dim) - VectorXd(g_x.array() * state.lambda.array());

        Eigen::JacobiSVD<MatrixXd> svd_solver(A, Eigen::ComputeFullU | Eigen::ComputeFullV);
        VectorXd delta_z(z_dim);
        delta_z = svd_solver.solve(b);

        double alpha = 1;
        for (std::size_t i = 0; i < state.lambda.size(); i++) {
            double delta_lambda = delta_z(state.x.size() + i);
            if (delta_lambda >= 0) continue;
            double next_lambda = state.lambda(i) + delta_lambda;
            if (next_lambda > 0) continue;
            alpha = std::min(alpha, state.lambda(i) / std::fabs(delta_lambda));
        }
        // NOTE: At the moment, this doesn't work if this is disabled, since it will overshoot
        if (config_.verify_constraints) {
            while (true) {
                X test_x = manifold_sum<X>(state.x, alpha * delta_z.template head<XDim>());
                auto test_g_x = g(test_x);
                if ((test_g_x.array() > 0).all()) break;
                alpha *= 0.5;
            }
        }
        delta_z *= alpha;

        manifold_add(state.x, delta_z.head(x_dim));
        state.lambda += delta_z.tail(g_dim);

        state.converged = (delta_z.norm() <= config_.convergence_delta_norm);
        if (!state.converged) {
            state.iteration++;
        }

        return !state.converged && state.iteration < config_.max_iterations;
    }

    Config& config() { return config_; }
    const Config& config() const { return config_; }

private:
    const ScalarDFunction2<X>& f;
    const DFunction2<X, VectorXd>& g;
    Config config_;
};

} // namespace mbox
