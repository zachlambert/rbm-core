
#include "rbm/opt/unconstrained.h"
#include "rbm/opt/least_squares.h"
#include "rbm/opt/interior_point.h"
#include "rbm/opt/solve.h"
#include <iostream>

#include "rbm/gui/window.h"

// TODO
#if 0
#include "gui/gui/window.h"


void plot_scalar_function(
    const std::string& label,
    const std::function<double(const rbm::Vector2d&)>& f,
    rbm::Vector2d x_lower,
    rbm::Vector2d x_upper,
    double resolution)
{
    rbm::Vector2d x_range = x_upper - x_lower;
    std::size_t rows = std::ceil(x_range.y() / resolution);
    std::size_t cols = std::ceil(x_range.x() / resolution);

    std::vector<double> fs(rows * cols);
    std::size_t k = 0;
    for (std::size_t i = 0; i < rows; i++) {
        for (std::size_t j = 0; j < cols; j++) {
            rbm::Vector2d x = x_lower;
            x.x() = x_lower.x() + ((double)j + 0.5) * x_range.x() / cols;
            x.y() = x_upper.y() - ((double)i + 0.5) * x_range.y() / rows;
            fs[k] = f(x);
            k++;
        }
    }

    ImPlot::PushColormap(ImPlotColormap_Viridis);
    ImPlot::PlotHeatmap(label.c_str(), fs.data(), rows, cols, 0, 0, "", ImPlotPoint(x_lower.x(), x_lower.y()), ImPlotPoint(x_upper.x(), x_upper.y()));
    ImPlot::PopColormap();
}

void create_plot(const std::function<void()>& plot) {
    gui::Window window;
    if (!window.init()) return;
    while (window.running()) {
        window.poll_events();
        window.render_start();
        ImGui::Begin("Result");
        if (ImPlot::BeginPlot("##plot", ImVec2(-1, 600))) {
            plot();
            ImPlot::EndPlot();
        }
        ImGui::End();
        window.render_end();
    }
}
#endif

template <typename... Lists>
std::tuple<rbm::Vector2d, rbm::Vector2d> get_min_max(double padding, const Lists&... xs_list) {
    rbm::Vector2d min = rbm::Vector2d::Constant(std::numeric_limits<double>::max());
    rbm::Vector2d max = rbm::Vector2d::Constant(-std::numeric_limits<double>::max());
    auto apply_xs = [&](const std::vector<rbm::Vector2d>& xs) {
        for (const auto& x: xs) {
            min.x() = std::min(min.x(), x.x());
            min.y() = std::min(min.y(), x.y());
            max.x() = std::max(max.x(), x.x());
            max.y() = std::max(max.y(), x.y());
        }
    };
    (apply_xs(xs_list),...);

    min -= rbm::Vector2d::Constant(padding);
    max += rbm::Vector2d::Constant(padding);
    return std::make_tuple(min, max);
}

template <int XDim>
std::function<double(const rbm::Vectord<XDim>&)> error_half_square_norm_f(const rbm::Vectord<XDim>& target) {
    return [target](const rbm::Vectord<XDim>& x) {
        return 0.5 * (target - x).squaredNorm();
    };
}

template <int XDim>
std::function<rbm::Vectord<XDim>(const rbm::Vectord<XDim>&)> error_half_square_norm_gradient(const rbm::Vectord<XDim>& target) {
    return [target](const rbm::Vectord<XDim>& x) {
        return x - target;
    };
}

template <int XDim>
std::function<rbm::Matrixd<XDim, XDim>(const rbm::Vectord<XDim>&)> error_half_square_norm_hessian(const rbm::Vectord<XDim>& target) {
    return [target](const rbm::Vectord<XDim>& x) {
        return rbm::Matrixd<XDim, XDim>::Identity(x.size(), x.size());
    };
}

template <int XDim>
std::function<rbm::Vectord<XDim>(const rbm::Vectord<XDim>&)> error_vector_f(const rbm::Vectord<XDim>& target) {
    return [target](const rbm::Vectord<XDim>& x) {
        return target - x;
    };
}

template <int XDim>
std::function<rbm::Matrixd<XDim, XDim>(const rbm::Vectord<XDim>&)> error_vector_dfdx(const rbm::Vectord<XDim>& target) {
    return [target](const rbm::Vectord<XDim>& x) {
        return -rbm::Matrixd<XDim, XDim>::Identity(x.rows(), x.rows());
    };
}

std::function<rbm::VectorXd(const rbm::Vector2d&)> polynomial_error_f(const std::vector<std::vector<double>> coefficients) {
    return [=](const rbm::Vector2d& x) {
        rbm::VectorXd result(coefficients.size());
        for (std::size_t n = 0; n < coefficients.size(); n ++) {
            double y = 0;
            for (std::size_t i = 0; i < coefficients[n].size(); i++) {
                y += coefficients[n][i] * std::pow(x.x(), i);
            }
            result(n, 0) = x.y() - y;
        }
        return result;
    };
}

std::function<rbm::MatrixX2d(const rbm::Vector2d&)> polynomial_error_dfdx(const std::vector<std::vector<double>> coefficients) {
    return [=](const rbm::Vector2d& x) {
        rbm::MatrixX2d result(coefficients.size(), 2);
        result.block(0, 1, coefficients.size(), 1).setOnes();
        for (std::size_t n = 0; n < coefficients.size(); n ++) {
            double dydx = 0;
            for (std::size_t i = 1; i < coefficients[n].size(); i++) {
                dydx += coefficients[n][i] * i * std::pow(x.x(), i - 1);
            }
            result(n, 0) = -dydx;
        }
        return result;
    };
}

std::function<rbm::MTensord<Eigen::Dynamic, 2, 2>(const rbm::Vector2d&)> polynomial_error_d2fdx2(const std::vector<std::vector<double>> coefficients) {
    return [=](const rbm::Vector2d& x) {
        rbm::MTensord<Eigen::Dynamic, 2, 2> result(coefficients.size(), 2, 2);
        result[1].setZero();
        result[0].block(0, 1, coefficients.size(), 1).setZero();
        for (std::size_t n = 0; n < coefficients.size(); n ++) {
            double d2ydx2 = 0;
            for (std::size_t i = 2; i < coefficients[n].size(); i++) {
                d2ydx2 += coefficients[n][i] * i * (i - 1) * std::pow(x.x(), i - 2);
            }
            result[0](n, 0) = -d2ydx2;
        }
        return result;
    };
}

std::vector<std::function<double(double)>> polynomial_error_solution(const std::vector<std::vector<double>>& coefficients) {
    std::vector<std::function<double(double)>> result;
    for (const auto& coefficients_n: coefficients) {
        // Copy coefficients_n
        result.push_back([=](double x) -> double {
            double y = 0;
            for (std::size_t i = 0; i < coefficients_n.size(); i++) {
                y += std::pow(x, i) * coefficients_n[i];
            }
            return y;
        });
    }
    return result;
};


int main() {
    sviz::Window window("Mathbox optimisation example");
    {
        rbm::Vector2d initial_x = rbm::Vector2d::Zero();
        rbm::Vector2d target = rbm::Vector2d(1, 2);

        auto f1 = rbm::ScalarDFunction1<rbm::Vector2d>(
            error_half_square_norm_f(target),
            error_half_square_norm_gradient(target)
        );
        auto f2 = rbm::ScalarDFunction2<rbm::Vector2d>(
            error_half_square_norm_f(target),
            error_half_square_norm_gradient(target),
            error_half_square_norm_hessian(target)
        );

        rbm::unconstrained::GradientDescent<rbm::Vector2d> solver1(f1);
        solver1.config().learning_rate = 1e-2;
        solver1.config().max_iterations = 1000;

        rbm::unconstrained::GaussNewton<rbm::Vector2d> solver2(f2);

        std::vector<rbm::Vector2d> xs1, xs2;
        rbm::solve(initial_x, solver1, xs1);
        rbm::solve(initial_x, solver2, xs2);
        rbm::Vector2d x_min, x_max;
        std::tie(x_min, x_max) = get_min_max(1, xs1, xs2);

#if 0
        create_plot([&](){
            plot_scalar_function("f", f1.f(), x_min, x_max, 0.01);
            ImPlot::SetNextMarkerStyle(ImPlotMarker_Circle);
            ImPlot::PlotLine("xs_grad_descent", &xs1[0].x(), &xs1[0].y(), xs1.size(), ImPlotFlags_None, 0, sizeof(rbm::Vector2d));
            ImPlot::SetNextMarkerStyle(ImPlotMarker_Circle);
            ImPlot::PlotLine("xs_gauss_newton", &xs2[0].x(), &xs2[0].y(), xs2.size(), ImPlotFlags_None, 0, sizeof(rbm::Vector2d));
        });
#endif
    }

    {
        std::vector<std::vector<double>> coefficients;
        rbm::Vector2d initial_x = rbm::Vector2d::Zero();
        coefficients.push_back({ 2, -2, 0.6 });
        auto g = rbm::DFunction1(polynomial_error_f(coefficients), polynomial_error_dfdx(coefficients));
        auto f = [&](const rbm::Vector2d& x) {
            return 0.5 * g(x).squaredNorm();
        };

        rbm::least_squares::GradientDescent<rbm::Vector2d> solver1(g);
        rbm::least_squares::NewtonsMethod<rbm::Vector2d> solver2(g);
        rbm::least_squares::LevenbergMarquardt<rbm::Vector2d> solver3(g);

        std::vector<rbm::Vector2d> xs1, xs2, xs3;
        rbm::solve(initial_x, solver1, xs1);
        rbm::solve(initial_x, solver2, xs2);
        rbm::solve(initial_x, solver3, xs3);
        rbm::Vector2d x_min, x_max;
        std::tie(x_min, x_max) = get_min_max(1, xs1, xs2, xs3);

        auto g_solution = polynomial_error_solution(coefficients);
        std::vector<std::vector<rbm::Vector2d>> g_xs(coefficients.size());
        for (std::size_t i = 0; i < coefficients.size(); i++) {
            const auto& gi_sol = g_solution[i];
            auto& gi_xs = g_xs[i];
            for (double x = x_min.x(); x <= x_max.x(); x += 0.05) {
                double y = gi_sol(x);
                if (y > x_max.y() || y < x_min.y()) continue;
                gi_xs.push_back(rbm::Vector2d(x, y));
            }
        }

#if 0
        create_plot([&](){
            plot_scalar_function("f", f, x_min, x_max, 0.01);
            ImPlot::SetNextMarkerStyle(ImPlotMarker_Circle);
            ImPlot::PlotLine("xs_steepest", &xs1[0].x(), &xs1[0].y(), xs1.size(), ImPlotFlags_None, 0, sizeof(rbm::Vector2d));
            ImPlot::SetNextMarkerStyle(ImPlotMarker_Circle);
            ImPlot::PlotLine("xs_newton", &xs2[0].x(), &xs2[0].y(), xs2.size(), ImPlotFlags_None, 0, sizeof(rbm::Vector2d));
            ImPlot::SetNextMarkerStyle(ImPlotMarker_Circle);
            ImPlot::PlotLine("xs_lma", &xs3[0].x(), &xs3[0].y(), xs3.size(), ImPlotFlags_None, 0, sizeof(rbm::Vector2d));
            for (std::size_t i = 0; i < g_xs.size(); i++) {
                char label[32];
                snprintf(label, sizeof(label), "g%zu", i);
                ImPlot::PlotLine(label, &g_xs[i][0].x(), &g_xs[i][0].y(), g_xs[i].size(), ImPlotFlags_None, 0, sizeof(rbm::Vector2d));
            }
        });
#endif
    }

    {
        rbm::Vector2d initial_x(2, 1);

        rbm::Vector2d target(1.5, 0);
        auto f = rbm::ScalarDFunction2<rbm::Vector2d>(
            error_half_square_norm_f(target),
            error_half_square_norm_gradient(target),
            error_half_square_norm_hessian(target)
        );

        std::vector<std::vector<double>> coefficients;
        coefficients.push_back({ 2, -2, 0.6 });
        coefficients.push_back({ 7, -4, 0 });
        auto g = rbm::DFunction2(
            polynomial_error_f(coefficients),
            polynomial_error_dfdx(coefficients),
            polynomial_error_d2fdx2(coefficients)
        );

        rbm::InteriorPointSolver<rbm::Vector2d> solver(f, g);
        solver.config().initial_mu = 1e-4;
        auto b = solver.make_barrier_function();

        std::vector<rbm::Vector2d> xs;
        rbm::solve(initial_x, solver, xs);
        rbm::Vector2d x_min, x_max;
        std::tie(x_min, x_max) = get_min_max(1, xs);

        auto g_solution = polynomial_error_solution(coefficients);
        std::vector<std::vector<rbm::Vector2d>> g_xs(coefficients.size());
        for (std::size_t i = 0; i < coefficients.size(); i++) {
            const auto& gi_sol = g_solution[i];
            auto& gi_xs = g_xs[i];
            for (double x = x_min.x(); x <= x_max.x(); x += 0.05) {
                double y = gi_sol(x);
                if (y > x_max.y() || y < x_min.y()) continue;
                gi_xs.push_back(rbm::Vector2d(x, y));
            }
        }

#if 0
        create_plot([&](){
            plot_scalar_function("f", f, x_min, x_max, 0.01);
            plot_scalar_function("b", b, x_min, x_max, 0.01);
            for (std::size_t i = 0; i < coefficients.size(); i++) {
                char label[32];
                snprintf(label, sizeof(label), "g%zu", i);
                ImPlot::PlotLine(label, &g_xs[i][0].x(), &g_xs[i][0].y(), g_xs[i].size(), ImPlotFlags_None, 0, sizeof(rbm::Vector2d));
            }
            ImPlot::SetNextMarkerStyle(ImPlotMarker_Circle);
            ImPlot::PlotLine("xs", &xs[0].x(), &xs[0].y(), xs.size(), ImPlotFlags_None, 0, sizeof(rbm::Vector2d));
        });
#endif
    }
    {
        rbm::VectorXd initial_x(10);
        initial_x.setZero();

        rbm::VectorXd target(10);
        target.setRandom();
        std::cout << "Target: " << target << std::endl;
        auto f1 = rbm::ScalarDFunction1<rbm::VectorXd>(
            error_half_square_norm_f(target),
            error_half_square_norm_gradient(target)
        );
        auto f2 = rbm::ScalarDFunction2<rbm::VectorXd>(
            error_half_square_norm_f(target),
            error_half_square_norm_gradient(target),
            error_half_square_norm_hessian(target)
        );

        rbm::unconstrained::GradientDescent<rbm::VectorXd> solver1(f1);
        solver1.config().learning_rate = 1e-2;
        solver1.config().max_iterations = 10000;

        auto result1 = rbm::solve<rbm::VectorXd>(initial_x, solver1);
        std::cout << "Grad descent result:\n" << result1 << std::endl;

        rbm::unconstrained::GaussNewton<rbm::VectorXd> solver2(f2);
        auto result2 = rbm::solve<rbm::VectorXd>(initial_x, solver2);
        std::cout << "Gauss newton result:\n" << result2 << std::endl;
    }
}
