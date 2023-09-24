
#include "math/opt/unconstrained.h"
#include "math/opt/least_squares.h"
#include "math/opt/interior_point.h"
#include "math/opt/solve.h"
#include <iostream>

#include "gui/gui/window.h"


void plot_scalar_function(
    const std::string& label,
    const std::function<double(const math::Vector2d&)>& f,
    math::Vector2d x_lower,
    math::Vector2d x_upper,
    double resolution)
{
    math::Vector2d x_range = x_upper - x_lower;
    std::size_t rows = std::ceil(x_range.y() / resolution);
    std::size_t cols = std::ceil(x_range.x() / resolution);

    std::vector<double> fs(rows * cols);
    std::size_t k = 0;
    for (std::size_t i = 0; i < rows; i++) {
        for (std::size_t j = 0; j < cols; j++) {
            math::Vector2d x = x_lower;
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

template <typename... Lists>
std::tuple<math::Vector2d, math::Vector2d> get_min_max(double padding, const Lists&... xs_list) {
    math::Vector2d min = math::Vector2d::Constant(std::numeric_limits<double>::max());
    math::Vector2d max = math::Vector2d::Constant(-std::numeric_limits<double>::max());
    auto apply_xs = [&](const std::vector<math::Vector2d>& xs) {
        for (const auto& x: xs) {
            min.x() = std::min(min.x(), x.x());
            min.y() = std::min(min.y(), x.y());
            max.x() = std::max(max.x(), x.x());
            max.y() = std::max(max.y(), x.y());
        }
    };
    (apply_xs(xs_list),...);

    min -= math::Vector2d::Constant(padding);
    max += math::Vector2d::Constant(padding);
    return std::make_tuple(min, max);
}

#if 0
template <int XDim>
math::DifferentiableFunctiond<double, 2, math::Vectord<XDim>> make_f_denominator_layout(const math::Vectord<XDim>& target) {
    return math::DifferentiableFunctiond<double, 2, math::Vectord<XDim>>(
        [=](const math::Vectord<XDim>& x) -> double {
            math::Vectord<XDim> e = x - target;
            return 0.5 * e.transpose() * e;
        },
        [=](const math::Vectord<XDim>& x) -> math::Vectord<XDim> {
            return x - target;
        },
        [=](const math::Vectord<XDim>& x) -> Matrixd<XDim, XDim> {
            if constexpr (XDim == Eigen::Dynamic) {
                MatrixXd d2f(x.size(), x.size());
                d2f.setIdentity();
                return d2f;
            }
            if constexpr (XDim != Eigen::Dynamic) {
                return Matrixd<XDim, XDim>::Identity();
            }
        }
    );
}

template <int XDim>
math::DifferentiableFunctionn<math::VectorXd, 1, math::Vectord<XDim>> make_e_numerator_layout(const math::Vectord<XDim>& target) {
    return math::DifferentiableFunctionn<math::VectorXd, 1, math::Vectord<XDim>>(
        [=](const math::Vectord<XDim>& x) -> math::VectorXd {
            math::VectorXd result(x.size());
            result = x - target;
            return result;
        },
        [=](const math::Vectord<XDim>& x) -> Matrixd<Eigen::Dynamic, XDim> {
            Matrixd<Eigen::Dynamic, XDim> result(x.size(), x.size());
            result.setIdentity();
            return result;
        }
    );
}

math::DifferentiableFunctiond<math::VectorXd, 2, math::Vector2d> make_g_denominator_layout(const std::vector<std::vector<double>>& coefficients) {
    return math::DifferentiableFunctiond<math::VectorXd, 2, math::Vector2d>(
        [=](const math::Vector2d& x) -> math::VectorXd {
            math::VectorXd result(coefficients.size(), 1);
            for (std::size_t n = 0; n < coefficients.size(); n++) {
                double y = 0;
                for (std::size_t i = 0; i < coefficients[n].size(); i++) {
                    y += coefficients[n][i] * std::pow(x.x(), i);
                }
                result(n, 0) = x.y() - y;
            }
            return result;
        },
        [=](const math::Vector2d& x) -> Matrixd<2, Eigen::Dynamic> {
            Matrixd<2, Eigen::Dynamic> result(2, coefficients.size());
            for (std::size_t n = 0; n < coefficients.size(); n++) {
                double dydx = 0;
                for (std::size_t i = 1; i < coefficients[n].size(); i++) {
                    dydx += coefficients[n][i] * i * std::pow(x.x(), i - 1);
                }
                result.template block<2, 1>(0, n) = math::Vector2d(-dydx, 1);
            }
            return result;
        },
        [=](const math::Vector2d& x) -> std::vector<Matrix2d> {
            std::vector<Matrix2d> result(coefficients.size());
            for (std::size_t n = 0; n < coefficients.size(); n++) {
                double d2ydx2 = 0;
                for (std::size_t i = 2; i < coefficients[n].size(); i++) {
                    d2ydx2 += coefficients[n][i] * (i - 1) * i * std::pow(x.x(), i - 2);
                }
                result[n].setZero();
                result[n](0, 0) = d2ydx2;
            }
            return result;
        }
    );
}
#endif

template <int XDim>
std::function<double(const math::Vectord<XDim>&)> error_half_square_norm_f(const math::Vectord<XDim>& target) {
    return [target](const math::Vectord<XDim>& x) {
        return 0.5 * (target - x).squaredNorm();
    };
}

template <int XDim>
std::function<math::Vectord<XDim>(const math::Vectord<XDim>&)> error_half_square_norm_gradient(const math::Vectord<XDim>& target) {
    return [target](const math::Vectord<XDim>& x) {
        return x - target;
    };
}

template <int XDim>
std::function<math::Matrixd<XDim, XDim>(const math::Vectord<XDim>&)> error_half_square_norm_hessian(const math::Vectord<XDim>& target) {
    return [target](const math::Vectord<XDim>& x) {
        return math::Matrixd<XDim, XDim>::Identity(x.size(), x.size());
    };
}

template <int XDim>
std::function<math::Vectord<XDim>(const math::Vectord<XDim>&)> error_vector_f(const math::Vectord<XDim>& target) {
    return [target](const math::Vectord<XDim>& x) {
        return target - x;
    };
}

template <int XDim>
std::function<math::Matrixd<XDim, XDim>(const math::Vectord<XDim>&)> error_vector_dfdx(const math::Vectord<XDim>& target) {
    return [target](const math::Vectord<XDim>& x) {
        return -math::Matrixd<XDim, XDim>::Identity(x.rows(), x.rows());
    };
}

std::function<math::VectorXd(const math::Vector2d&)> polynomial_error_f(const std::vector<std::vector<double>> coefficients) {
    return [=](const math::Vector2d& x) {
        math::VectorXd result(coefficients.size());
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

std::function<math::MatrixX2d(const math::Vector2d&)> polynomial_error_dfdx(const std::vector<std::vector<double>> coefficients) {
    return [=](const math::Vector2d& x) {
        math::MatrixX2d result(coefficients.size(), 2);
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

std::function<math::MTensord<Eigen::Dynamic, 2, 2>(const math::Vector2d&)> polynomial_error_d2fdx2(const std::vector<std::vector<double>> coefficients) {
    return [=](const math::Vector2d& x) {
        math::MTensord<Eigen::Dynamic, 2, 2> result(coefficients.size(), 2, 2);
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
    {
        math::Vector2d initial_x = math::Vector2d::Zero();
        math::Vector2d target = math::Vector2d(1, 2);

        auto f1 = math::ScalarDFunction1<math::Vector2d>(
            error_half_square_norm_f(target),
            error_half_square_norm_gradient(target)
        );
        auto f2 = math::ScalarDFunction2<math::Vector2d>(
            error_half_square_norm_f(target),
            error_half_square_norm_gradient(target),
            error_half_square_norm_hessian(target)
        );

        math::unconstrained::GradientDescent<math::Vector2d> solver1(f1);
        solver1.config().learning_rate = 1e-2;
        solver1.config().max_iterations = 1000;

        math::unconstrained::GaussNewton<math::Vector2d> solver2(f2);

        std::vector<math::Vector2d> xs1, xs2;
        math::solve(initial_x, solver1, xs1);
        math::solve(initial_x, solver2, xs2);
        math::Vector2d x_min, x_max;
        std::tie(x_min, x_max) = get_min_max(1, xs1, xs2);

        create_plot([&](){
            plot_scalar_function("f", f1.f(), x_min, x_max, 0.01);
            ImPlot::SetNextMarkerStyle(ImPlotMarker_Circle);
            ImPlot::PlotLine("xs_grad_descent", &xs1[0].x(), &xs1[0].y(), xs1.size(), ImPlotFlags_None, 0, sizeof(math::Vector2d));
            ImPlot::SetNextMarkerStyle(ImPlotMarker_Circle);
            ImPlot::PlotLine("xs_gauss_newton", &xs2[0].x(), &xs2[0].y(), xs2.size(), ImPlotFlags_None, 0, sizeof(math::Vector2d));
        });
    }

    {
        std::vector<std::vector<double>> coefficients;
        math::Vector2d initial_x = math::Vector2d::Zero();
        coefficients.push_back({ 2, -2, 0.6 });
        auto g = math::DFunction1(polynomial_error_f(coefficients), polynomial_error_dfdx(coefficients));
        auto f = [&](const math::Vector2d& x) {
            return 0.5 * g(x).squaredNorm();
        };

        math::least_squares::GradientDescent<math::Vector2d> solver1(g);
        math::least_squares::NewtonsMethod<math::Vector2d> solver2(g);
        math::least_squares::LevenbergMarquardt<math::Vector2d> solver3(g);

        std::vector<math::Vector2d> xs1, xs2, xs3;
        math::solve(initial_x, solver1, xs1);
        math::solve(initial_x, solver2, xs2);
        math::solve(initial_x, solver3, xs3);
        math::Vector2d x_min, x_max;
        std::tie(x_min, x_max) = get_min_max(1, xs1, xs2, xs3);

        auto g_solution = polynomial_error_solution(coefficients);
        std::vector<std::vector<math::Vector2d>> g_xs(coefficients.size());
        for (std::size_t i = 0; i < coefficients.size(); i++) {
            const auto& gi_sol = g_solution[i];
            auto& gi_xs = g_xs[i];
            for (double x = x_min.x(); x <= x_max.x(); x += 0.05) {
                double y = gi_sol(x);
                if (y > x_max.y() || y < x_min.y()) continue;
                gi_xs.push_back(math::Vector2d(x, y));
            }
        }

        create_plot([&](){
            plot_scalar_function("f", f, x_min, x_max, 0.01);
            ImPlot::SetNextMarkerStyle(ImPlotMarker_Circle);
            ImPlot::PlotLine("xs_steepest", &xs1[0].x(), &xs1[0].y(), xs1.size(), ImPlotFlags_None, 0, sizeof(math::Vector2d));
            ImPlot::SetNextMarkerStyle(ImPlotMarker_Circle);
            ImPlot::PlotLine("xs_newton", &xs2[0].x(), &xs2[0].y(), xs2.size(), ImPlotFlags_None, 0, sizeof(math::Vector2d));
            ImPlot::SetNextMarkerStyle(ImPlotMarker_Circle);
            ImPlot::PlotLine("xs_lma", &xs3[0].x(), &xs3[0].y(), xs3.size(), ImPlotFlags_None, 0, sizeof(math::Vector2d));
            for (std::size_t i = 0; i < g_xs.size(); i++) {
                char label[32];
                snprintf(label, sizeof(label), "g%zu", i);
                ImPlot::PlotLine(label, &g_xs[i][0].x(), &g_xs[i][0].y(), g_xs[i].size(), ImPlotFlags_None, 0, sizeof(math::Vector2d));
            }
        });
    }

    {
        math::Vector2d initial_x(2, 1);

        math::Vector2d target(1.5, 0);
        auto f = math::ScalarDFunction2<math::Vector2d>(
            error_half_square_norm_f(target),
            error_half_square_norm_gradient(target),
            error_half_square_norm_hessian(target)
        );

        std::vector<std::vector<double>> coefficients;
        coefficients.push_back({ 2, -2, 0.6 });
        coefficients.push_back({ 7, -4, 0 });
        auto g = math::DFunction2(
            polynomial_error_f(coefficients),
            polynomial_error_dfdx(coefficients),
            polynomial_error_d2fdx2(coefficients)
        );

        math::InteriorPointSolver<math::Vector2d> solver(f, g);
        solver.config().initial_mu = 1e-4;
        auto b = solver.make_barrier_function();

        std::vector<math::Vector2d> xs;
        math::solve(initial_x, solver, xs);
        math::Vector2d x_min, x_max;
        std::tie(x_min, x_max) = get_min_max(1, xs);

        auto g_solution = polynomial_error_solution(coefficients);
        std::vector<std::vector<math::Vector2d>> g_xs(coefficients.size());
        for (std::size_t i = 0; i < coefficients.size(); i++) {
            const auto& gi_sol = g_solution[i];
            auto& gi_xs = g_xs[i];
            for (double x = x_min.x(); x <= x_max.x(); x += 0.05) {
                double y = gi_sol(x);
                if (y > x_max.y() || y < x_min.y()) continue;
                gi_xs.push_back(math::Vector2d(x, y));
            }
        }

        create_plot([&](){
            plot_scalar_function("f", f, x_min, x_max, 0.01);
            plot_scalar_function("b", b, x_min, x_max, 0.01);
            for (std::size_t i = 0; i < coefficients.size(); i++) {
                char label[32];
                snprintf(label, sizeof(label), "g%zu", i);
                ImPlot::PlotLine(label, &g_xs[i][0].x(), &g_xs[i][0].y(), g_xs[i].size(), ImPlotFlags_None, 0, sizeof(math::Vector2d));
            }
            ImPlot::SetNextMarkerStyle(ImPlotMarker_Circle);
            ImPlot::PlotLine("xs", &xs[0].x(), &xs[0].y(), xs.size(), ImPlotFlags_None, 0, sizeof(math::Vector2d));
        });
    }
    {
        math::VectorXd initial_x(10);
        initial_x.setZero();

        math::VectorXd target(10);
        target.setRandom();
        std::cout << "Target: " << target << std::endl;
        auto f1 = math::ScalarDFunction1<math::VectorXd>(
            error_half_square_norm_f(target),
            error_half_square_norm_gradient(target)
        );
        auto f2 = math::ScalarDFunction2<math::VectorXd>(
            error_half_square_norm_f(target),
            error_half_square_norm_gradient(target),
            error_half_square_norm_hessian(target)
        );

        math::unconstrained::GradientDescent<math::VectorXd> solver1(f1);
        solver1.config().learning_rate = 1e-2;
        solver1.config().max_iterations = 10000;

        auto result1 = math::solve<math::VectorXd>(initial_x, solver1);
        std::cout << "Grad descent result:\n" << result1 << std::endl;

        math::unconstrained::GaussNewton<math::VectorXd> solver2(f2);
        auto result2 = math::solve<math::VectorXd>(initial_x, solver2);
        std::cout << "Gauss newton result:\n" << result2 << std::endl;
    }
}
