
#include "math/algebra/dfunction.h"
#include <iostream>

template <typename X>
void test_expansion(const math::ScalarDFunction1<X>& f, const X& x_1, const X& x_2) {
    std::cout << "f(" << x_1 << ") = " << f(x_1) << std::endl;
    std::cout << "f(" << x_2 << ") = " << f(x_2) << std::endl;

    auto delta_x = math::manifold_difference(x_1, x_2);
    double y_approx = f(x_1) + f.gradient(x_1).transpose() * delta_x;
    std::cout << "f(" << x_2 << ") ~= " << y_approx << " to first order approximation" << std::endl;
}

template <typename X>
void test_expansion(const math::ScalarDFunction2<X>& f, const X& x_1, const X& x_2) {
    std::cout << "f(" << x_1 << ") = " << f(x_1) << std::endl;
    std::cout << "f(" << x_2 << ") = " << f(x_2) << std::endl;

    auto delta_x = math::manifold_difference(x_1, x_2);
    double y_approx = f(x_1)
        + f.gradient(x_1).transpose() * delta_x
        + 0.5 * (delta_x.transpose() * f.hessian(x_1) * delta_x).value();
    std::cout << "f(" << x_2 << ") ~= " << y_approx << " to second order approximation" << std::endl;
}

template <typename X, typename Y>
void test_expansion(const math::DFunction1<X, Y>& f, const X& x_1, const X& x_2) {
    std::cout << "f(" << x_1 << ") = " << f(x_1) << std::endl;
    std::cout << "f(" << x_2 << ") = " << f(x_2) << std::endl;

    auto delta_x = math::manifold_difference(x_1, x_2);
    Y y_approx = math::manifold_sum(f(x_1), f.dfdx(x_1) * delta_x);
    std::cout << "f(" << x_2 << ") ~= " << y_approx << " to first order approximation" << std::endl;
}

template <typename X, typename Y>
void test_expansion(const math::DFunction2<X, Y>& f, const X& x_1, const X& x_2) {
    std::cout << "f(" << x_1 << ") = " << f(x_1) << std::endl;
    std::cout << "f(" << x_2 << ") = " << f(x_2) << std::endl;

    auto delta_x = math::manifold_difference(x_1, x_2);
    Y y_approx = math::manifold_sum(f(x_1), f.dfdx(x_1) * delta_x + f.d2fdx2(x_1) * delta_x * delta_x);
    std::cout << "f(" << x_2 << ") ~= " << y_approx << " to second order approximation" << std::endl;
}

int main() {
    std::cout << "f(x) = 0.5 * x^2" << std::endl;

    double x_1 = 1.2;
    double x_2 = 1.4;

    math::ScalarDFunction1<double> fa_1(
        [](const double& x) { return 0.5 * std::pow(x, 2); },
        [](const double& x) { return math::Vector1d(x); }
    );
    test_expansion(fa_1, x_1, x_2);

    math::ScalarDFunction2<double> fa_2(
        [](const double& x) { return 0.5 * std::pow(x, 2); },
        [](const double& x) { return math::Vector1d(x); },
        [](const double& x) { return math::Matrix1d::Identity(); }
    );
    test_expansion(fa_2, x_1, x_2);

    std::cout << std::endl << "f(x) = [cos(x) // sin(x)]" << std::endl;

    math::DFunction1<double, math::Vector2d> fb_1(
        [](const double& x) { return math::Vector2d(std::cos(x), std::sin(x)); },
        [](const double& x) { return math::Vector2d(-std::sin(x), std::cos(x)); }
    );
    test_expansion(fb_1, x_1, x_2);

    math::DFunction2<double, math::Vector2d> fb_2(
        [](const double& x) { return math::Vector2d(std::cos(x), std::sin(x)); },
        [](const double& x) { return math::Vector2d(-std::sin(x), std::cos(x)); },
        [](const double& x) -> math::MTensord<2, 1, 1> {
        math::MTensord<2, 1, 1> result;
            result[0] = math::Vector2d(-std::cos(x), -std::sin(x));
            return result;
        }
    );
    test_expansion(fb_2, x_1, x_2);

    std::cout << std::endl << "f(x) = 0.5 * (x - b)^T * A * (x - b)" << std::endl;

    math::Matrix2d A;
    A << 1, 0, 0, 1;
    math::Vector2d b;
    b << 1, 1;

    math::Vector2d x_3(1, 1);
    math::Vector2d x_4(1.2, 0.8);

    math::ScalarDFunction1<math::Vector2d> fc_1(
        [&](const math::Vector2d& x) -> double {
            return 0.5 * (x - b).transpose() * A * (x - b);
        },
        [&](const math::Vector2d& x) -> math::Vector2d {
            return A.transpose() * (x - b);
        }
    );
    test_expansion(fc_1, x_3, x_4);

    math::ScalarDFunction2<math::Vector2d> fc_2(
        [&](const math::Vector2d& x) -> double {
            return 0.5 * (x - b).transpose() * A * (x - b);
        },
        [&](const math::Vector2d& x) -> math::Vector2d {
            return A.transpose() * (x - b);
        },
        [&](const math::Vector2d& x) {
            return A;
        }
    );
    test_expansion(fc_2, x_3, x_4);

    std::cout << std::endl << "f(x) = [e^(x_1) * cos(x_2) // e^(x_1) * sin(x_2) ]" << std::endl;

    math::DFunction1<math::Vector2d, math::Vector2d> fd_1(
        [&](const math::Vector2d& x) -> math::Vector2d {
            math::Vector2d result;
            result(0) = std::exp(x(0)) * std::cos(x(1));
            result(1) = std::exp(x(0)) * std::sin(x(1));
            return result;
        },
        [&](const math::Vector2d& x) -> math::Matrix2d {
            math::Matrix2d result;
            result(0, 0) = std::exp(x(0)) * std::cos(x(1));
            result(1, 0) = std::exp(x(0)) * std::sin(x(1));
            result(0, 1) = -std::exp(x(0)) * std::sin(x(1));
            result(1, 1) = std::exp(x(0)) * std::cos(x(1));
            return result;
        }
    );
    test_expansion(fd_1, x_3, x_4);

    math::DFunction2<math::Vector2d, math::Vector2d> fd_2(
        [&](const math::Vector2d& x) -> math::Vector2d {
            math::Vector2d result;
            result(0) = std::exp(x(0)) * std::cos(x(1));
            result(1) = std::exp(x(0)) * std::sin(x(1));
            return result;
        },
        [&](const math::Vector2d& x) -> math::Matrix2d {
            math::Matrix2d result;
            result(0, 0) = std::exp(x(0)) * std::cos(x(1));
            result(1, 0) = std::exp(x(0)) * std::sin(x(1));
            result(0, 1) = -std::exp(x(0)) * std::sin(x(1));
            result(1, 1) = std::exp(x(0)) * std::cos(x(1));
            return result;
        },
        [&](const math::Vector2d& x) -> math::MTensord<2, 2, 2> {
            math::MTensord<2, 2, 2> result;
            result[0](0, 0) = std::exp(x(0)) * std::cos(x(1));
            result[0](1, 0) = std::exp(x(0)) * std::sin(x(1));
            result[0](0, 1) = -std::exp(x(0)) * std::sin(x(1));
            result[0](1, 1) = std::exp(x(0)) * std::cos(x(1));
            result[1](0, 0) = -std::exp(x(0)) * std::sin(x(1));
            result[1](1, 0) = std::exp(x(0)) * std::cos(x(1));
            result[1](0, 1) = -std::exp(x(0)) * std::cos(x(1));
            result[1](1, 1) = -std::exp(x(0)) * std::sin(x(1));
            return result;
        }
    );
    test_expansion(fd_2, x_3, x_4);
}
