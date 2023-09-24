
#include "math/transform/angle.h"
#include "math/transform/transform.h"
#include "math/transform/adjoint.h"
#include "math/transform/json.h"
#include "math/transform/io.h"
#include <chrono>


void test_angle() {
    std::cout << "=== Angle ===" << std::endl;

    math::Angled a = -M_PI;
    math::Angled b = 0.5 * M_PI;
    auto c = a - b;
    std::cout << a << " - " << b << " = " << c << "\n";
}

void test_serialize() {
    std::cout << "=== Serialize ===" << std::endl;

    math::EulerTransform3d pose;
    pose.translation() = math::Vector3d(1, 2, 3);
    pose.rotation() = math::EulerRotation3d(0, 0.5, -0.5);

    core::Json json;
    json.get() = pose;
    std::cout << json << std::endl;
}

void test_transform2() {
    std::cout << "=== TRANSFORM 2 ===" << std::endl;

    math::Vector3d dx_coords;
    dx_coords << M_PI/2, 1, 0;
    math::LogTransform2d dx(dx_coords);

    {
        math::Transform2d X = math::Transform2d::Identity();
        X = X * dx.exp();
        std::cout << math::EulerTransform2d(X) << std::endl;
    }
}

void test_transform3() {
    std::cout << "=== TRANSFORM 3 ===" << std::endl;

    math::Vector6d dx_coords;
    dx_coords << 0, 0, 0.5, 1, 0, 0;
    math::LogTransform3d dx(dx_coords);

    {
        math::Transform3d X = math::Transform3d::Identity();
        X = X * dx.exp();
        std::cout << math::EulerTransform3d(X) << std::endl;
    }

    {
        math::CompactTransform3d X = math::CompactTransform3d::Identity();
        X = X * dx.exp_compact();
        std::cout << math::EulerTransform3d(X) << std::endl;
    }

    auto benchmark = [&](size_t N, const std::function<void()>& func) {
        std::chrono::high_resolution_clock clock;
        double s = 0;
        for (size_t i = 0; i < N; i++) {
            auto start = clock.now();
            func();
            auto end = clock.now();
            s += static_cast<double>(duration_cast<std::chrono::nanoseconds>(end - start).count());
        }
        s /= N;
        return s;
    };

    {
        math::Transform3d X1, X2, X3;
        X1.translation() << 1, 2, 3;
        X1.rotation() = math::EulerRotation3d(0, 0.5, 1).toRotationMatrix();
        X2.translation() << 3, -2, 0.5;
        X2.rotation() = math::EulerRotation3d(-0.1, 0.3, 0.2).toRotationMatrix();
        std::cout << "Tranform product: " << benchmark(1000, [&]() {
            X3 = X1 * X2;
        }) << "ns" << std::endl;
    }

    {
        math::CompactTransform3d X1, X2, X3;
        X1.translation() << 1, 2, 3;
        X1.rotation() = math::EulerRotation3d(0, 0.5, 1).toCompact();
        X2.translation() << 3, -2, 0.5;
        X2.rotation() = math::EulerRotation3d(-0.1, 0.3, 0.2).toCompact();
        std::cout << "Tranform product (compact): " << benchmark(1000, [&]() {
            X3 = X1 * X2;
        }) << "ns" << std::endl;
    }

    {
        math::Transform3d X;
        math::Vector3d x, y;
        X.translation() << 1, 2, 3;
        X.rotation() = math::EulerRotation3d(0, 0.5, 1).toRotationMatrix();
        x << 0, 1, 2;
        std::cout << "Tranform vector product: " << benchmark(1000, [&]() {
            y = X * x;
        }) << "ns" << std::endl;
    }

    {
        math::CompactTransform3d X;
        math::Vector3d x, y;
        X.translation() << 1, 2, 3;
        X.rotation() = math::EulerRotation3d(0, 0.5, 1).toRotationMatrix();
        x << 0, 1, 2;
        std::cout << "Tranform vector product (compact): " << benchmark(1000, [&]() {
            y = X * x;
        }) << "ns" << std::endl;
    }

    {
        math::Vector6d x_coords;
        x_coords << 1, 2, 3, 4, 5, 6;
        math::Transform3d X;
        math::CompactTransform3d X_compact;

        math::LogTransform3d x(x_coords);
        std::cout << "Transform exp: " << benchmark(1000, [&]() {
            X = x.exp();
        }) << "ns" << std::endl;
        std::cout << "Transform exp (compact): " << benchmark(1000, [&]() {
            X_compact = x.exp_compact();
        }) << "ns" << std::endl;
    }

    {
        math::Transform3d X;
        X.translation() << 1, 2, 3;
        X.rotation() = math::EulerRotation3d(0, 0.5, 1).toRotationMatrix();
        math::CompactTransform3d X_min;
        X_min.translation() << 1, 2, 3;
        X_min.rotation() = math::EulerRotation3d(0, 0.5, 1).toCompact();

        math::Vector6d x_coords;

        std::cout << "Transform log: " << benchmark(1000, [&]() {
            x_coords = math::LogTransform3d(X).coords();
        }) << "ns" << std::endl;
        std::cout << "Transform log (compact): " << benchmark(1000, [&]() {
            x_coords = math::LogTransform3d(X).coords();
        }) << "ns" << std::endl;
    }

    {
        math::Transform3d X;
        X.translation() << 1, 2, 3;
        X.rotation() = math::EulerRotation3d(0, 0.5, 1).toRotationMatrix();
        math::CompactTransform3d X_min;
        X_min.translation() << 1, 2, 3;
        X_min.rotation() = math::EulerRotation3d(0, 0.5, 1).toRotationMatrix();

        Eigen::Matrix<double, 6, 6> adjX_matrix;

        std::cout << "Adjoint matrix: " << benchmark(1000, [&]() {
            math::AdjointTransform3d adjX(X);
            adjX_matrix = adjX.matrix();
        }) << "ns" << std::endl;
        std::cout << "Adjoint matrix (compact): " << benchmark(1000, [&]() {
            math::AdjointTransform3d adjX(X_min);
            adjX_matrix = adjX.matrix();
        }) << "ns" << std::endl;
    }
}

int main() {
    test_angle();
    test_serialize();
    test_transform3();
    test_transform2();
}
