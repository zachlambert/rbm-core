#include <owl/transform/transform.h>
#include <owl/transform/euler.h>
#include <owl/transform/io.h>
#include <owl/transform/adjoint.h>
#include <iostream>
#include <chrono>


double benchmark(std::size_t N, const std::function<void()>& func)
{
    std::chrono::high_resolution_clock clock;
    double s = 0;
    for (std::size_t i = 0; i < N; i++) {
        auto start = clock.now();
        func();
        auto end = clock.now();
        s += static_cast<double>(duration_cast<std::chrono::nanoseconds>(end - start).count());
    }
    s /= N;
    return s;
};

int main()
{
    {
        owl::Transform3d X1, X2, X3;
        X1.translation() << 1, 2, 3;
        X1.rotation() = owl::EulerRotation3d(0, 0.5, 1).toRotationMatrix();
        X2.translation() << 3, -2, 0.5;
        X2.rotation() = owl::EulerRotation3d(-0.1, 0.3, 0.2).toRotationMatrix();
        std::cout << "Tranform product: " << benchmark(1000, [&]() {
            X3 = X1 * X2;
        }) << "ns" << std::endl;
    }

    {
        owl::CompactTransform3d X1, X2, X3;
        X1.translation() << 1, 2, 3;
        X1.rotation() = owl::EulerRotation3d(0, 0.5, 1).toCompact();
        X2.translation() << 3, -2, 0.5;
        X2.rotation() = owl::EulerRotation3d(-0.1, 0.3, 0.2).toCompact();
        std::cout << "Tranform product (compact): " << benchmark(1000, [&]() {
            X3 = X1 * X2;
        }) << "ns" << std::endl;
    }

    {
        owl::Transform3d X;
        owl::Vector3d x, y;
        X.translation() << 1, 2, 3;
        X.rotation() = owl::EulerRotation3d(0, 0.5, 1).toRotationMatrix();
        x << 0, 1, 2;
        std::cout << "Tranform vector product: " << benchmark(1000, [&]() {
            y = X * x;
        }) << "ns" << std::endl;
    }

    {
        owl::CompactTransform3d X;
        owl::Vector3d x, y;
        X.translation() << 1, 2, 3;
        X.rotation() = owl::EulerRotation3d(0, 0.5, 1).toRotationMatrix();
        x << 0, 1, 2;
        std::cout << "Tranform vector product (compact): " << benchmark(1000, [&]() {
            y = X * x;
        }) << "ns" << std::endl;
    }

    {
        owl::Vector6d x_coords;
        x_coords << 1, 2, 3, 4, 5, 6;
        owl::Transform3d X;
        owl::CompactTransform3d X_compact;

        owl::LogTransform3d x(x_coords);
        std::cout << "Transform exp: " << benchmark(1000, [&]() {
            X = x.exp();
        }) << "ns" << std::endl;
        std::cout << "Transform exp (compact): " << benchmark(1000, [&]() {
            X_compact = x.exp_compact();
        }) << "ns" << std::endl;
    }

    {
        owl::Transform3d X;
        X.translation() << 1, 2, 3;
        X.rotation() = owl::EulerRotation3d(0, 0.5, 1).toRotationMatrix();
        owl::CompactTransform3d X_min;
        X_min.translation() << 1, 2, 3;
        X_min.rotation() = owl::EulerRotation3d(0, 0.5, 1).toCompact();

        owl::Vector6d x_coords;

        std::cout << "Transform log: " << benchmark(1000, [&]() {
            x_coords = owl::LogTransform3d(X).coords();
        }) << "ns" << std::endl;
        std::cout << "Transform log (compact): " << benchmark(1000, [&]() {
            x_coords = owl::LogTransform3d(X).coords();
        }) << "ns" << std::endl;
    }

    {
        owl::Transform3d X;
        X.translation() << 1, 2, 3;
        X.rotation() = owl::EulerRotation3d(0, 0.5, 1).toRotationMatrix();
        owl::CompactTransform3d X_min;
        X_min.translation() << 1, 2, 3;
        X_min.rotation() = owl::EulerRotation3d(0, 0.5, 1).toRotationMatrix();

        Eigen::Matrix<double, 6, 6> adjX_matrix;

        std::cout << "Adjoint matrix: " << benchmark(1000, [&]() {
            owl::AdjointTransform3d adjX(X);
            adjX_matrix = adjX.matrix();
        }) << "ns" << std::endl;
        std::cout << "Adjoint matrix (compact): " << benchmark(1000, [&]() {
            owl::AdjointTransform3d adjX(X_min);
            adjX_matrix = adjX.matrix();
        }) << "ns" << std::endl;
    }
}

