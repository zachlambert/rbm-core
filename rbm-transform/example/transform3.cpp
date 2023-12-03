#include <rbm/transform/transform.h>
#include <rbm/transform/euler.h>
#include <rbm/transform/io.h>
#include <rbm/transform/adjoint.h>
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
        rbm::Transform3d X1, X2, X3;
        X1.translation() << 1, 2, 3;
        X1.rotation() = rbm::EulerRotation3d(0, 0.5, 1).toRotationMatrix();
        X2.translation() << 3, -2, 0.5;
        X2.rotation() = rbm::EulerRotation3d(-0.1, 0.3, 0.2).toRotationMatrix();
        std::cout << "Tranform product: " << benchmark(1000, [&]() {
            X3 = X1 * X2;
        }) << "ns" << std::endl;
    }

    {
        rbm::CompactTransform3d X1, X2, X3;
        X1.translation() << 1, 2, 3;
        X1.rotation() = rbm::EulerRotation3d(0, 0.5, 1).toCompact();
        X2.translation() << 3, -2, 0.5;
        X2.rotation() = rbm::EulerRotation3d(-0.1, 0.3, 0.2).toCompact();
        std::cout << "Tranform product (compact): " << benchmark(1000, [&]() {
            X3 = X1 * X2;
        }) << "ns" << std::endl;
    }

    {
        rbm::Transform3d X;
        rbm::Vector3d x, y;
        X.translation() << 1, 2, 3;
        X.rotation() = rbm::EulerRotation3d(0, 0.5, 1).toRotationMatrix();
        x << 0, 1, 2;
        std::cout << "Tranform vector product: " << benchmark(1000, [&]() {
            y = X * x;
        }) << "ns" << std::endl;
    }

    {
        rbm::CompactTransform3d X;
        rbm::Vector3d x, y;
        X.translation() << 1, 2, 3;
        X.rotation() = rbm::EulerRotation3d(0, 0.5, 1).toRotationMatrix();
        x << 0, 1, 2;
        std::cout << "Tranform vector product (compact): " << benchmark(1000, [&]() {
            y = X * x;
        }) << "ns" << std::endl;
    }

    {
        rbm::Vector6d x_coords;
        x_coords << 1, 2, 3, 4, 5, 6;
        rbm::Transform3d X;
        rbm::CompactTransform3d X_compact;

        rbm::LogTransform3d x(x_coords);
        std::cout << "Transform exp: " << benchmark(1000, [&]() {
            X = x.exp();
        }) << "ns" << std::endl;
        std::cout << "Transform exp (compact): " << benchmark(1000, [&]() {
            X_compact = x.exp_compact();
        }) << "ns" << std::endl;
    }

    {
        rbm::Transform3d X;
        X.translation() << 1, 2, 3;
        X.rotation() = rbm::EulerRotation3d(0, 0.5, 1).toRotationMatrix();
        rbm::CompactTransform3d X_min;
        X_min.translation() << 1, 2, 3;
        X_min.rotation() = rbm::EulerRotation3d(0, 0.5, 1).toCompact();

        rbm::Vector6d x_coords;

        std::cout << "Transform log: " << benchmark(1000, [&]() {
            x_coords = rbm::LogTransform3d(X).coords();
        }) << "ns" << std::endl;
        std::cout << "Transform log (compact): " << benchmark(1000, [&]() {
            x_coords = rbm::LogTransform3d(X).coords();
        }) << "ns" << std::endl;
    }

    {
        rbm::Transform3d X;
        X.translation() << 1, 2, 3;
        X.rotation() = rbm::EulerRotation3d(0, 0.5, 1).toRotationMatrix();
        rbm::CompactTransform3d X_min;
        X_min.translation() << 1, 2, 3;
        X_min.rotation() = rbm::EulerRotation3d(0, 0.5, 1).toRotationMatrix();

        Eigen::Matrix<double, 6, 6> adjX_matrix;

        std::cout << "Adjoint matrix: " << benchmark(1000, [&]() {
            rbm::AdjointTransform3d adjX(X);
            adjX_matrix = adjX.matrix();
        }) << "ns" << std::endl;
        std::cout << "Adjoint matrix (compact): " << benchmark(1000, [&]() {
            rbm::AdjointTransform3d adjX(X_min);
            adjX_matrix = adjX.matrix();
        }) << "ns" << std::endl;
    }
}

