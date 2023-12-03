#include <gtest/gtest.h>
#include <rbm/transform/transform.h>
#include <rbm/transform/euler.h>
#include <rbm/transform/io.h>

TEST(Transform, Transform2)
{
    rbm::Vector3d dx_coords;
    dx_coords << M_PI/2, 1, 0;
    rbm::LogTransform2d dx(dx_coords);

    {
        rbm::Transform2d X = rbm::Transform2d::Identity();
        X = X * dx.exp();
        std::cout << rbm::EulerTransform2d(X) << std::endl;
    }
}

