#include <gtest/gtest.h>
#include <owl/transform/transform.h>
#include <owl/transform/euler.h>
#include <owl/transform/io.h>

TEST(Transform, Transform2)
{
    owl::Vector3d dx_coords;
    dx_coords << M_PI/2, 1, 0;
    owl::LogTransform2d dx(dx_coords);

    {
        owl::Transform2d X = owl::Transform2d::Identity();
        X = X * dx.exp();
        std::cout << owl::EulerTransform2d(X) << std::endl;
    }
}

