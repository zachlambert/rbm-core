#include <gtest/gtest.h>
#include <mathbox/transform/transform.h>
#include <mathbox/transform/euler.h>
#include <mathbox/transform/io.h>

TEST(Transform, Transform2)
{
    mbox::Vector3d dx_coords;
    dx_coords << M_PI/2, 1, 0;
    mbox::LogTransform2d dx(dx_coords);

    {
        mbox::Transform2d X = owl::Transform2d::Identity();
        X = X * dx.exp();
        std::cout << mbox::EulerTransform2d(X) << std::endl;
    }
}

