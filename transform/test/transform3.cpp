#include <gtest/gtest.h>
#include <owl/transform/transform.h>
#include <owl/transform/euler.h>
#include <owl/transform/io.h>


TEST(Transform, Transform3)
{
    std::cout << "=== TRANSFORM 3 ===" << std::endl;

    owl::Vector6d dx_coords;
    dx_coords << 0, 0, 0.5, 1, 0, 0;
    owl::LogTransform3d dx(dx_coords);

    {
        owl::Transform3d X = owl::Transform3d::Identity();
        X = X * dx.exp();
        std::cout << owl::EulerTransform3d(X) << std::endl;
    }

    {
        owl::CompactTransform3d X = owl::CompactTransform3d::Identity();
        X = X * dx.exp_compact();
        std::cout << owl::EulerTransform3d(X) << std::endl;
    }
}

