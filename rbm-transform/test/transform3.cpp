#include <gtest/gtest.h>
#include <rbm/transform/transform.h>
#include <rbm/transform/euler.h>
#include <rbm/transform/io.h>


TEST(Transform, Transform3)
{
    std::cout << "=== TRANSFORM 3 ===" << std::endl;

    rbm::Vector6d dx_coords;
    dx_coords << 0, 0, 0.5, 1, 0, 0;
    rbm::LogTransform3d dx(dx_coords);

    {
        rbm::Transform3d X = rbm::Transform3d::Identity();
        X = X * dx.exp();
        std::cout << rbm::EulerTransform3d(X) << std::endl;
    }

    {
        rbm::CompactTransform3d X = rbm::CompactTransform3d::Identity();
        X = X * dx.exp_compact();
        std::cout << rbm::EulerTransform3d(X) << std::endl;
    }
}

