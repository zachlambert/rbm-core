#include <gtest/gtest.h>
#include <mathbox/transform/transform.h>
#include <mathbox/transform/euler.h>
#include <mathbox/transform/io.h>


TEST(Transform, Transform3)
{
    std::cout << "=== TRANSFORM 3 ===" << std::endl;

    mbox::Vector6d dx_coords;
    dx_coords << 0, 0, 0.5, 1, 0, 0;
    mbox::LogTransform3d dx(dx_coords);

    {
        mbox::Transform3d X = mbox::Transform3d::Identity();
        X = X * dx.exp();
        std::cout << mbox::EulerTransform3d(X) << std::endl;
    }

    {
        mbox::CompactTransform3d X = mbox::CompactTransform3d::Identity();
        X = X * dx.exp_compact();
        std::cout << mbox::EulerTransform3d(X) << std::endl;
    }
}

