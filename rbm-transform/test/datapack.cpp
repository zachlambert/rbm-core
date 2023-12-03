#include <gtest/gtest.h>
#include <mbox/transform/euler.h>
#include <mbox/transform/datapack.h>

TEST(Transform, Datapack)
{
    std::cout << "=== Datapack ===" << std::endl;

    mbox::EulerTransform3d pose;
    pose.translation() = mbox::Vector3d(1, 2, 3);
    pose.rotation() = mbox::EulerRotation3d(0, 0.5, -0.5);

    datapack::Json json;
    json.get() = pose;
    std::cout << json << std::endl;
}

