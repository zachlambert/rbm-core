#include <gtest/gtest.h>
#include <mathbox/transform/euler.h>
#include <mathbox/transform/serialize.h>

TEST(Transform, Serialize)
{
    std::cout << "=== Serialize ===" << std::endl;

    mbox::EulerTransform3d pose;
    pose.translation() = mbox::Vector3d(1, 2, 3);
    pose.rotation() = mbox::EulerRotation3d(0, 0.5, -0.5);

    datapack::Json json;
    json.get() = pose;
    std::cout << json << std::endl;
}

