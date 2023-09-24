#include <gtest/gtest.h>
#include <owl/transform/euler.h>
#include <owl/transform/serialize.h>

TEST(Transform, Serialize)
{
    std::cout << "=== Serialize ===" << std::endl;

    owl::EulerTransform3d pose;
    pose.translation() = owl::Vector3d(1, 2, 3);
    pose.rotation() = owl::EulerRotation3d(0, 0.5, -0.5);

    parrot::Json json;
    json.get() = pose;
    std::cout << json << std::endl;
}

