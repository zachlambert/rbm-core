#include <gtest/gtest.h>
#include <rbm/transform/euler.h>
#include <rbm/transform/datapack.h>

TEST(Transform, Datapack)
{
    std::cout << "=== Datapack ===" << std::endl;

    rbm::EulerTransform3d pose;
    pose.translation() = rbm::Vector3d(1, 2, 3);
    pose.rotation() = rbm::EulerRotation3d(0, 0.5, -0.5);

    datapack::Json json;
    json.get() = pose;
    std::cout << json << std::endl;
}

