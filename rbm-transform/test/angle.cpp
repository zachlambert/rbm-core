#include <gtest/gtest.h>
#include <rbm/transform/angle.h>

TEST(Transform, Angle)
{
    std::cout << "=== Angle ===" << std::endl;

    rbm::Angled a = -M_PI;
    rbm::Angled b = 0.5 * M_PI;
    auto c = a - b;
    std::cout << a << " - " << b << " = " << c << "\n";
}

