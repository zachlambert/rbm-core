#include <gtest/gtest.h>
#include <owl/transform/angle.h>

TEST(Transform, Angle)
{
    std::cout << "=== Angle ===" << std::endl;

    owl::Angled a = -M_PI;
    owl::Angled b = 0.5 * M_PI;
    auto c = a - b;
    std::cout << a << " - " << b << " = " << c << "\n";
}

