#include <gtest/gtest.h>
#include <mbox/transform/angle.h>

TEST(Transform, Angle)
{
    std::cout << "=== Angle ===" << std::endl;

    mbox::Angled a = -M_PI;
    mbox::Angled b = 0.5 * M_PI;
    auto c = a - b;
    std::cout << a << " - " << b << " = " << c << "\n";
}

