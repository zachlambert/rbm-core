#include <gtest/gtest.h>
#include <rbm/cpp/optional_tuple.h>
#include <string>

TEST (OptionalTuple, Main) {
    rbm::OptionalTuple<int, double, std::string> x;

    x.set<0>(20);
    EXPECT_EQ(x.value<0>(), 20);

    EXPECT_THROW(x.value<1>(), std::bad_optional_access);
    EXPECT_THROW(x.set<2>("Hello"), std::bad_optional_access);

    x.set<1>(1.234);
    EXPECT_EQ(x.value<1>(), 1.234);

    x.set<2>("Hello!");
    EXPECT_EQ(x.value<2>(), "Hello!");

    x.resize<1>();
    ASSERT_THROW(x.value<1>(), std::bad_optional_access);

    x.resize<0>();
    ASSERT_EQ(x.size(), 0);
    ASSERT_EQ(x.capacity(), 3);

    auto in_range = std::make_tuple<int, double>(12, -10.2);
    x.set_range<2>(in_range);
    ASSERT_EQ(x.value<0>(), 12);
    ASSERT_EQ(x.value<1>(), -10.2);

    rbm::OptionalTuple<int, double> y;
    y.set_range<2>(x);
    ASSERT_EQ(y.value<0>(), x.value<0>());
    ASSERT_EQ(y.value<1>(), x.value<1>());

    y.value<0>() = 100;
    y.value<1>() = 100;
    x.set_range<2>(y);
    ASSERT_EQ(x.value<0>(), 100);
    ASSERT_EQ(x.value<1>(), 100);

    auto z = x.subset<int, double>();
    ASSERT_EQ(z.value<0>(), 100);
    ASSERT_EQ(z.value<1>(), 100);
}
