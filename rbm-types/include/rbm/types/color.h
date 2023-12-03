#pragma once

#include <array>


namespace rbm {

template <typename Scalar>
struct ColorRGB {
    union {
        std::array<Scalar, 3> data;
        struct {
            Scalar r;
            Scalar g;
            Scalar b;
        };
    };
    ColorRGB(): r(0), g(0), b(0) {}
    ColorRGB(Scalar r, Scalar g, Scalar b): r(r), g(g), b(b) {}

    template <typename Scalar_>
    ColorRGB<Scalar_> cast() const {
        ColorRGB<Scalar_> result;
        for (std::size_t i = 0; i < data.size(); i++) {
            result.data[i] = static_cast<Scalar_>(data[i]);
        }
        return result;
    }

    static ColorRGB White() {
        return ColorRGB(1, 1, 1);
    }
    static ColorRGB Red() {
        return ColorRGB(1, 0, 0);
    }
    static ColorRGB Green() {
        return ColorRGB(0, 1, 0);
    }
    static ColorRGB Blue() {
        return ColorRGB(0, 0, 1);
    }
};
typedef ColorRGB<float> ColorRGBf;
typedef ColorRGB<double> ColorRGBd;

} // namespace rbm
