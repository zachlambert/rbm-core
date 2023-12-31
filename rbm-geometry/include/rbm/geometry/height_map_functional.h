#pragma once

#include "rbm/geometry/height_map.h"


namespace rbm {

class HeightMapGaussianMixture: public HeightMapInterface {
public:
    struct Component {
        Vector2d centre;
        double spread;
        double height;
        Component():
            centre(Vector2d::Zero()),
            spread(0),
            height(0)
        {}
        Component(const Vector2d& centre, double spread, double height):
            centre(centre),
            spread(spread),
            height(height)
        {}
        double z(const Vector2d& x) const {
            return height * std::exp(-0.5 * (x - centre).squaredNorm() / std::pow(spread, 2));
        }
        Vector2d gradient(const Vector2d& x) const {
            return -(x - centre) * z(x) / std::pow(spread, 2);
        }
    };

    const std::vector<Component>& components() const { return components_; }
    std::vector<Component>& components() { return components_; }

    bool contains(const Vector2d& position) const override {
        return true;
    }

    double height(const Vector2d& position) const override;
    Vector3d normal(const Vector2d& position) const override;

    std::unique_ptr<HeightMapInterface> clone() const override {
        return std::make_unique<HeightMapGaussianMixture>(*this);
    }

private:
    std::vector<Component> components_;
};

} // namespace rbm
