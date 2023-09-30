#pragma once

#include "owl/geometry/height_map.h"


namespace owl {

class HeightMapGaussianMixture: public HeightMapInterface {
public:
    struct Component {
        math::Vector2d centre;
        double spread;
        double height;
        Component():
            centre(math::Vector2d::Zero()),
            spread(0),
            height(0)
        {}
        Component(const math::Vector2d& centre, double spread, double height):
            centre(centre),
            spread(spread),
            height(height)
        {}
        double z(const math::Vector2d& x) const {
            return height * std::exp(-0.5 * (x - centre).squaredNorm() / std::pow(spread, 2));
        }
        math::Vector2d gradient(const math::Vector2d& x) const {
            return -(x - centre) * z(x) / std::pow(spread, 2);
        }
    };

    const std::vector<Component>& components() const { return components_; }
    std::vector<Component>& components() { return components_; }

    bool contains(const math::Vector2d& position) const override {
        return true;
    }

    double height(const math::Vector2d& position) const override;
    math::Vector3d normal(const math::Vector2d& position) const override;

    std::unique_ptr<HeightMapInterface> clone() const override {
        return std::make_unique<HeightMapGaussianMixture>(*this);
    }

private:
    std::vector<Component> components_;
};

} // namespace owl
