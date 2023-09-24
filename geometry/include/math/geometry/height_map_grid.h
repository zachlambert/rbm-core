#pragma once

#include "math/geometry/height_map.h"
#include "math/geometry/grid.h"


namespace math {

class HeightMapGrid: public HeightMapInterface {
public:

    math::Vector2d& origin() { return dimensions_.origin; }
    const math::Vector2d& origin() const { return dimensions_.origin; }
    math::Vector2d& size() { return dimensions_.size; }
    const math::Vector2d& size() const { return dimensions_.size; }

    void resize(const math::Vector2i& divisions) {
        height_data_.resize(divisions, GridStorageType::CORNER);
        normal_data_.resize(divisions, GridStorageType::CELL);
    }

    double& get_height(const math::Vector2i& index) {
        return height_data_.get(index);
    }
    double get_height(const math::Vector2i& index) const {
        return height_data_.get(index);
    }
    math::Vector3d& get_normal(const math::Vector2i& index) {
        return normal_data_.get(index);
    }
    math::Vector3d get_normal(const math::Vector2i& index) const {
        return normal_data_.get(index);
    }

    bool contains(const math::Vector2d& position) const override {
        return true;
    }

    double height(const math::Vector2d& position) const override;
    math::Vector3d normal(const math::Vector2d& position) const override;

    bool has_normal() const override {
        return true;
    }

    math::Vector2d cell_size() const {
        return dimensions_.size.array() / height_data_.divisions().cast<double>().array();
    }

    virtual std::unique_ptr<HeightMapInterface> clone() const override {
        return std::make_unique<HeightMapGrid>(*this);
    }

private:
    GridData2d<double> height_data_;
    GridData2d<math::Vector3d> normal_data_;
    GridDimensions2d dimensions_;
};

} // namespace math
