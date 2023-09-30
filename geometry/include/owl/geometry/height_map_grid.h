#pragma once

#include "owl/geometry/height_map.h"
#include "owl/geometry/grid.h"


namespace owl {

class HeightMapGrid: public HeightMapInterface {
public:

    Vector2d& origin() { return dimensions_.origin; }
    const Vector2d& origin() const { return dimensions_.origin; }
    Vector2d& size() { return dimensions_.size; }
    const Vector2d& size() const { return dimensions_.size; }

    void resize(const Vector2i& divisions) {
        height_data_.resize(divisions, GridStorageType::CORNER);
        normal_data_.resize(divisions, GridStorageType::CELL);
    }

    double& get_height(const Vector2i& index) {
        return height_data_.get(index);
    }
    double get_height(const Vector2i& index) const {
        return height_data_.get(index);
    }
    Vector3d& get_normal(const Vector2i& index) {
        return normal_data_.get(index);
    }
    Vector3d get_normal(const Vector2i& index) const {
        return normal_data_.get(index);
    }

    bool contains(const Vector2d& position) const override {
        return true;
    }

    double height(const Vector2d& position) const override;
    Vector3d normal(const Vector2d& position) const override;

    bool has_normal() const override {
        return true;
    }

    Vector2d cell_size() const {
        return dimensions_.size.array() / height_data_.divisions().cast<double>().array();
    }

    virtual std::unique_ptr<HeightMapInterface> clone() const override {
        return std::make_unique<HeightMapGrid>(*this);
    }

private:
    GridData2d<double> height_data_;
    GridData2d<Vector3d> normal_data_;
    GridDimensions2d dimensions_;
};

} // namespace owl
