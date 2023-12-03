#pragma once

#include <vector>
#include <iostream>
#include <rbm/types/matrix.h>

namespace rbm {

template <int N>
void step_grid_index(Vectori<N>& index, const Vectori<N>& end) {
    int i = 0;
    while (true) {
        index(i)++;
        if (index(i) != end(i)){
            break;
        }
        i++;
        if (i == N) {
            return;
        }
    }
    while (i != 0) {
        i--;
        index(i) = 0;
    }
}

enum class GridStorageType {
    CELL,
    CORNER
};

template <typename Scalar, int N>
struct GridDimensions {
    Vector<Scalar, N> origin;
    Vector<Scalar, N> size;

    GridDimensions():
        origin(Vector<Scalar, N>::Zero()),
        size(Vector<Scalar, N>::Zero())
    {}
    Vector<Scalar, N> to_normalized(const Vector<Scalar, N>& position) const {
        return (position - origin).array() / size.array();
    }
    Vector<Scalar, N> from_normalized(const Vector<Scalar, N>& position) const {
        return origin + size.array() * position.array();
    }
    bool contains(const Vector<Scalar, N>& position) const {
        return (position >= origin).any() && (position <= origin + size).any();
    }
};

using GridDimensions2f = GridDimensions<float, 2>;
using GridDimensions2d = GridDimensions<double, 2>;
using GridDimensions3f = GridDimensions<float, 3>;
using GridDimensions3d = GridDimensions<double, 3>;

template <typename Scalar, int N>
struct GridPoint {
    Vectori<N> a;
    Vectori<N> b;
    Vector<Scalar, N> fraction;
};

using GridPoint2f = GridPoint<float, 2>;
using GridPoint2d = GridPoint<double, 2>;
using GridPoint3f = GridPoint<float, 3>;
using GridPoint3d = GridPoint<double, 3>;

template <typename T, typename Scalar, int N>
class GridData {
public:
    void resize(const Vectori<N>& divisions, GridStorageType storage_type = GridStorageType::CELL) {
        divisions_ = divisions;
        storage_type_ = storage_type;

        Vectori<N> num_cells;
        if (storage_type == GridStorageType::CELL) {
            num_cells = divisions;
        } else {
            num_cells = divisions.array() + 1;
        }

        std::size_t data_size = 1;
        for (std::size_t i = 0; i < N; i++) {
            data_size *= num_cells(i);
        }
        data_.resize(data_size);
    }

    T& get(const Vectori<N>& index) {
        std::size_t data_i = 0;
        std::size_t stride = 1;
        for (std::size_t i = 0; i < N; i++) {
            data_i += stride * index(i);
            stride *= divisions_(i);
        }
        return data_[data_i];
    }
    const T& get(const Vectori<N>& index) const {
        return const_cast<GridData&>(*this).get(index);
    }

    T* get_if(const Vectori<N>& index) {
        if (!contains_index(index)) {
            return nullptr;
        }
        return &get(index);
    }
    const T* get_if(const Vectori<N>& index) const {
        return const_cast<GridData&>(*this).get_if(index);
    }

    bool contains_index(const Vectori<N>& index) const {
        if ((index.array() < 0).any() || (index.array() >= divisions_.array()).any()) {
            return false;
        }
        return true;
    }

    // Note: Query methods depend on whether data is stored at the cells
    // or corners

    GridPoint<Scalar, N> query(const Vector<Scalar, N>& normalized_position) const {
        Vector<Scalar, N> scaled = scaled_position(normalized_position);
        GridPoint<Scalar, N> result;
        result.a = scaled.template cast<int>(); // Rounds down on casting
        result.b = result.a.array() + 1;
        result.fraction = scaled.array() - scaled.array().floor();
        return result;
    }

    Vectori<N> query_closest(const Vector<Scalar, N>& normalized_position) const {
        Vector<Scalar, N> scaled = scaled_position(normalized_position);
        return scaled.array().round().template cast<int>();
    }

    Vector<Scalar, N> cell_size(const GridDimensions<Scalar, N>& dimensions) const {
        return dimensions.size.array() / divisions_;
    }

    const Vectori<N>& divisions() const { return divisions_; }

private:
    Vector<Scalar, N> scaled_position(const Vector<Scalar, N>& normalized_position) const {
        Vector<Scalar, N> result =
            divisions_.template cast<double>().array() * normalized_position.array();
        if (storage_type_ == GridStorageType::CELL) {
            result.array() -= 0.5;
        }
        return result;
    }

    Vectori<N> divisions_;
    std::vector<T> data_;
    GridStorageType storage_type_;
};

template <typename T>
using GridData2f = GridData<T, float, 2>;
template <typename T>
using GridData2d = GridData<T, double, 2>;
template <typename T>
using GridData3f = GridData<T, float, 3>;
template <typename T>
using GridData3d = GridData<T, double, 3>;

template <typename T, typename Scalar, int N>
class Grid {
public:

    Vector<Scalar, N>& origin() { return dimensions_.origin; }
    const Vector<Scalar, N>& origin() const { return dimensions_.origin; }
    Vector<Scalar, N>& size() { return dimensions_.size; }
    const Vector<Scalar, N>& size() const { return dimensions_.size; }

    void resize(const Vectori<N>& divisions, GridStorageType storage_type = GridStorageType::CELL) {
        data_.resize(divisions, storage_type);
    }

    T& get(const Vectori<N>& index) {
        return data_.get(index);
    }
    const T& get(const Vectori<N>& index) const {
        return data_.get(index);
    }

    T* get_if(const Vectori<N>& index) {
        return data_.get_if(index);
    }
    const T* get_if(const Vectori<N>& index) const {
        return data_.get_if(index);
    }

    bool contains_index(const Vectori<N>& index) const {
        return data_.contains_index(index);
    }

    bool contains(const Vector<Scalar, N>& position) const {
        return data_.contains_index(index);
    }

    GridPoint<Scalar, N> query(const Vector<Scalar, N>& position) {
        return data_.query(dimensions_.to_normalized(position));
    }

    Vectori<N> query_closest(const Vector<Scalar, N>& position) {
        return data_.query_closest(dimensions_.to_normalized(position));
    }

    Vector<Scalar, N> cell_size() const {
        return dimensions_.size_.array() / data_.divisions_;
    }

private:
    GridData<T, Scalar, N> data_;
    GridDimensions<Scalar, N> dimensions_;
};

template <typename T>
using Grid2f = Grid<T, float, 2>;
template <typename T>
using Grid2d = Grid<T, double, 2>;
template <typename T>
using Grid3f = Grid<T, float, 3>;
template <typename T>
using Grid3d = Grid<T, double, 3>;

} // namespace rbm
