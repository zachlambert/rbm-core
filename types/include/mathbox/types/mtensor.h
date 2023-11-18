#pragma once

#include "mathbox/types/matrix.h"
#include "cpp_utils/darray.h"


namespace mbox {

// MTensor = Matrix tensor
// ie: Uses eigen matrix types
// Very lightweight, only really a convenience class to return a list of matrices.

template <typename Scalar, int Rows, int Cols, int Depth>
class MTensor {
    using matrix_t = Matrix<Scalar, Rows, Cols>;
public:
    MTensor():
        rows_(Rows == -1 ? 0 : Rows),
        cols_(Cols == -1 ? 0 : Cols),
        depth_(Depth == -1 ? 0 : Depth)
    {}
    MTensor(int rows, int cols, int depth):
        rows_(rows),
        cols_(cols),
        depth_(depth)
    {
        if constexpr(Rows != -1) {
            assert(rows_ == Rows);
        }
        if constexpr(Cols != -1) {
            assert(cols_ == Cols);
        }
        if constexpr(Depth != -1) {
            assert(depth == Depth);
        }

        if constexpr(Depth == -1) {
            data_.resize(depth);
        }
        for (auto& matrix: data_) {
            matrix.resize(rows, cols);
        }
    }
    void setZero()
    {
        for (auto& matrix: data_) {
            matrix.setZero();
        }
    }
    static MTensor Zero()
    {
        MTensor result;
        result.setZero();
        return result;
    }
    static MTensor Zero(int rows, int cols, int depth)
    {
        MTensor result(rows, cols, depth);
        result.setZero();
        return result;
    }

    const matrix_t& operator[](std::size_t index) const
    {
        return data_[index];
    }
    matrix_t& operator[](std::size_t index)
    {
        return data_[index];
    }

    MTensor<Scalar, Depth, Cols, Rows> reverse() const
    {
        MTensor<Scalar, Depth, Cols, Rows> reverse(depth_, cols_, rows_);
        for (std::size_t i = 0; i < rows_; i++) {
            for (std::size_t j = 0; j < cols_; j++) {
                for (std::size_t k = 0; k < depth_; k++) {
                    reverse[i](k, j) = (*this)[k](i, j);
                }
            }
        }
        return reverse;
    }

    friend matrix_t operator*(const MTensor& lhs, const Vector<Scalar, Depth>& rhs)
    {
        assert(rhs.rows() > 0);
        assert(rhs.rows() == lhs.depth());
        matrix_t result = matrix_t::Zero(lhs.rows(), lhs.cols());
        for (std::size_t i = 0; i < lhs.depth(); i++) {
            result += lhs[i] * rhs(i);
        }
        return result;
    }

    template <int OtherCols>
    friend MTensor<Scalar, Rows, Cols, OtherCols> operator*(const MTensor& lhs, const Matrix<Scalar, Depth, OtherCols>& rhs)
    {
        assert(rhs.rows() > 0);
        assert(rhs.cols() > 0);
        assert(rhs.rows() == lhs.depth());
        MTensor<Scalar, Rows, Cols, OtherCols> result(lhs.rows(), lhs.cols(), rhs.cols());
        for (std::size_t j = 0; j < result.depth(); j++) {
            result[j].setZero();
            for (std::size_t i = 0; i < lhs.depth(); i++) {
                result[j] += lhs[i] * rhs(i, j);
            }
        }
        return result;
    }

    template <int OtherRows>
    friend MTensor<Scalar, OtherRows, Cols, Depth> operator*(const Matrix<Scalar, OtherRows, Rows>& lhs, const MTensor& rhs)
    {
        assert(lhs.rows() > 0);
        assert(lhs.rows() > 0);
        assert(lhs.cols() == rhs.rows());
        MTensor<Scalar, OtherRows, Cols, Depth> result(lhs.rows(), rhs.cols(), rhs.depth());
        for (std::size_t k = 0; k < rhs.depth(); k++) {
            result[k].setZero();
            for (std::size_t i = 0; i < lhs.rows(); i++) {
                result[k].block(i, 0, 1, rhs.cols()) = lhs.block(i, 0, 1, lhs.cols()) * rhs[k];
            }
        }
        return result;
    }

    std::size_t rows() const { return rows_; }
    std::size_t cols() const { return cols_; }
    std::size_t depth() const { return depth_; }

private:
    std::size_t rows_;
    std::size_t cols_;
    std::size_t depth_;
    cpp_utils::darray<Matrix<Scalar, Rows, Cols>, Depth> data_;
};

template <int Rows, int Cols, int Depth>
using MTensorf = MTensor<float, Rows, Cols, Depth>;
template <int Rows, int Cols, int Depth>
using MTensord = MTensor<double, Rows, Cols, Depth>;

} // namespace mbox
