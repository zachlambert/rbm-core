#pragma once

#include "mathbox/transform/transform.h"

namespace mbox {

template <typename Scalar, int Dim>
class AdjointTransform {
public:
    typedef Eigen::Matrix<Scalar, dim_se_adj<Dim>(), dim_se_adj<Dim>()> matrix_t;
    typedef Transform<Scalar, Dim> transform_t;

    AdjointTransform(const transform_t& transform):
        transform(transform)
    {}
    AdjointTransform(const CompactTransform<Scalar, Dim>& transform):
        transform(transform.toTransform())
    {}
    AdjointTransform inverse()
    {
        return AdjointTransform(transform.inverse());
    }
    AdjointTransform& operator*=(const AdjointTransform& rhs)
    {
        transform = transform * rhs.transform;
    }

    matrix_t matrix()
    {
        matrix_t result;
        result.template block<Dim, Dim>(0, 0) = transform.rotation();
        result.template block<Dim, Dim>(0, Dim).setZero();
        result.template block<Dim, Dim>(Dim, 0) = cross_product_matrix(transform.translation().eval()) * transform.rotation();
        result.template block<Dim, Dim>(Dim, Dim) = transform.rotation();
        return result;
    }

    template <typename Scalar_>
    AdjointTransform<Scalar_, Dim> cast() const {
        return AdjointTransform<Scalar_, Dim>(transform.template cast<Scalar_>());
    }

    static AdjointTransform Identity() {
        return AdjointTransform(Transform<Scalar, Dim>::Identity());
    }

private:
    transform_t transform;
};
typedef AdjointTransform<float, 2> AdjointTransform2f;
typedef AdjointTransform<double, 2> AdjointTransform2d;
typedef AdjointTransform<float, 3> AdjointTransform3f;
typedef AdjointTransform<double, 3> AdjointTransform3d;

template <typename Scalar, int Dim>
AdjointTransform<Scalar, Dim> operator*(AdjointTransform<Scalar, Dim> lhs, const AdjointTransform<Scalar, Dim>& rhs)
{
    lhs *= rhs;
    return lhs;
}

template <typename Scalar, int Dim>
class AdjointLieBracket {
public:
    typedef Eigen::Matrix<Scalar, dim_se_adj<Dim>(), dim_se_adj<Dim>()> matrix_t;
    typedef Eigen::Vector<Scalar, dim_se_adj<Dim>()> coords_t;

    AdjointLieBracket(const coords_t& coords):
        coords(coords)
    {}
    matrix_t matrix()
    {
        matrix_t result;
        result.template block<Dim, Dim>(0, 0) = cross_product_matrix(coords.template head<dim_so<Dim>()>().eval());
        result.template block<Dim, Dim>(0, Dim).setZero();
        result.template block<Dim, Dim>(Dim, 0) = cross_product_matrix(coords.template tail<Dim>().eval());
        result.template block<Dim, Dim>(Dim, Dim) = cross_product_matrix(coords.template head<dim_so<Dim>()>().eval());
        return result;
    }

    template <typename Scalar_>
    AdjointLieBracket<Scalar_, Dim> cast() const {
        return AdjointLieBracket<Scalar_, Dim>(coords.template cast<Scalar_>());
    }
private:
    coords_t coords;
};
typedef AdjointLieBracket<float, 2> AdjointLieBracket2f;
typedef AdjointLieBracket<double, 2> AdjointLieBracket2d;
typedef AdjointLieBracket<float, 3> AdjointLieBracket3f;
typedef AdjointLieBracket<double, 3> AdjointLieBracket3d;

} // namespace mbox
