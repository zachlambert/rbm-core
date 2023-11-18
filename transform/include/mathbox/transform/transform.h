#pragma once

#include <Eigen/Geometry>
#include "mathbox/transform/group_dimensions.h"
#include "mathbox/transform/angle.h"
#include "mathbox/transform/cross_product_matrix.h"
#include "mathbox/transform/screw_transform.h"
#include "mathbox/transform/rotation.h"
#include <cmath>


namespace mbox {

template <typename Scalar, int Dim>
class Transform {
    typedef Eigen::Matrix<Scalar, Dim + 1, Dim + 1> matrix_t;
    typedef Eigen::Vector<Scalar, Dim> translation_t;
    typedef Rotation<Scalar, Dim> rotation_t;
public:
    Transform() {
        matrix_.setIdentity();
    }
    Transform(const rotation_t& rotation_, const translation_t& translation_) {
        matrix_.setIdentity();
        rotation() = rotation_;
        translation() = translation_;
    }
    Transform(const Eigen::Transform<Scalar, Dim, Eigen::Isometry>& transform) {
        matrix_.setIdentity();
        rotation() = transform.rotation();
        translation() = transform.translation();
    }
    Transform& operator=(const Eigen::Transform<Scalar, Dim, Eigen::Isometry>& transform) {
        matrix_.setIdentity();
        rotation() = transform.rotation();
        translation() = transform.translation();
        return *this;
    }
private:
    Transform(const matrix_t& matrix): matrix_(matrix) {}
    template <typename Scalar_, int Dim_>
    friend class Transform;
public:

    Eigen::Block<matrix_t, Dim, 1> translation() { return matrix_.template block<Dim, 1>(0, Dim); }
    Eigen::Block<const matrix_t, Dim, 1> translation() const { return matrix_.template block<Dim, 1>(0, Dim); }
    Eigen::Block<matrix_t, Dim, Dim> rotation() { return matrix_.template block<Dim, Dim>(0, 0); }
    Eigen::Block<const matrix_t, Dim, Dim> rotation() const { return matrix_.template block<Dim, Dim>(0, 0); }

    const matrix_t& matrix() const { return matrix_; }

    static Transform Identity()
    {
        return Transform();
    }
    void setIdentity() {
        matrix_.setIdentity();
    }

    Transform inverse() const {
        Transform result;
        result.rotation() = rotation().transpose();
        result.translation() = -result.rotation() * translation();
        return result;
    }

    template <typename Scalar_>
    Transform<Scalar_, Dim> cast() const {
        return Transform<Scalar_, Dim>(matrix_.template cast<Scalar_>());
    }

    friend Eigen::Vector<Scalar, Dim> operator*(const Transform<Scalar, Dim>& lhs, const Eigen::Vector<Scalar, Dim>& rhs) {
        return lhs.translation() + lhs.rotation() * rhs;
    }

    Transform& operator*=(const Transform<Scalar, Dim>& rhs) {
        translation() += rotation() * rhs.translation();
        rotation() *= rhs.rotation();
        return *this;
    }

    friend Transform operator*(Transform lhs, const Transform& rhs) {
        lhs *= rhs;
        return lhs;
    }

private:
    matrix_t matrix_;
};
typedef Transform<float, 2> Transform2f;
typedef Transform<float, 3> Transform3f;
typedef Transform<double, 2> Transform2d;
typedef Transform<double, 3> Transform3d;
template <int Dim>
using Transformf = Transform<float, Dim>;
template <int Dim>
using Transformd = Transform<double, Dim>;


template <typename Scalar, int Dim>
class CompactTransform {
    typedef CompactRotation<Scalar, Dim> rotation_t;
    typedef Eigen::Vector<Scalar, Dim> translation_t;

public:
    CompactTransform() {} // Default initialised
    CompactTransform(const rotation_t& rotation, const translation_t& translation):
        rotation_(rotation), translation_(translation)
    {}

    rotation_t& rotation() { return rotation_; }
    const rotation_t& rotation() const { return rotation_; }
    translation_t& translation() { return translation_; }
    const translation_t& translation() const { return translation_; }

    static CompactTransform Identity()
    {
        CompactTransform result;
        result.rotation_.setIdentity();
        result.translation_.setZero();
        return result;
    }

    CompactTransform inverse() const {
        return CompactTransform(rotation_.inverse(), -(rotation_.inverse() * translation_));
    }

    template <typename Scalar_>
    CompactTransform<Scalar_, Dim> cast() const {
        return CompactTransform<Scalar_, Dim>(rotation_.template cast<Scalar_>(), translation_.template cast<Scalar_>());
    }

    Transform<Scalar, Dim> toTransform() const {
        return Transform<Scalar, Dim>(rotation_.toRotationMatrix(), translation_);
    }
    CompactTransform& operator=(const Transform<Scalar, Dim>& transform) {
        rotation_ = transform.rotation();
        translation_ = transform.translation();
        return *this;
    }

    friend Eigen::Vector<Scalar, Dim> operator*(const CompactTransform<Scalar, Dim>& lhs, const Eigen::Vector<Scalar, Dim>& rhs) {
        return lhs.translation_ + lhs.rotation_ * rhs;
    }

    CompactTransform& operator*=(const CompactTransform<Scalar, Dim>& rhs) {
        translation_ += rotation_ * rhs.translation();
        rotation_ *= rhs.rotation();
        return *this;
    }

    friend CompactTransform operator*(CompactTransform lhs, const CompactTransform& rhs) {
        lhs *= rhs;
        return lhs;
    }

private:
    rotation_t rotation_;
    translation_t translation_;
};
typedef CompactTransform<float, 2> CompactTransform2f;
typedef CompactTransform<float, 3> CompactTransform3f;
typedef CompactTransform<double, 2> CompactTransform2d;
typedef CompactTransform<double, 3> CompactTransform3d;


template <typename Scalar, int Dim>
class LogTransform {
public:
    static_assert(Dim == 2 || Dim == 3);
    typedef Eigen::Vector<Scalar, dim_se<Dim>()> coords_t;
    typedef Eigen::Matrix<Scalar, Dim + 1, Dim + 1> matrix_t;

    LogTransform(const coords_t& coords):
        coords_(coords)
    {}

    LogTransform(const Transform<Scalar, Dim>& transform) {
        coords_angular() = LogRotation<Scalar, Dim>(Rotation<Scalar, Dim>(transform.rotation())).coords();
        coords_linear() = screw_displacement_matrix<Scalar, Dim>(coords_angular()).inverse() * transform.translation();
    }
    Transform<Scalar, Dim> exp() const {
        Transform<Scalar, Dim> result;
        result.rotation() = LogRotation<Scalar, Dim>(coords_angular()).exp();
        result.translation() = screw_displacement_matrix<Scalar, Dim>(coords_angular()) * coords_linear();
        return result;
    }

    LogTransform(const CompactTransform<Scalar, Dim>& transform) {
        coords_angular() = LogRotation<Scalar, Dim>(transform.rotation()).coords();
        coords_linear() = screw_displacement_matrix<Scalar, Dim>(coords_angular()).inverse() * transform.translation();
    }
    CompactTransform<Scalar, Dim> exp_compact() const {
        CompactTransform<Scalar, Dim> result;
        result.rotation() = LogRotation<Scalar, Dim>(coords_angular()).exp_compact();
        result.translation() = screw_displacement_matrix<Scalar, Dim>(coords_angular()) * coords_linear();
        return result;
    }

    coords_t& coords() { return coords_; }
    const coords_t& coords()const { return coords_; }

    Eigen::VectorBlock<coords_t, dim_so<Dim>()> coords_angular() { return coords_.template head<dim_so<Dim>()>(); }
    Eigen::Vector<Scalar, dim_so<Dim>()> coords_angular() const { return coords_.template head<dim_so<Dim>()>(); }
    Eigen::VectorBlock<coords_t, Dim> coords_linear() { return coords_.template tail<Dim>(); }
    Eigen::Vector<Scalar, Dim> coords_linear() const { return coords_.template tail<Dim>(); }

    matrix_t matrix() {
        matrix_t matrix;
        // [x] = |[ang] lin |
        //       | 0^T   0  |
        matrix.template block<Dim, Dim>(0, 0) = cross_product_matrix(coords_angular()); // = LogRotation(coords_angular()).matrix()
        matrix.template block<Dim, 1>(0, Dim) = coords_linear();
        matrix.template block<1, Dim + 1>(Dim, 0).setZero();
        return matrix;
    }
private:
    coords_t coords_;
};
typedef LogTransform<float, 2> LogTransform2f;
typedef LogTransform<double, 2> LogTransform2d;
typedef LogTransform<float, 3> LogTransform3f;
typedef LogTransform<double, 3> LogTransform3d;

} // namespace mbox
