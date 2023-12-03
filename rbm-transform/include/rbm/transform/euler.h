#pragma once

#include "rbm/transform/transform.h"
#include "rbm/transform/angle.h"


namespace rbm {

template <typename Scalar, int Dim>
class EulerRotation {
    typedef Vector<Scalar, dim_so<Dim>()> coords_t;
public:
    EulerRotation() {} // Default initialised
    EulerRotation(const coords_t& coords):
        coords_(coords)
    {}
    EulerRotation(Scalar yaw):
        coords_(yaw)
    {
        static_assert(Dim == 2);
    }
    EulerRotation(Scalar roll, Scalar pitch, Scalar yaw):
        coords_(roll, pitch, yaw)
    {
        static_assert(Dim == 3);
    }

    coords_t& coords() { return coords_; }
    const coords_t& coords() const { return coords_; }

    Scalar& roll() {
        static_assert(Dim == 3);
        return coords_.x();
    }
    const Scalar& roll() const {
        static_assert(Dim == 3);
        return coords_.x();
    }
    Scalar& pitch() {
        static_assert(Dim == 3);
        return coords_.y();
    }
    const Scalar& pitch() const {
        static_assert(Dim == 3);
        return coords_.y();
    }
    Scalar& yaw() {
        if constexpr (Dim == 2) {
            return coords_(0);
        }
        if constexpr (Dim == 3) {
            return coords_.z();
        }
    }
    const Scalar& yaw() const {
        if constexpr (Dim == 2) {
            return coords_(0);
        }
        if constexpr (Dim == 3) {
            return coords_.z();
        }
    }

    void setIdentity() {
        coords_.setZero();
    }
    static EulerRotation Identity() {
        return EulerRotation(coords_t::Zero().eval());
    }
    template <typename Scalar_>
    EulerRotation<Scalar_, Dim> cast() const {
        return EulerRotation<Scalar_, Dim>(coords_.template cast<Scalar_>());
    }

    EulerRotation(const CompactRotation<Scalar, Dim>& q) {
        (*this) = q;
    }
    EulerRotation& operator=(const CompactRotation<Scalar, Dim>& q) {
        if constexpr (Dim == 2) {
            yaw() = q;
        }
        if constexpr (Dim == 3) {
            roll() = atan2(2*(q.w()*q.x() + q.y()*q.z()), 1 - 2*(pow(q.x(), 2) + pow(q.y(), 2)));
            pitch() = asin(2*(q.w()*q.y() - q.x()*q.z()));
            yaw() = atan2(2*(q.w()*q.z() + q.x()*q.y()), 1 - 2*(pow(q.y(), 2) + pow(q.z(), 2)));
        }
        return *this;
    }
    CompactRotation<Scalar, Dim> toCompact() const {
        if constexpr (Dim == 2) {
            return Angle<Scalar>(yaw());
        }
        if constexpr (Dim == 3) {
            typedef Vector<Scalar, 3> Axis;
            typedef Eigen::AngleAxis<Scalar> AngleAxis;
            Eigen::Quaternion<Scalar> q;
            q =
                AngleAxis(yaw(), Axis::UnitZ())
                * AngleAxis(pitch(), Axis::UnitY())
                * AngleAxis(roll(), Axis::UnitX());
            return q;
        }
    }

    EulerRotation(const Matrix<Scalar, Dim, Dim>& R) {
        (*this) = R;
    }
    EulerRotation& operator=(const Matrix<Scalar, Dim, Dim>& R) {
        if constexpr (Dim == 3) {
            pitch() = std::atan2(-R(2, 0), std::hypot(R(0, 0), R(1, 0)));
            roll() = std::atan2(R(2, 1), R(2, 2));
        }
        yaw() = std::atan2(R(1, 0), R(0, 0));
        return *this;
    }
    Matrix<Scalar, 3, 3> toRotationMatrix() const {
        return toCompact().toRotationMatrix();
    }

    void normalize() {
        for (std::size_t i = 0; i < coords_.size(); i++) {
            coords_(i) = Angle<Scalar>(coords_(i));
        }
        if constexpr (Dim == 3) {
            if (std::cos(pitch()) > 0)  return;

            std::array<Angle<Scalar>, Dim> angles;
            for (std::size_t i = 0; i < Dim; i++) angles[i] = coords_(i);
            angles[0] += M_PI;
            angles[1] = Angle<Scalar>(M_PI) - angles[1];
            angles[2] += M_PI;
            for (std::size_t i = 0; i < Dim; i++) coords_(i) = angles[i];
        }
    }

private:
    Vector<Scalar, Dim> coords_;
};
typedef EulerRotation<float, 2> EulerRotation2f;
typedef EulerRotation<double, 2> EulerRotation2d;
typedef EulerRotation<float, 3> EulerRotation3f;
typedef EulerRotation<double, 3> EulerRotation3d;


template <typename Scalar, int Dim>
class EulerTransform {
    typedef EulerRotation<Scalar, Dim> rotation_t;
    typedef Vector<Scalar, Dim> translation_t;
public:
    EulerTransform():
        rotation_(rotation_t::Identity()),
        translation_(translation_t::Zero())
    {}
    EulerTransform(const rotation_t& rotation, const translation_t& translation):
        rotation_(rotation),
        translation_(translation)
    {}

    rotation_t& rotation() { return rotation_; }
    const rotation_t& rotation() const { return rotation_; }
    translation_t& translation() { return translation_; }
    const translation_t& translation() const { return translation_; }

    EulerTransform(const CompactTransform<Scalar, Dim>& transform):
        rotation_(transform.rotation()),
        translation_(transform.translation())
    {}
    EulerTransform& operator=(const CompactTransform<Scalar, Dim>& transform) {
        rotation() = transform.rotation();
        translation() = transform.translation();
        return *this;
    }
    CompactTransform<Scalar, Dim> toCompact() const {
        return CompactTransform<Scalar, Dim>(rotation(), translation());
    }

    EulerTransform(const Transform<Scalar, Dim>& transform):
        rotation_(transform.rotation().eval()),
        translation_(transform.translation())
    {}
    EulerTransform& operator=(const Transform<Scalar, Dim>& transform) {
        rotation() = transform.rotation();
        translation() = transform.translation();
        return *this;
    }
    Transform<Scalar, Dim> toTransform() const {
        Transform<Scalar, Dim> transform;
        transform.translation() = translation();
        transform.rotation() = rotation().toRotationMatrix();
        return transform;
    }

    void setIdentity() {
        rotation_.setIdentity();
        translation_.setZero();
    }
    static EulerTransform Identity() {
        return EulerTransform(rotation_t::Identity(), translation_t::Zero());
    }

    template <typename Scalar_>
    EulerTransform<Scalar_, Dim> cast() const {
        return EulerTransform<Scalar_, Dim>(rotation_.template cast<Scalar_>(), translation_.template cast<Scalar_>());
    }

private:
    rotation_t rotation_;
    translation_t translation_;
};
typedef EulerTransform<float, 2> EulerTransform2f;
typedef EulerTransform<float, 3> EulerTransform3f;
typedef EulerTransform<double, 2> EulerTransform2d;
typedef EulerTransform<double, 3> EulerTransform3d;

} // namespace rbm
