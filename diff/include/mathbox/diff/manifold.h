#pragma once

#include "mathbox/types/matrix.h"
#include "mathbox/transform/transform.h"
#include "mathbox/transform/angle.h"
#include <type_traits>


namespace mbox {

template <typename T>
struct manifold_details {};


template <typename T>
static constexpr int manifold_dim = manifold_details<T>::dim;
template <typename T>
using manifold_delta = Vectord<manifold_dim<T>>;

template <typename T>
using manifold_hessian = Matrixd<manifold_dim<T>, manifold_dim<T>>;

template <typename T>
int manifold_dynamic_dim(const T& x) {
    if constexpr (manifold_dim<T> == Eigen::Dynamic) {
        return manifold_details<T>::dynamic_dim(x);
    }
    if constexpr (manifold_dim<T> != Eigen::Dynamic) {
        return manifold_details<T>::dim;
    }
}

template <typename T>
manifold_delta<T> zero_manifold_delta(const T& x) {
    if constexpr(manifold_dim<T> == -1) {
        int dim = manifold_dynamic_dim(x);
        return VectorXd::Zero(dim);
    }
    if constexpr(manifold_dim<T> != -1) {
        return manifold_delta<T>::Zero();
    }
}

template <typename T>
manifold_hessian<T> zero_manifold_hessian(const T& x) {
    if constexpr(manifold_dim<T> == -1) {
        int dim = manifold_dynamic_dim(x);
        return MatrixXd::Zero(dim, dim);
    }
    if constexpr(manifold_dim<T> != -1) {
        return manifold_hessian<T>::Zero();
    }
}

template <typename T>
void manifold_add(T& t, const manifold_delta<T>& delta) {
    manifold_details<T>::add(t, delta);
}
template <typename T>
T manifold_sum(T t, const manifold_delta<T>& delta) {
    manifold_add(t, delta);
    return t;
}

template <typename T>
manifold_delta<T> manifold_difference(const T& a, const T& b) {
    return manifold_details<T>::difference(a, b);
}

template <>
struct manifold_details<double> {
    static constexpr int dim = 1;
    using coord_type = double;
    static int dynamic_dim(double) {
        return 1;
    }
    static void add(double& x, const Vector1d& delta) {
        x += delta.value();
    }
    static Vector1d difference(double a, double b) {
        Vector1d delta;
        delta(0) = b - a;
        return delta;
    }
    static double zero() {
        return 0;
    }
};

template <int N>
struct manifold_details<Vectord<N>> {
    static constexpr int dim = N;
    using coord_type = Vectord<N>;
    static int dynamic_dim(const Vectord<N>& x) {
        return x.size(); // May or may not be dynamic
    }
    static void add(Vectord<N>& x, const Vectord<N>& delta) {
        x += delta;
    }
    static Vectord<N> difference(const Vectord<N>& a, const Vectord<N>& b) {
        return b - a;
    }
    static Vectord<N> zero() {
        return Vectord<N>::Zero();
    }
};

template <int Rows, int Cols>
struct manifold_details<Matrixd<Rows, Cols>> {
    static constexpr int dim = Rows + Cols;
    using coord_type = Vectord<dim>;
    static void add(Matrixd<Rows, Cols>& x, const Vectord<dim>& delta) {
        x.reshaped() += delta;
    }
    static Vectord<dim> difference(const Matrixd<Rows, Cols>& a, const Matrixd<Rows, Cols>& b) {
        return b.reshaped() - a.reshaped();
    }
    static Matrixd<Rows, Cols> zero() {
        return Matrixd<Rows, Cols>::Zero();
    }
};

template <int Dim>
struct manifold_details<Rotationd<Dim>> {
    static constexpr int dim = dim_so<Dim>();
    using coord_type = Vectord<dim>;
    static void add(Rotationd<Dim>& x, const Vectord<dim>& delta) {
        x = x * LogRotation<double, Dim>(delta).exp();
    }
    static Vectord<dim> difference(const Rotationd<Dim>& a, const Rotationd<Dim>& b) {
        return LogRotation<double, Dim>(a.transpose() * b).coords();
    }
    static Rotationd<Dim> zero() {
        return Rotationd<Dim>::Identity();
    }
};

template <int Dim>
struct manifold_details<Transformd<Dim>> {
    static constexpr int dim = dim_se<Dim>();
    using coord_type = Vectord<dim>;
    static void add(Transformd<Dim>& x, const Vectord<dim>& delta) {
        x = x * LogTransform<double, Dim>(delta).exp();
    }
    static Vectord<dim> difference(const Transformd<Dim>& a, const Transformd<Dim>& b) {
        return LogTransform<double, Dim>(a.inverse() * b).coords();
    }
    static Transformd<Dim> zero() {
        return Transformd<Dim>::Identity();
    }
};

template <>
struct manifold_details<Angled> {
    static constexpr int dim = 1;
    using coord_type = Vectord<dim>;
    static void add(Angled& x, const Vectord<dim>& delta) {
        x = x + delta(0);
    }
    static Vectord<dim> difference(const Angled& a, const Angled& b) {
        return Vector1d(b - a);
    }
};

struct Empty {};
template <>
struct manifold_details<Empty> {
    static constexpr int dim = 0;
    static int dynamic_dim(const Empty& x) {
        return 0;
    }
    static void add(Empty& x, const Vectord<0>& delta) {}
    static Vectord<0> difference(const Empty& a, const Empty& b) { return Vectord<0>(); }
};

} // namespace mbox
