#pragma once

#include "mathbox/types/matrix.h"


namespace mbox {

template <typename T, int Dim>
Matrix<T, Dim+1, Dim+1> make_homogeneous_translation(const Vector<T, Dim>& translation) {
    Matrix<T, Dim+1, Dim+1> result;
    result.setIdentity();
    result.template block<Dim, 1>(0, Dim) = translation;
    return result;
}

template <typename T, int Dim>
Matrix<T, Dim+1, Dim+1> make_homogeneous_scaling(const Vector<T, Dim>& scaling) {
    Matrix<T, Dim+1, Dim+1> result;
    result.setIdentity();
    for (int i = 0; i < Dim; i++) result(i, i) = scaling(i);
    return result;
}

template <typename T, int Dim>
Matrix<T, Dim+1, Dim+1> make_homogeneous_scaling(T scaling) {
    Matrix<T, Dim+1, Dim+1> result;
    result.setIdentity();
    for (int i = 0; i < Dim; i++) result(i, i) = scaling;
    return result;
}

} // namespace mbox
