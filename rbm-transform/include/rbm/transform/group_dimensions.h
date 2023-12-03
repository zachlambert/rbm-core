#pragma once


namespace rbm {

template <int Dim>
constexpr int dim_so() {
    // Note: Can probably generalise this
    static_assert(Dim == 2 || Dim == 3);
    return (Dim == 2 ? 1 : 3);
}

template <int Dim>
constexpr int dim_se() {
    return dim_so<Dim>() + Dim;
}

template <int Dim>
constexpr int dim_se_adj() {
    return Dim * 2;
}

} // namespace rbm
