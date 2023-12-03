#pragma once

#include <vector>
#include <array>
#include <type_traits>
#include <assert.h>


namespace rbm {

template <typename T, int Dim, typename Allocator = std::allocator<T>>
using darray = std::conditional_t<Dim != -1, std::array<T, Dim>, std::vector<T, Allocator>>;

template <typename T, int Dim>
void resize_darray(darray<T, Dim>& darray, std::size_t size) {
    if constexpr (Dim != -1) {
        assert(size == darray.size());
    }
    if constexpr (Dim == -1) {
        darray.resize(size);
    }
}

} // namespace rbm
