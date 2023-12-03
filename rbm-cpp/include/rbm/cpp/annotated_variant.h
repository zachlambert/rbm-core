#pragma once

#include <concepts>
#include <string_view>


namespace rbm {

template <typename T>
struct variant_details {};

template <typename T>
constexpr size_t variant_count = variant_details<T>::count;

template <typename T>
T variant_construct(std::size_t index) { return variant_details<T>::construct(index); }

template <typename T>
std::size_t variant_index(const T& v) { return variant_details<T>::index(v); }

template <typename T>
std::string_view variant_label(std::size_t index) { return variant_details<T>::label(index); }

template <typename T>
std::string_view variant_label(const T& v) { return variant_details<T>::label(variant_details<T>::index(v)); }

template <typename T>
concept annotated_variant = requires(const T& v, std::size_t i) {
    { variant_details<T>::count } -> std::convertible_to<std::size_t>;
    { variant_details<T>::construct(i) } -> std::convertible_to<T>;
    { variant_details<T>::index(v) } -> std::convertible_to<std::size_t>;
    { variant_details<T>::label(i) } -> std::convertible_to<std::string_view>;
};

} // namespace rbm
