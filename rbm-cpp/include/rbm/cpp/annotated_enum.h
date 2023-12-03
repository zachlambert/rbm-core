#pragma once

#include <concepts>
#include <string_view>

namespace rbm {

template <typename T>
struct enum_details {};

template <typename T>
constexpr size_t enum_count = enum_details<T>::count;

template <typename T>
T enum_get(size_t index) { return enum_details<T>::get(index); }

template <typename T>
std::size_t enum_index(T v) { return enum_details<T>::index(v); }

template <typename T>
const char* enum_label(std::size_t index) { return enum_details<T>::label(index); }

template <typename T>
const char* enum_label(T v) { return enum_details<T>::label(enum_details<T>::index(v)); }

template <typename T>
concept annotated_enum = requires(T v, size_t i) {
    std::is_enum<T>();
    { enum_details<T>::count } -> std::convertible_to<std::size_t>;
    { enum_details<T>::get(i) } -> std::convertible_to<T>;
    { enum_details<T>::index(v) } -> std::convertible_to<std::size_t>;
    { enum_details<T>::enum_label(i) } -> std::convertible_to<const char*>;
};

} // namespace rbm
