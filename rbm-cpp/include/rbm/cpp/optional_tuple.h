#pragma once

#include <tuple>
#include <exception>
#include <optional>


namespace rbm {

template <typename In, typename Out, std::size_t InSize>
struct tuple_range_iterator {
    static void set(const In& in, Out& out) {
        static_assert(InSize > 0);
        std::get<InSize-1>(out) = std::get<InSize-1>(in);
        tuple_range_iterator<In, Out, InSize-1>::set(in, out);
    }
};

template <typename In, typename Out>
struct tuple_range_iterator<In, Out, 0> {
    static void set(const In& in, Out& out) {}
};

template <typename ...T>
struct OptionalTuple {
    typedef std::tuple<T...> Tuple;
    static constexpr std::size_t capacity_ = std::tuple_size_v<Tuple>;
public:
    OptionalTuple():
        size_(0)
    {}

    template <std::size_t Index>
    std::tuple_element_t<Index, Tuple>& value() {
        static_assert(Index >= 0 && Index < capacity_);
        if (Index >= size_) {
            throw std::bad_optional_access();
        }
        return std::get<Index>(tuple_);
    }

    template <std::size_t Index>
    const std::tuple_element_t<Index, Tuple>& value() const {
        static_assert(Index < capacity_);
        if (Index >= size_) {
            throw std::bad_optional_access();
        }
        return std::get<Index>(tuple_);
    }

    template <std::size_t Index>
    void set(const std::tuple_element_t<Index, Tuple>& value) {
        std::size_t new_size = Index + 1;
        if (size_ + 1 < new_size) {
            throw std::bad_optional_access();
        }
        size_ = std::max(size_, new_size);
        std::get<Index>(tuple_) = value;
    }

    template <std::size_t InSize, typename ...Args>
    void set_range(const std::tuple<Args...>& tuple_in) {
        typedef std::tuple<Args...> TupleIn;
        static_assert(InSize <= std::tuple_size_v<TupleIn>);
        static_assert(InSize <= capacity_);

        tuple_range_iterator<TupleIn, Tuple, InSize>::set(tuple_in, tuple_);
        size_ = std::max(InSize, size_);
    }

    template <std::size_t InSize, typename ...Args>
    void set_range(const OptionalTuple<Args...>& tuple_in) {
        typedef OptionalTuple<Args...> TupleIn;
        static_assert(InSize <= TupleIn::capacity_);
        static_assert(InSize <= capacity_);

        if (InSize > tuple_in.size_) {
            throw std::bad_optional_access();
        }

        tuple_range_iterator<typename TupleIn::Tuple, Tuple, InSize>::set(tuple_in.tuple_, tuple_);
        size_ = std::max(InSize, size_);
    }

    const Tuple& tuple() const {
        if (size_ != capacity_) {
            throw std::bad_optional_access();
        }
        return tuple_;
    }

    template <typename ...Args>
    OptionalTuple<Args...> subset() const {
        typedef OptionalTuple<Args...> TupleOut;
        if (TupleOut::capacity_ > size_) {
            throw std::bad_optional_access();
        }
        TupleOut out;
        out.template set_range<TupleOut::capacity_>(*this);
        return out;
    }

    template <std::size_t Index>
    bool has_value() const {
        static_assert(Index < capacity_);
        return Index < size_;
    }

    int size() const {
        return size_;
    }

    template <std::size_t NewSize>
    void resize() {
        static_assert(NewSize <= capacity_);
        size_ = NewSize;
    }

    static constexpr std::size_t capacity() {
        return capacity_;
    }

private:
    Tuple tuple_;
    std::size_t size_;

    template <typename ...T_>
    friend class OptionalTuple;
};

} // namespace rbm
