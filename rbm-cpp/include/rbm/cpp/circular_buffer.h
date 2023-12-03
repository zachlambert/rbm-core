#pragma once

#include <vector>
#include <optional>


namespace rbm {

template <typename T>
class CircularBuffer {
public:
    CircularBuffer(size_t max_size = 100):
        tail(0),
        max_size_(max_size)
    {}

    void set_max_size(size_t max_size) {
        new_max_size = max_size;
    }
    size_t max_size()const {
        return max_size_;
    }

    void push_back(const T& value)
    {
        // If "new size" > "current size", immediately change the
        // size since we can continue growing
        if (new_max_size.has_value() && new_max_size.value() >= max_size_) {
            max_size_ = new_max_size.value();
            new_max_size = std::nullopt;
        }

        // If new_max_size still has a value, then it is less than the
        // current max_value. Wait until the buffer wraps around to
        // tail = new_max_size, then resize, discarding the back portion.
        if (new_max_size.has_value() && tail == new_max_size.value()) {
            data_.resize(new_max_size.value());
            max_size_ = new_max_size.value();
            new_max_size = std::nullopt;
            tail = 0;
        }

        if (max_size_ == 0) return;

        if (tail == 0 && data_.size() < max_size_) {
            // If the buffer size is less than max size, push to the end if
            // tail == 0, meaning we are at the end. If the buffer is resized
            // (to a larger size), then there will be some delay before it wraps
            // around to the end again.
            data_.push_back(value);

        } else {
            // Otherwise, set data[tail] and increment tail
            data_[tail] = value;
            tail = (tail + 1) % data_.size();
        }
    }

    // TODO: Create an iterator

    const T& operator[](size_t index)const
    {
        size_t i = (tail + index) % data_.size();
        return data_[i];
    }

    T& operator[](size_t index)
    {
        size_t i = (tail + index) % data_.size();
        return data_[i];
    }

    size_t size()const {
        return data_.size();
    }

    bool empty()const {
        return data_.size() == 0;
    }

    // Access raw data (used for plotting)

    const std::vector<T>& data()const {
        return data_;
    }
    std::vector<T>& data() {
        return data_;
    }
    size_t offset()const {
        return tail;
    }

private:
    std::vector<T> data_;
    size_t tail;
    size_t max_size_;
    std::optional<size_t> new_max_size;
};

} // namespace rbm
