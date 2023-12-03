#pragma once

namespace rbm {

inline bool flag_or(bool& flag, bool condition) {
    flag |= condition;
    return condition;
}

inline bool flag_and(bool& flag, bool condition) {
    flag &= condition;
    return condition;
}

} // namespace rbm
