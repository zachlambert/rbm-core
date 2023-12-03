#pragma once

#include <string>

namespace rbm {

class Widget {
public:
    virtual bool init() { return true; }
    virtual void shutdown() {}
    virtual void render() {}
};

} // namespace rbm
