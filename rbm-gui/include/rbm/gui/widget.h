#pragma once

#include <string>

namespace sviz {

class Widget {
public:
    virtual bool init() { return true; }
    virtual void shutdown() {}
    virtual void render() {}
};

} // namespace sviz
