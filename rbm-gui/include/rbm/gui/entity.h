#pragma once

#include "sviz/render/renderer.h"

namespace sviz {

class Entity {
public:
    virtual bool init(Renderers& renderers) = 0;
    virtual void render(Renderers& renderers) = 0;
};

} // namespace sviz
