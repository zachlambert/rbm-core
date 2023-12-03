#pragma once

#include "rbm/gui/renderer.h"


namespace rbm {

class Entity {
public:
    virtual bool init(Renderers& renderers) = 0;
    virtual void render(Renderers& renderers) = 0;
};

} // namespace rbm
