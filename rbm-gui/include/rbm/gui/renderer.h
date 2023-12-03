#pragma once

#include <mbox/types/matrix.h>
#include <memory>
#include <type_traits>
// Not required in this header, but required for any rendering code, so
// conventient to leave here
#include <GL/glew.h>


namespace sviz {

class Renderer {
public:
    virtual void render(const mbox::Matrix4f& view, const mbox::Matrix4f& projection) = 0;
};

class Renderers {
public:
    template <typename T>
    std::shared_ptr<T> get()
    {
        static_assert(std::is_base_of_v<Renderer, T>);
        for (auto& renderer: renderers) {
            if (auto downcasted = std::dynamic_pointer_cast<T>(renderer)) {
                return downcasted;
            }
        }
        auto renderer = std::make_shared<T>();
        renderers.push_back(renderer);
        return renderer;
    }

    void render(const mbox::Matrix4f& view, const mbox::Matrix4f& projection)
    {
        for (auto& renderer: renderers) {
            renderer->render(view, projection);
        }
    }

private:
    std::vector<std::shared_ptr<Renderer>> renderers;
};

} // namespace sviz
