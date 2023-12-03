#pragma once

#include <unordered_map>
#include <memory>
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <rbm/types/color.h>
#include "rbm/gui/widget.h"
#include "rbm/gui/entity.h"
#include "rbm/gui/renderer.h"
#include "rbm/gui/camera_controller.h"


namespace rbm {

class Viewport: public Widget {
public:
    Viewport(int width, int height);

    bool init() override;
    void render() override;
    void add_entity(const std::string& name, const std::shared_ptr<Entity>& entity);

private:
    int width;
    int height;
    ColorRGBd bg_color;

    // Frame buffer data
    struct {
        unsigned int fbo;
        unsigned int rbo;
        unsigned int texture_id;
    } glData;

    Renderers renderers;
    std::unordered_map<std::string, std::shared_ptr<Entity>> entities;
    Camera camera;
    CameraController camera_controller;
};

} // namespace rbm
