#include "rbm/gui/viewport.h"

#include <imgui.h>
#include <iostream>
#include "rbm/gui/marker_renderer.h"


namespace rbm {

Viewport::Viewport(int width, int height):
    width(width),
    height(height),
    bg_color(mbox::ColorRGBd(1.0, 1.0, 1.0))
{}

bool Viewport::init()
{
    // Generate and bind frame buffer
    glGenFramebuffers(1, &glData.fbo);
    glBindFramebuffer(GL_FRAMEBUFFER, glData.fbo);

    // Generate texture
    glGenTextures(1, &glData.texture_id);
    glBindTexture(GL_TEXTURE_2D, glData.texture_id);
    glTexImage2D(
        GL_TEXTURE_2D, 0, GL_RGB,
        width, height, 0, GL_RGB,
        GL_UNSIGNED_BYTE, NULL);

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glBindTexture(glData.texture_id, 0);

    // Bind texture to fbo
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, glData.texture_id, 0);

    // Generate a render buffer object
    glGenRenderbuffers(1, &glData.rbo);
    glBindRenderbuffer(GL_RENDERBUFFER, glData.rbo);
    glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH24_STENCIL8, width, height);
    glBindRenderbuffer(GL_RENDERBUFFER, 0);

    // Bind rbo to fbo
    glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_STENCIL_ATTACHMENT, GL_RENDERBUFFER, glData.rbo);

    // Verify that the frame buffer is complete
    if(glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE) {
        return false;
    }

    // Bind default frame buffer again
    glBindFramebuffer(GL_FRAMEBUFFER, 0);

    for (auto& [name, entity]: entities) {
        if (!entity->init(renderers)) {
            std::cerr << "Failed to initialise entity '" << name << "'" << std::endl;
            return false;
        }
    }

    mbox::Transform3d initial_camera_pose = mbox::EulerTransform3d(
        mbox::EulerRotation3d(0, 0.9, 0),
        mbox::Vector3d(-10, 0, 10)).toTransform();
    camera_controller.init(initial_camera_pose, camera);

    renderers.get<MarkerRenderer>();

    return true;
}

void Viewport::render()
{
    // Enable frame buffer

    glBindFramebuffer(GL_FRAMEBUFFER, glData.fbo);

    glClearColor(bg_color.r, bg_color.g, bg_color.b, 1);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_CULL_FACE);

    glViewport(0, 0, width, height);
    float aspect_ratio = (float)width / height;
    camera.update_projection(aspect_ratio);

    // Render items

    {
        mbox::Vector3d click_pos;
        if (camera_controller.is_clicked(click_pos)) {
            mbox::Cylinder3d cylinder;
            cylinder.pose.translation() = click_pos;
            cylinder.length = 0.01;
            cylinder.radius = 0.12;
            renderers.get<MarkerRenderer>()->queue_primitive(cylinder, mbox::ColorRGBd::White());
        }
    }

    for (auto& [_, entity]: entities) {
        entity->render(renderers);
    }
    renderers.render(camera.get_view(), camera.get_projection());

    // Disable frame buffer

    glBindFramebuffer(GL_FRAMEBUFFER, 0);

    // Render image to imgui window

    ImGui::BeginChild("##frame");

    // Must call this within the imgui frame
    camera_controller.update_pose(camera);

    ImVec2 viewport_size = ImGui::GetWindowSize();
    ImVec2 uv1, uv2;
    if (viewport_size.x > aspect_ratio*viewport_size.y) {
        uv1.x = 0;
        uv2.x = 1;
        float yScale = aspect_ratio*viewport_size.y/viewport_size.x;
        uv1.y = 0.5*(1 + yScale);
        uv2.y = 0.5*(1 - yScale);
    } else {
        uv1.y = 1;
        uv2.y = 0;
        float xScale = (float)viewport_size.x/(aspect_ratio*viewport_size.y);
        uv1.x = 0.5*(1 - xScale);
        uv2.x = 0.5*(1 + xScale);
    }

    ImGui::Image((ImTextureID)glData.texture_id, viewport_size, uv1, uv2);
    ImGui::EndChild();
}

void Viewport::add_entity(const std::string& name, const std::shared_ptr<Entity>& entity)
{
    entities.emplace(name, entity);
}

} // namespace rbm
