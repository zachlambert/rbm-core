#pragma once

#include <mbox/transform/transform.h>


namespace sviz {

class Camera {
public:
    Camera();
    void update_view(const mbox::Transform3d& pose);
    void update_projection(float aspect_ratio);

    const mbox::Matrix4f& get_view()const { return view; }
    const mbox::Matrix4f& get_projection()const { return projection; }

private:
    float clipping_near;
    float clipping_far;
    float fov_degrees;

    mbox::Matrix4f view;
    mbox::Matrix4f projection;

    friend class CameraController;
};

} // namespace sviz
