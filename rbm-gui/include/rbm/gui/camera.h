#pragma once

#include <rbm/transform/transform.h>


namespace rbm {

class Camera {
public:
    Camera();
    void update_view(const Transform3d& pose);
    void update_projection(float aspect_ratio);

    const Matrix4f& get_view()const { return view; }
    const Matrix4f& get_projection()const { return projection; }

private:
    float clipping_near;
    float clipping_far;
    float fov_degrees;

    Matrix4f view;
    Matrix4f projection;

    friend class CameraController;
};

} // namespace rbm
