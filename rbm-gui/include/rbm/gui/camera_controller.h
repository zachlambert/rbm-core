#pragma once

#include <chrono>
#include <rbm/transform/euler.h>
#include "rbm/gui/camera.h"


namespace rbm {

class CameraController {
public:
    CameraController();

    void init(const EulerTransform3d& start_pose, Camera& camera);
    void update_pose(Camera& camera);
    bool is_clicked(Vector3d& click_pos);

private:
    Transform3d pose;

    enum class State {
        INACTIVE,
        INVALID,
        PAN,
        ROTATE,
        SCROLL
    } state;

    Vector3d click_pos_cs;
    Vector3d click_pos_ws;
    Transform3d prev_camera_pose;

    float scroll_amount;
    typedef std::chrono::high_resolution_clock clock_t;
    clock_t::time_point scroll_start;
};

} // namespace rbm
