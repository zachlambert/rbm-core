#pragma once

#include <chrono>
#include <mbox/transform/euler.h>
#include "sviz/render/camera.h"


namespace sviz {

class CameraController {
public:
    CameraController();

    void init(const mbox::EulerTransform3d& start_pose, Camera& camera);
    void update_pose(Camera& camera);
    bool is_clicked(mbox::Vector3d& click_pos);

private:
    mbox::Transform3d pose;

    enum class State {
        INACTIVE,
        INVALID,
        PAN,
        ROTATE,
        SCROLL
    } state;

    mbox::Vector3d click_pos_cs;
    mbox::Vector3d click_pos_ws;
    mbox::Transform3d prev_camera_pose;

    float scroll_amount;
    typedef std::chrono::high_resolution_clock clock_t;
    clock_t::time_point scroll_start;
};

} // namespace sviz
