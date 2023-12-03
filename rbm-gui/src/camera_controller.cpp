#include "rbm/gui/camera_controller.h"
#include <imgui.h>
#include <rbm/transform/miscellaneous.h>
#include <rbm/transform/euler.h>
#include <rbm/transform/radians.h>


namespace rbm {

CameraController::CameraController() {}

void CameraController::init(const EulerTransform3d& start_pose, Camera& camera)
{
    pose = start_pose.toTransform();
    camera.update_view(pose);
}

void CameraController::update_pose(Camera& camera)
{
    State current_state = State::INACTIVE;

    if (ImGui::IsWindowHovered()) {
        if (ImGui::IsMouseDown(ImGuiMouseButton_Left)) {
            current_state = State::PAN;
        } else if (ImGui::IsMouseDown(ImGuiMouseButton_Right)) {
            current_state = State::ROTATE;
        } else {
            float scroll = ImGui::GetIO().MouseWheel;
            if (scroll != 0) {
                current_state = State::SCROLL;
                scroll_amount += scroll;
                scroll_start = clock_t::now();
            } else if (state == State::SCROLL) {
                if ((double)std::chrono::duration_cast<std::chrono::nanoseconds>(clock_t::now() - scroll_start).count() * 1e-9 < 0.1) {
                    current_state = State::SCROLL;
                }
            }
        }
    }

    Vector2d pixel_coords;
    {
        ImVec2 window_size = ImGui::GetWindowSize();
        ImVec2 window_pos = ImGui::GetWindowPos();
        ImVec2 click_pos = ImGui::GetMousePos();

        pixel_coords.x() = -1.0 + 2*((double)click_pos.x - (double)window_pos.x) / window_size.x;
        pixel_coords.y() = -1.0 + 2*((double)click_pos.y - (double)window_pos.y) / window_size.y;
        if (window_size.x > window_size.y) {
            pixel_coords.y() *= (double)window_size.y / window_size.x;
        } else {
            pixel_coords.x() *= (double)window_size.x / window_size.y;
        }
    }
    if (abs(pixel_coords.x()) > 1 || abs(pixel_coords.y()) > 1) {
        current_state = State::INACTIVE;
    }

    if (current_state == State::INACTIVE) {
        state = State::INACTIVE;
        scroll_amount = 0;
        return;
    }

    Vector3d dir_cs; // Direction cameraspace
    {
        double focal_length = 1 / tan(0.5 * to_radians(camera.fov_degrees));
        dir_cs.y() = -pixel_coords.x();
        dir_cs.z() = -pixel_coords.y();
        dir_cs.x() = focal_length;
        dir_cs.normalize();
    }

    if (state == State::INACTIVE) {
        // New click, so calculate the position of the click in the world
        // and save the current camera pose.

        // Project the ray from the camera onto the world
        // For now, just project onto plane at z = 0. Perhaps can project onto geometry later.
        Vector3d pos_ws = pose.translation();
        Vector3d dir_ws = pose.rotation() * dir_cs;

        double distance = (0 - pos_ws.z()) / dir_ws.z();

        if (distance > 100 || distance < 0) {
            // Not looking towards the plane significantly, we get a large distance
            // so avoid moving the camera
            state = State::INVALID;
        } else {
            prev_camera_pose = pose;
            state = current_state;

            // Store the cameraspace position for calculation of camera motion
            click_pos_cs = dir_cs * distance;
            // Store the worldspace position to display a marker
            click_pos_ws = pos_ws + dir_ws * distance;
        }
    }
    else if (state == State::PAN)
    {
        pose = prev_camera_pose;
        pose.translation() += pose.rotation() * (click_pos_cs - dir_cs * (click_pos_cs.x() / dir_cs.x()));
    }
    else if (state == State::ROTATE)
    {
        pose = prev_camera_pose;
        EulerRotation3d euler(pose.rotation().eval());

        // Make sure roll is kept at zero
        euler.roll() = 0;

        Eigen::AngleAxisd pitch_only = Eigen::AngleAxisd(euler.pitch(), Vector3d::UnitY());

        Vector3d from = pitch_only * click_pos_cs;
        Vector3d to = pitch_only * dir_cs;
        // Solve:
        // Rz(yaw_change) * Ry(pitch_change) * to = from
        // (scale invariant)

        // 1. [Ry(pitch_change) * to] defines a cone about the y axis.
        // Find the Rz(yaw_change) that rotates the line defined by from
        // onto that cone.
        double yaw_change;
        {
            double cone_angle = asin(hypot(to.x(), to.z()) / to.norm());
            double inclination = asin(from.z() / from.norm());
            double angle_to_line_on_cone_projected_onto_xy_plane = acos(cos(cone_angle) / cos(inclination));
            double yaw_on_cone;
            if (to.y() > 0) {
                yaw_on_cone = M_PI/2 - angle_to_line_on_cone_projected_onto_xy_plane;
            } else {
                yaw_on_cone = -M_PI/2 + angle_to_line_on_cone_projected_onto_xy_plane;
            }
            yaw_change = atan2(from.y(), from.x()) - yaw_on_cone;
        }
        euler.yaw() += yaw_change;

        // Let from' = Rz(-yaw_change) * from
        // Solve: Ry(pitch_change) * to = from'
        from = Eigen::AngleAxisd(-yaw_change, Vector3d::UnitZ()) * from;
        double pitch_change = (-atan2(from.z(), from.x())) - (-atan2(to.z(), to.x()));
        euler.pitch() += pitch_change;

        // Avoid crossing -pi/2, pi/2 pitch boundary
        constexpr double pitch_limit = M_PI * 0.49;
        euler.pitch() = std::clamp(euler.pitch(), -pitch_limit, pitch_limit);

        pose.rotation() = euler.toRotationMatrix();
    }
    else if (state == State::SCROLL)
    {
        // Mapping of scroll_amount -> distance is:
        // distance = exp(scroll_amount)
        // d' = exp(s') = exp(ln(d) + delta_s)
        // where delta_s = scroll_amount * scroll_rate;

        float distance = click_pos_cs.norm();
        const float scroll_rate = -0.1;
        const float exp_dist_cutoff = 20;
        float new_distance;
        if (distance < exp_dist_cutoff) {
            new_distance = exp(log(distance) + scroll_amount * scroll_rate);
        } else {
            // Maintain gradient at exp_dist_cutoff = exp(x)
            new_distance = distance + scroll_amount * scroll_rate * exp_dist_cutoff;
        }

        pose = prev_camera_pose;
        pose.translation() -= pose.rotation() * click_pos_cs.normalized() * (new_distance - distance);
    }

    camera.update_view(pose);
}

bool CameraController::is_clicked(Vector3d& click_pos)
{
    if (state != State::INVALID && state != State::INACTIVE)
    {
        click_pos = click_pos_ws;
        return true;
    }
    return false;
}

} // namespace rbm
