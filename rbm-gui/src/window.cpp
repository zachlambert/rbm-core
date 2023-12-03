#include "sviz/gui/window.h"
#include "sviz/gui/font.h"
#include "imgui.h"
#include "backends/imgui_impl_glfw.h"
#include "backends/imgui_impl_opengl3.h"

#include <iostream>
#include <unordered_map>
#include <filesystem>


namespace sviz {

bool Window::init()
{
    auto font_path = find_font_path(config.font);
    if (!font_path.has_value()) {
        std::cerr << "Failed to find font" << std::endl;
        return false;
    }

    if (!glfwInit()) {
        std::cerr << "Failed to initialise glfw" << std::endl;
        return false;
    }

    // GL 3.0 + GLSL 130
    const char* glsl_version = "#version 130";
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 0);

    // Create window with graphics context
    window = glfwCreateWindow(config.width, config.height, title.c_str(), nullptr, nullptr);
    if (window == nullptr) {
        std::cerr << "Failed to create window" << std::endl;
        return false;
    }

    glfwMakeContextCurrent(window);

    if (config.vsync) {
        glfwSwapInterval(1);
    } else {
        glfwSwapInterval(0);
    }

    if (glewInit() != GLEW_OK) {
        std::cerr << "Failed to initialise glew" << std::endl;
        return false;
    }

    // Setup Dear ImGui context
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO &io = ImGui::GetIO();

    // Setup Dear ImGui style
    ImGui::StyleColorsDark();

    // Setup Platform/Renderer backends
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init(glsl_version);

    io.Fonts->AddFontFromFileTTF(font_path.value().c_str(), config.font_size);

    io.KeyRepeatDelay = 0.12; // Default 0.275
    io.KeyRepeatRate = 0.005; // Default 0.05
    // io.IniFilename = nullptr; // Disable saving ini file
    io.ConfigFlags |= ImGuiConfigFlags_DockingEnable;
    io.ConfigWindowsMoveFromTitleBarOnly = true;

    for (auto& [name, widget]: widgets) {
        if (!widget->init()) {
            std::cerr << "Failed to initialise widget '" << name << "'" << std::endl;
            return false;
        }
    }

    initialised = true;

    return true;
}

void Window::shutdown()
{
    if (!initialised) return;
    initialised = false;

    for (auto& [name, widget]: widgets) {
        widget->shutdown();
    }

    glfwMakeContextCurrent(window);

    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();

    glfwDestroyWindow(window);
    glfwTerminate();
}

void Window::poll_events()
{
    glfwPollEvents();

    GLFWgamepadstate gamepad_state;
    user_input.valid = glfwGetGamepadState(GLFW_JOYSTICK_1, &gamepad_state);
    if (user_input.valid) {
        auto threshold = [](double axis) -> double {
            const double min = 0.2;
            if (axis < -min) return (axis + min) / (1.0 - min);
            if (axis > min) return (axis - min) / (1.0 - min);
            return 0;
        };
        user_input.axes[(size_t)UserInput::Axis::LEFT_X] = -threshold(gamepad_state.axes[GLFW_GAMEPAD_AXIS_LEFT_X]);
        user_input.axes[(size_t)UserInput::Axis::LEFT_Y] = -threshold(gamepad_state.axes[GLFW_GAMEPAD_AXIS_LEFT_Y]);
        user_input.axes[(size_t)UserInput::Axis::RIGHT_X] = -threshold(gamepad_state.axes[GLFW_GAMEPAD_AXIS_RIGHT_X]);
        user_input.axes[(size_t)UserInput::Axis::RIGHT_Y] = -threshold(gamepad_state.axes[GLFW_GAMEPAD_AXIS_RIGHT_Y]);

        user_input.buttons[(size_t)UserInput::Button::A] = gamepad_state.buttons[GLFW_GAMEPAD_BUTTON_A] == GLFW_PRESS;
        user_input.buttons[(size_t)UserInput::Button::B] = gamepad_state.buttons[GLFW_GAMEPAD_BUTTON_B] == GLFW_PRESS;
        user_input.buttons[(size_t)UserInput::Button::X] = gamepad_state.buttons[GLFW_GAMEPAD_BUTTON_X] == GLFW_PRESS;
        user_input.buttons[(size_t)UserInput::Button::Y] = gamepad_state.buttons[GLFW_GAMEPAD_BUTTON_Y] == GLFW_PRESS;
        user_input.buttons[(size_t)UserInput::Button::LB] = gamepad_state.buttons[GLFW_GAMEPAD_BUTTON_LEFT_BUMPER] == GLFW_PRESS;
        user_input.buttons[(size_t)UserInput::Button::RB] = gamepad_state.buttons[GLFW_GAMEPAD_BUTTON_RIGHT_BUMPER] == GLFW_PRESS;
    }
}

void Window::render_start()
{
    glClearColor(0, 0, 0, 0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // Start the ImGui frame
    ImGui_ImplOpenGL3_NewFrame();
    ImGui_ImplGlfw_NewFrame();
    ImGui::NewFrame();
    ImGui::DockSpaceOverViewport();
}

void Window::render_end()
{
    ImGui::Render();
    int display_w, display_h;
    glfwGetFramebufferSize(window, &display_w, &display_h);
    glViewport(0, 0, display_w, display_h);
    ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

    config.width = display_w;
    config.height = display_h;
    glfwSwapBuffers(window);
}

} // namespace sviz
