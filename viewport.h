//
// Created by Matej on 11/20/2024.
//
#pragma once
#include <string>

#include "camera.h"
#include "geometry.h"
#include "GLFW/glfw3.h"

// A convenient wrapper around GLFW window and ImGui
class Viewport {
public:
    Viewport(const std::string& title, const Vector2i& size, Camera& camera);
    ~Viewport();

    Vector2i getSize();
    [[nodiscard]] float getDpiScale() const;
    void close() const;
    [[nodiscard]] bool shouldClose() const;
    void updateInput();
    void refresh() const;
    [[nodiscard]] bool isKeyPressed(int key) const;
    [[nodiscard]] bool isMouseButtonPressed(int button) const;
    [[nodiscard]] Vector2f getMousePosition() const;

    // using KeyCallback = std::function<void(int key, int scancode, int action, int mods)>;
    // using MouseButtonCallback = std::function<void(int button, int action, int mods)>;
    // using MouseMoveCallback = std::function<void(double posX, double posY)>;
    // using ScrollCallback = std::function<void(double offsetX, double offsetY)>;

private:
    static void keyCallback(GLFWwindow* window, int key, int scancode, int action, int mods);
    static void mouseButtonCallback(GLFWwindow* window, int button, int action, int mods);
    static void mouseMoveCallback(GLFWwindow* window, double posX, double posY);
    static void scrollCallback(GLFWwindow* window, double offsetX, double offsetY);

    GLFWwindow* window;
    Camera& camera;
    GLuint renderedImage;
    Vector2i size;
    float dpiScale = 1.f;
    Vector2f prevMousePosition {};
};