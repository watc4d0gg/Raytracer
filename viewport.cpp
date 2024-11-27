//
// Created by Matej on 11/20/2024.
//

#include "viewport.h"

#include <algorithm>
#include <iostream>
#include <numeric>

#include "geometry.h"
#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"
#include "imgui_internal.h"

static constexpr float rotationSpeedFactor = .3f;
static constexpr float translationSpeedFactor = .005f;
static constexpr float zoomSpeedFactor = .5f;

static void errorCallback(const int error, const char* description) {
    std::cerr << "GLFW error: " << error << std::endl;
    std::cerr << description << std::endl;
    exit(1);
}

static void framebufferSizeCallback(GLFWwindow* window, const int width, const int height) {
    glViewport(0, 0, width, height);
}

Viewport::Viewport(const std::string& title, const Vector2i& size, Camera& camera): camera(camera) {
    glfwSetErrorCallback(errorCallback);
    if (!glfwInit()) {
        std::cerr << "Failed to initialize GLFW" << std::endl;
        exit(1);
    }

    glfwWindowHint(GLFW_VISIBLE, GLFW_TRUE);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 5);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_COMPAT_PROFILE);

    GLFWmonitor* monitor = glfwGetPrimaryMonitor();
    if (!monitor) {
        std::cerr << "Failed to get GLFW monitor" << std::endl;
        exit(1);
    }

    float scaleX, scaleY;
    glfwGetMonitorContentScale(monitor, &scaleX, &scaleY);
    if (scaleX > 1.f || scaleY > 1.f) {
        dpiScale = scaleX;
        glfwWindowHint(GLFW_SCALE_TO_MONITOR, GLFW_TRUE);
    }

    window = glfwCreateWindow(size.x, size.y, title.c_str(), nullptr, nullptr);
    if (window == nullptr) {
        std::cerr << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        exit(1);
    }

    glfwMakeContextCurrent(window);
    glfwSwapInterval(1); // Enable v-sync

    // Keep the aspect ratio the same as the render
    const Vector2i resolution = camera.getRenderOutput().getResolution();
    const int divisor = std::gcd(resolution.x, resolution.y);
    glfwSetWindowAspectRatio(window, resolution.x / divisor, resolution.y / divisor);
    glfwGetWindowSize(window, &this->size.x, &this->size.y);
    glfwSetWindowUserPointer(window, this);

    glfwSetKeyCallback(window, keyCallback);
    glfwSetMouseButtonCallback(window, mouseButtonCallback);
    glfwSetCursorPosCallback(window, mouseMoveCallback);
    glfwSetScrollCallback(window, scrollCallback);
    glfwSetFramebufferSizeCallback(window, framebufferSizeCallback);

    // Create a texture for rendered output to view in the viewport
    glGenTextures(1, &renderedImage);
    glBindTexture(GL_TEXTURE_2D, renderedImage);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glBindTexture(GL_TEXTURE_2D, 0);

    // Initialize ImGui
    ImGui::CreateContext();
    const ImGuiIO &io = ImGui::GetIO(); (void)io;
    ImGui::StyleColorsDark();
    ImGuiStyle &style = ImGui::GetStyle();
    style.ScaleAllSizes(dpiScale);
    ImFontConfig config;
    config.SizePixels = 13 * dpiScale;
    io.Fonts->AddFontDefault(&config);

    if (!ImGui_ImplGlfw_InitForOpenGL(window, true)) {
        std::cerr << "Failed to initialize imgui" << std::endl;
        exit(1);
    }

    if (!ImGui_ImplOpenGL3_Init()) {
        std::cerr << "Failed to initialize imgui" << std::endl;
        exit(1);
    }
}

Viewport::~Viewport() {
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();
    glfwDestroyWindow(window);
    glfwTerminate();
}

Vector2i Viewport::getSize() {
    return size;
}

float Viewport::getDpiScale() const {
    return dpiScale;
}

void Viewport::close() const {
    glfwSetWindowShouldClose(window, true);
}

bool Viewport::shouldClose() const {
    return glfwWindowShouldClose(window);
}

void Viewport::updateInput() {
    glfwPollEvents();
    ImGui_ImplOpenGL3_NewFrame();
    ImGui_ImplGlfw_NewFrame();
    ImGui::NewFrame();
}

void Viewport::refresh() const {
    // Render the texture
    const auto renderOutput = camera.getRenderOutput();
    const auto& resolution = renderOutput.getResolution();

    glPushAttrib(GL_ALL_ATTRIB_BITS);

    glBindTexture(GL_TEXTURE_2D, renderedImage);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, resolution.x, resolution.y, 0, GL_RGB, GL_FLOAT, renderOutput.getPixelData().data());

    glDisable(GL_LIGHTING);
    glDisable(GL_LIGHT0);
    glDisable(GL_COLOR_MATERIAL);
    glDisable(GL_NORMALIZE);
    glColor3f(1.0f, 1.0f, 1.0f);

    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, renderedImage);

    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();

    glBegin(GL_QUADS);
    glTexCoord2f(0.f, 1.f);
    glVertex3f(-1.f, -1.f, 0.f);
    glTexCoord2f(1.f, 1.f);
    glVertex3f(1.f, -1.f, 0.f);
    glTexCoord2f(1.f, 0.f);
    glVertex3f(1.f, 1.f, 0.f);
    glTexCoord2f(0.f, 0.f);
    glVertex3f(-1.f, 1.f, 0.f);
    glEnd();

    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();

    glPopAttrib();

    ImGui::Render();
    ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
    glfwSwapBuffers(window);
}

bool Viewport::isKeyPressed(const int key) const {
    return glfwGetKey(window, key) == GLFW_PRESS;
}

bool Viewport::isMouseButtonPressed(const int button) const {
    return glfwGetMouseButton(window, button) == GLFW_PRESS;
}

Vector2f Viewport::getMousePosition() const {
    double x, y;
    glfwGetCursorPos(window, &x, &y);
    return { static_cast<float>(x), static_cast<float>(y) };
}

void Viewport::keyCallback(GLFWwindow* window, const int key, const int scancode, const int action, const int mods) {
    ImGui_ImplGlfw_KeyCallback(window, key, scancode, action, mods);

    if (ImGui::GetIO().WantCaptureKeyboard) {
        return;
    }

    const auto* viewport = static_cast<Viewport*>(glfwGetWindowUserPointer(window));

    // Register keys
    if (action == GLFW_PRESS) {
        switch (key) {
            // TODO: add more?
            case GLFW_KEY_ESCAPE: {
                viewport->close();
                std::cout << "here" << std::endl;
            } break;
        }
    }
}

void Viewport::mouseButtonCallback(GLFWwindow* window, const int button, const int action, const int mods) {

    if (ImGui::GetIO().WantCaptureMouse) {
        return;
    }

    auto* viewport = static_cast<Viewport*>(glfwGetWindowUserPointer(window));

    // Register rotation and translation clicks
    if ((button == GLFW_MOUSE_BUTTON_LEFT || button == GLFW_MOUSE_BUTTON_RIGHT) && action == GLFW_PRESS) {
        viewport->prevMousePosition = viewport->getMousePosition();
    }
}

void Viewport::mouseMoveCallback(GLFWwindow* window, const double posX, const double posY) {

    if (ImGui::GetIO().WantCaptureMouse) {
        return;
    }

    auto* viewport = static_cast<Viewport*>(glfwGetWindowUserPointer(window));

    // Register the computation of rotation and translation from mouse movements
    const Vector2f currentMousePosition = { static_cast<float>(posX), static_cast<float>(posY) };
    const auto rotate = viewport->isMouseButtonPressed(GLFW_MOUSE_BUTTON_LEFT);
    const auto translate = viewport->isMouseButtonPressed(GLFW_MOUSE_BUTTON_RIGHT);
    if (rotate || translate) {
        const auto diff = currentMousePosition - viewport->prevMousePosition;

        if (rotate) {
            const auto prevAngles = viewport->camera.getEulerAngles();
            const Vector3f angles = {
                std::clamp(prevAngles.x + radians(diff.y * rotationSpeedFactor), -pi / 2.f, pi / 2.f),
                prevAngles.y + radians(diff.x * rotationSpeedFactor),
                prevAngles.z,
            };
            viewport->camera.setEulerAngles(angles);

        } else {
            const auto lookAt = viewport->camera.getLookAt()
                - diff.x * translationSpeedFactor * viewport->camera.right()
                + diff.y * translationSpeedFactor * viewport->camera.up();
            viewport->camera.setLookAt(lookAt);
        }

        viewport->prevMousePosition = currentMousePosition;
    }
}

void Viewport::scrollCallback(GLFWwindow* window, const double offsetX, const double offsetY) {

    if (ImGui::GetIO().WantCaptureMouse) {
        return;
    }

    const auto* viewport = static_cast<Viewport*>(glfwGetWindowUserPointer(window));

    // Register distance changes from mouse scroll
    const auto distance = viewport->camera.getDistanceFromLookAt() - static_cast<float>(offsetY) * zoomSpeedFactor;
    viewport->camera.setDistanceFromLookAt(std::clamp(distance, 0.1f, 100.0f));
}







