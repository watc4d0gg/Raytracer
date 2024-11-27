//
// Created by Matej on 11/23/2024.
//

#include "camera.h"

#include "geometry.h"


Camera::Camera(
    const Vector2i& renderResolution,
    const float fieldOfView,
    const float exposure,
    const Point3f& position,
    const Point3f& lookAt,
    const Vector3f& up):
    renderOutput(Image(renderResolution)),
    fieldOfView(fieldOfView),
    halfScreenHeight(std::tan(radians(fieldOfView) / 2.f)),
    halfScreenWidth(getAspectRatio() * halfScreenHeight),
    exposure(exposure),
    lookAt(lookAt),
    distanceFromLookAt(lookAt.distance(position)),
    eulerAngles(Quaternion::lookAt(position, lookAt, up).toEulerAngles()) {
}

Image& Camera::getRenderOutput() {
    return renderOutput;
}

float Camera::getFieldOfView() const {
    return fieldOfView;
}

float Camera::getAspectRatio() const {
    const auto resolution = renderOutput.getResolution();

    if (resolution.x == 0 || resolution.y == 0) {
        return 1.f;
    }

    return static_cast<float>(resolution.x) / static_cast<float>(resolution.y);
}

float Camera::getExposure() const {
    return exposure;
}

Point3f Camera::getPosition() const {
    return lookAt - distanceFromLookAt * forward();
}

Vector3f Camera::forward() const {
    return Quaternion::fromEulerAngles(eulerAngles) * Vector3f(0.f, 0.f, 1.f);
}

Vector3f Camera::up() const {
    return Quaternion::fromEulerAngles(eulerAngles) * Vector3f(0.f, 1.f, 0.f);
}

Vector3f Camera::right() const {
    return Quaternion::fromEulerAngles(eulerAngles) * Vector3f(1.f, 0.f, 0.f);
}

Point3f Camera::getLookAt() const {
    return lookAt;
}

void Camera::setLookAt(const Point3f& lookAt) {
    if (this->lookAt != lookAt) {
        this->lookAt = lookAt;
    }
}

float Camera::getDistanceFromLookAt() const {
    return distanceFromLookAt;
}

void Camera::setDistanceFromLookAt(const float distance) {
    if (this->distanceFromLookAt != distance) {
        this->distanceFromLookAt = distance;
    }
}

Vector3f Camera::getEulerAngles() const {
    return eulerAngles;
}

void Camera::setEulerAngles(const Vector3f& angles) {
    if (this->eulerAngles != angles) {
        this->eulerAngles = angles;
    }
}

Ray Camera::castRay(const Vector2i& pixel) const {
    const auto resolution = renderOutput.getResolution();
    const Vector2f normalizedPixel = {
        static_cast<float>(pixel.x) / static_cast<float>(resolution.x) * 2.f - 1.f,
        static_cast<float>(pixel.y) / static_cast<float>(resolution.y) * 2.f - 1.f,
    };
    const Vector3f ndc = {
        normalizedPixel.x * halfScreenWidth,
        -normalizedPixel.y * halfScreenHeight,
        1.f,
    };
    return { getPosition(), Quaternion::fromEulerAngles(eulerAngles) * ndc.normalized() };
}



