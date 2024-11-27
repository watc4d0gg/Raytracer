//
// Created by Matej on 11/23/2024.
//
#pragma once

#include "geometry.h"
#include "image.h"
#include "ray.h"


class Camera {
public:
    Camera(
        const Vector2i& renderResolution,
        float fieldOfView,
        float exposure,
        const Point3f& position,
        const Point3f& lookAt,
        const Vector3f& up = Vector3f(0.f, 1.f, 0.f)
    );
    [[nodiscard]] Image& getRenderOutput();
    [[nodiscard]] float getFieldOfView() const;
    [[nodiscard]] float getAspectRatio() const;
    [[nodiscard]] float getExposure() const;
    [[nodiscard]] Point3f getPosition() const;
    [[nodiscard]] Vector3f forward() const;
    [[nodiscard]] Vector3f up() const;
    [[nodiscard]] Vector3f right() const;
    [[nodiscard]] Point3f getLookAt() const;
    void setLookAt(const Point3f& lookAt);
    [[nodiscard]] float getDistanceFromLookAt() const;
    void setDistanceFromLookAt(float distance);
    [[nodiscard]] Vector3f getEulerAngles() const;
    void setEulerAngles(const Vector3f& angles);
    [[nodiscard]] Ray castRay(const Vector2i& pixel) const;
private:
    Image renderOutput;

    float fieldOfView;
    float halfScreenHeight;
    float halfScreenWidth;

    float exposure;

    Point3f lookAt;
    float distanceFromLookAt;
    Vector3f eulerAngles;
};