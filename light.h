//
// Created by Matej on 11/25/2024.
//
#pragma once
#include "geometry.h"

class Light {
public:
    virtual ~Light() = default;

    [[nodiscard]] virtual const Point3f& position() const;
    [[nodiscard]] virtual const Vector3f& intensity() const;
    // TODO: sample methods
};

class PointLight final : public Light {
public:
    PointLight(const Point3f& position, const Vector3f& intensity) : m_position(position), m_intensity(intensity) {}

    [[nodiscard]] const Point3f& position() const override {
        return m_position;
    }
    [[nodiscard]] const Vector3f& intensity() const override {
        return m_intensity;
    }

    Point3f m_position;
    Vector3f m_intensity;
};
