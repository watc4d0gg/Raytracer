//
// Created by Matej on 11/23/2024.
//
#pragma once

#include <iostream>
#include <optional>

#include "geometry.h"
#include "material.h"

class Ray {
public:
    Ray(const Point3f& origin, const Vector3f& direction, const float time = std::numeric_limits<float>::max()):
        origin(origin), direction(direction.normalized()), time(time) {};

    Point3f operator()(const float time) const {
        return origin + direction * time;
    }

    // https://en.wikipedia.org/wiki/Slab_method
    [[nodiscard]] std::optional<float> intersects(const Bounds& bounds) const {
        const auto tMin = (bounds.minPoint - origin) / direction;
        const auto tMax = (bounds.maxPoint - origin) / direction;
        const auto t1 = min(tMin, tMax);
        const auto t2 = max(tMin, tMax);
        const auto tNear = std::max(std::max(t1.x, t1.y), t1.z);
        const auto tFar = std::min(std::min(t2.x, t2.y), t2.z);
        return tNear <= tFar ? std::optional(tNear) : std::nullopt;
    }

    Point3f origin;
    Vector3f direction;
    float time;
};

class RayInteraction {
public:
    [[nodiscard]] Ray spawnRay(const Vector3f& direction) const {
        Ray newRay = { intersection, direction };
        // Offset slightly to not cause interactions with the point it is spawned from
        newRay.origin += .0001f * newRay.direction;
        return newRay;
    }

    Point3f intersection {};
    Normal3f normal {};
    Point2f uv {};
    float tHit {};
    Material material {};
};

