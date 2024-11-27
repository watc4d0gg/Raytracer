//
// Created by Matej on 11/25/2024.
//
#pragma once
#include <vector>

#include "light.h"
#include "primitive.h"


class Scene {
public:
    Scene(const Vector3f& ambientColor, const std::vector<PointLight>& lights, std::vector<std::shared_ptr<Primitive>>& primitives):
        ambientColor(ambientColor),
        primitives(primitives),
        aggregate(BoundingVolumeHierarchy(primitives)),
        lights(lights) {}

    Vector3f ambientColor;
    std::vector<std::shared_ptr<Primitive>>& primitives;
    BoundingVolumeHierarchy aggregate;
    std::vector<PointLight> lights;
};
