//
// Created by Matej on 11/25/2024.
//
#pragma once
#include <memory>
#include <variant>

#include "geometry.h"
#include "image.h"


class Material {
public:
    Material() = default;
    virtual ~Material() = default;

    [[nodiscard]] Vector3f diffuseColor(const Point2f& uv = Point2f { 0.f }) const {
        if (std::holds_alternative<std::shared_ptr<Image>>(diffuse)) {
            return kd * std::get<std::shared_ptr<Image>>(diffuse)->texel(uv);
        }
        return std::get<Vector3f>(diffuse);
    }

    [[nodiscard]] Vector3f specularColor(const Point2f& uv = Point2f { 0.f }) const {
        if (std::holds_alternative<std::shared_ptr<Image>>(specular)) {
            return std::get<std::shared_ptr<Image>>(specular)->texel(uv);
        }
        return std::get<Vector3f>(specular);
    }

    float kd;
    float ks;
    float ke;
    float alpha;
    std::variant<Vector3f, std::shared_ptr<Image>> diffuse;
    std::variant<Vector3f, std::shared_ptr<Image>> specular;
    // Vector3f diffuseColor;
    // Vector3f specularColor;
    bool isReflective;
    float reflectivity;
    bool isRefractive;
    float refractiveIndex;
};
