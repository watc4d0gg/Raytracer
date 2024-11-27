//
// Created by Matej on 11/25/2024.
//
#pragma once

#include <memory>
#include <vector>

#include "geometry.h"
#include "light.h"
#include "primitive.h"
#include "raytracer.h"

#include <nlohmann/json.hpp>

using json = nlohmann::json;


struct CameraData {
    Vector2i renderResolution;
    float fieldOfView;
    float exposure;
    Point3f position;
    Point3f lookAt;
    Vector3f up;
};

struct SceneData {
    Vector3f ambientColor;
    std::vector<PointLight> lights;
    std::vector<std::shared_ptr<Primitive>> primitives;
};

struct MeshData {
    std::vector<Point3f> vertices;
    std::vector<Normal3f> normals;
    std::vector<Point2f> uvs;
};

void loadSceneFromJson(const std::string& filename, int& maxBounces, RenderMode& renderMode, CameraData& cameraData, SceneData& sceneData);
