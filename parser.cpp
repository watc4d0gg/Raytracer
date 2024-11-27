//
// Created by Matej on 11/25/2024.
//

#include "parser.h"

#include <fstream>
#include <memory>
#include "material.h"

Material parseMaterial(json shape) {
    if (shape["material"].is_null()) {
        return {};
    }
    auto material = Material {};
    material.ks = shape["material"]["ks"].get<float>();
    material.kd = shape["material"]["kd"].get<float>();
    material.alpha = shape["material"]["specularexponent"].get<float>();

    std::variant<Vector3f, std::shared_ptr<Image>> diffuse;
    if (shape["material"]["diffusecolor"].is_string()) {
        diffuse = std::make_shared<Image>(shape["material"]["diffusecolor"].get<std::string>());
    } else {
        Vector3f color;
        color.x = shape["material"]["diffusecolor"][0].get<float>();
        color.y = shape["material"]["diffusecolor"][1].get<float>();
        color.z = shape["material"]["diffusecolor"][2].get<float>();
        diffuse = color;
    }
    material.diffuse = diffuse;

    std::variant<Vector3f, std::shared_ptr<Image>> specular;
    if (shape["material"]["specularcolor"].is_string()) {
        specular = std::make_shared<Image>(shape["material"]["specularcolor"].get<std::string>());
    } else {
        Vector3f color;
        color.x = shape["material"]["specularcolor"][0].get<float>();
        color.y = shape["material"]["specularcolor"][1].get<float>();
        color.z = shape["material"]["specularcolor"][2].get<float>();
        specular = color;
    }
    material.specular = specular;

    // material.diffuseColor.x = shape["material"]["diffusecolor"][0].get<float>();
    // material.diffuseColor.y = shape["material"]["diffusecolor"][1].get<float>();
    // material.diffuseColor.z = shape["material"]["diffusecolor"][2].get<float>();
    // material.specularColor.x = shape["material"]["specularcolor"][0].get<float>();
    // material.specularColor.y = shape["material"]["specularcolor"][1].get<float>();
    // material.specularColor.z = shape["material"]["specularcolor"][2].get<float>();

    material.isReflective = shape["material"]["isreflective"].get<bool>();
    material.reflectivity = shape["material"]["reflectivity"].get<float>();
    material.isRefractive = shape["material"]["isrefractive"].get<bool>();
    material.refractiveIndex = shape["material"]["refractiveindex"].get<float>();
    return material;
}

void loadSceneFromJson(const std::string& filename, int& maxBounces, RenderMode& renderMode, CameraData& cameraData, SceneData& sceneData) {
    std::ifstream sceneJson(filename);
    json data = json::parse(sceneJson);

    if (!data["nbounces"].is_null()) {
        maxBounces = data["nbounces"].get<int>();
    }

    const auto mode = data["rendermode"].get<std::string>();
    if (mode == "binary") {
        renderMode = Binary;
    } else if (mode == "phong") {
        renderMode = Phong;
    } else if (mode == "pathtracing") {
        renderMode = Pathtracing;
    }

    // Initialize camera data
    cameraData.renderResolution.x = data["camera"]["width"].get<int>();
    cameraData.renderResolution.y = data["camera"]["height"].get<int>();
    cameraData.position.x = data["camera"]["position"][0].get<float>();
    cameraData.position.y = data["camera"]["position"][1].get<float>();
    cameraData.position.z = data["camera"]["position"][2].get<float>();
    cameraData.lookAt.x = data["camera"]["lookAt"][0].get<float>();
    cameraData.lookAt.y = data["camera"]["lookAt"][1].get<float>();
    cameraData.lookAt.z = data["camera"]["lookAt"][2].get<float>();
    cameraData.up.x = data["camera"]["upVector"][0].get<float>();
    cameraData.up.y = data["camera"]["upVector"][1].get<float>();
    cameraData.up.z = data["camera"]["upVector"][2].get<float>();
    cameraData.fieldOfView = data["camera"]["fov"].get<float>();
    cameraData.exposure = data["camera"]["exposure"].get<float>();

    // Initialize scene data
    const auto scene = data["scene"];

    sceneData.ambientColor.x = scene["backgroundcolor"][0].get<float>();
    sceneData.ambientColor.y = scene["backgroundcolor"][1].get<float>();
    sceneData.ambientColor.z = scene["backgroundcolor"][2].get<float>();

    for (const auto& light : scene["lightsources"]) {
        const auto type = light["type"].get<std::string>();
        if (type == "pointlight") {
            auto pointLight = PointLight {
                { light["position"][0].get<float>(), light["position"][1].get<float>(), light["position"][2].get<float>() },
                { light["intensity"][0].get<float>(), light["intensity"][1].get<float>(), light["intensity"][2].get<float>() },
            };
            sceneData.lights.emplace_back(pointLight);
        }
    }

    for (const auto& shape : scene["shapes"]) {
        const auto type = shape["type"].get<std::string>();
        if (type == "sphere") {
            auto sphere = new Sphere {
                { shape["center"][0].get<float>(), shape["center"][1].get<float>(), shape["center"][2].get<float>() },
                shape["radius"].get<float>(),
                parseMaterial(shape),
            };

            sceneData.primitives.emplace_back(sphere);
        } else if (type == "cylinder") {
            auto cylinder = new Cylinder {
                { shape["center"][0].get<float>(), shape["center"][1].get<float>(), shape["center"][2].get<float>() },
                { shape["axis"][0].get<float>(), shape["axis"][1].get<float>(), shape["axis"][2].get<float>() },
                shape["radius"].get<float>(),
                shape["height"].get<float>(),
                parseMaterial(shape),
            };
            sceneData.primitives.emplace_back(cylinder);
        } else if (type == "triangle") {
            auto triangle = new Triangle {
                    { shape["v0"][0].get<float>(), shape["v0"][1].get<float>(), shape["v0"][2].get<float>() },
                    { shape["v1"][0].get<float>(), shape["v1"][1].get<float>(), shape["v1"][2].get<float>() },
                    { shape["v2"][0].get<float>(), shape["v2"][1].get<float>(), shape["v2"][2].get<float>() },
                    parseMaterial(shape),
                };
            sceneData.primitives.emplace_back(triangle);
        }
    }

    sceneJson.close();
}
