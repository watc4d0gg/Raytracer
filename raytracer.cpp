#include "raytracer.h"

#include <array>
#include <iostream>
#include <memory>
#include <thread>

#include "camera.h"
#include "imgui.h"
#include "light.h"
#include "parser.h"
#include "scene.h"
#include "viewport.h"


void gammaCorrect(Vector3f& color) {
    color.x = std::pow(color.x, 2.2f);
    color.y = std::pow(color.y, 2.2f);
    color.z = std::pow(color.z, 2.2f);
}

// https://frictionalgames.blogspot.com/2012/09/tech-feature-hdr-lightning.html
Vector3f toneFilmic(const Vector3f& color) {
    constexpr auto A = .15f;
    constexpr auto B = .50f;
    constexpr auto C = .10f;
    constexpr auto D = .20f;
    constexpr auto E = .02f;
    constexpr auto F = .30f;

    return (color * (A * color + C * B) + D * E) / (color * (A * color + B) + D * F) - E / F;
}

Vector3f blinnPhong(Ray& ray, Scene& scene, const int bounces = 5) {
    auto color = scene.ambientColor;

    if (const auto& interaction = scene.aggregate.intersect(ray)) {
        const auto& material = interaction->material;
        color *= 0.1f;

        // Light contribution
        {
            for (const auto& light : scene.lights) {
                auto lightDirection = light.position() - interaction->intersection;
                auto distance = lightDirection.length();
                auto visibilityRay = interaction->spawnRay(lightDirection);

                if (const auto& possibleOcclusion = scene.aggregate.intersect(visibilityRay)) {
                    if (possibleOcclusion->tHit < distance) {
                        continue;
                    }
                }

                // Diffuse light contribution
                const auto diffuseIntensity = std::max(0.f, interaction->normal.dot(lightDirection.normalized()));
                color += light.intensity() * material.kd * material.diffuseColor(interaction->uv) * diffuseIntensity;

                // Specular light contribution
                if (diffuseIntensity > 0.f) {
                    const auto halfVector = (lightDirection - ray.direction).normalized();
                    const auto specularIntensity = std::max(0.f, interaction->normal.dot(halfVector));
                    color += light.intensity() * material.ks * material.specularColor(interaction->uv) * std::pow(specularIntensity, material.alpha);
                }
            }
        }

        if (bounces > 0) {

            // Reflective contribution
            {
                if (material.isReflective) {
                    color *= 1.f - material.reflectivity;
                    auto reflectedRay = interaction->spawnRay(ray.direction - 2.f * interaction->normal.dot(ray.direction) * interaction->normal);
                    color += material.reflectivity * blinnPhong(reflectedRay, scene, bounces - 1);
                }
            }

            // Refractive contribution
            // Fundamentals of Computer Graphics, Fourth Edition
            {
                if (material.isRefractive) {
                    // Make the rest of the material transparent
                    color *= 0.f;

                    // The ray will reflect if there is total internal reflection, bypassing the isReflective flag
                    auto reflectedRay = interaction->spawnRay(ray.direction - 2.f * interaction->normal.dot(ray.direction) * interaction->normal);

                    const auto refractiveIndex = interaction->normal.dot(ray.direction) < 0 ? 1.f / material.refractiveIndex : material.refractiveIndex;
                    auto incomingAngle = -interaction->normal.dot(ray.direction);
                    const auto refractedAngleSquared = 1.f - std::pow(refractiveIndex, 2.f) * (1.f - std::pow(interaction->normal.dot(ray.direction), 2.f));
                    Vector3f attenuation;

                    if (incomingAngle < 0.f) {
                        attenuation = material.specularColor(interaction->uv) * std::exp(-interaction->tHit);
                        // Total internal reflection
                        if (refractedAngleSquared < 0.f) {
                            color += attenuation * blinnPhong(reflectedRay, scene, bounces - 1);
                            return color;
                        }
                    } else {
                        attenuation = Vector3f(1.f);
                    }

                    const auto refractedPerpendicular = refractiveIndex * (ray.direction + incomingAngle * interaction->normal);
                    const auto refractedParallel = std::sqrt(refractedAngleSquared) * interaction->normal;

                    // Schlick's approximation
                    const auto reflectance = std::pow((refractiveIndex - 1) / (refractiveIndex + 1), 2.f);
                    const auto reflectivity = reflectance + (1 - reflectance) * std::pow(1 - incomingAngle, 5.f);

                    auto refractedRay = interaction->spawnRay(refractedPerpendicular - refractedParallel);
                    color += attenuation * (reflectivity * blinnPhong(reflectedRay, scene, bounces - 1) + (1.f - reflectivity) * blinnPhong(refractedRay, scene, bounces - 1));
                }
            }
        }
    }

    return color;
}

int main() {
    std::cout << "Hello world!" << std::endl;

    // Parse the raytracer's scene
    int maxBounces;
    RenderMode renderMode;
    CameraData cameraData {};
    SceneData sceneData {};

    loadSceneFromJson("../scenes/mirror_image_refract.json", maxBounces, renderMode, cameraData, sceneData);

    Camera camera {
        cameraData.renderResolution,
        cameraData.fieldOfView,
        cameraData.exposure,
        cameraData.position,
        cameraData.lookAt,
        cameraData.up
    };

    Scene scene {
        sceneData.ambientColor,
        sceneData.lights,
        sceneData.primitives,
    };

    Image& output = camera.getRenderOutput();

    bool keepRendering = true;

    // TODO: allow CLI rendering
    Viewport* viewport;

    viewport = new Viewport { "Viewport Render", { 800, 600 }, camera };


    // Initialize multithreaded rendering
    std::vector<std::thread> renderThreads;
    {
        // Setup render threads per image tile
        const auto nThreads = static_cast<int>(std::thread::hardware_concurrency());
        const auto hTiles = static_cast<int>(std::floor(std::sqrt(static_cast<float>(nThreads) / camera.getAspectRatio())));
        const auto wTiles = nThreads / hTiles;
        const auto& resolution = camera.getRenderOutput().getResolution();
        renderThreads.reserve(wTiles * hTiles);

        // Compute tile sizes
        const auto tileWidth = resolution.x / wTiles;
        const auto tileHeight = resolution.y / hTiles;
        const auto lastWidth = tileWidth + (resolution.x % wTiles);
        const auto lastHeight = tileHeight + (resolution.y % hTiles);

        for (int i = 0; i < wTiles; i++) {
            for (int j = 0; j < hTiles; j++) {
                const auto startX = tileWidth * i;
                const auto startY = tileHeight * j;
                const auto width = i == wTiles - 1 ? lastWidth : tileWidth;
                const auto height = j == hTiles - 1 ? lastHeight : tileHeight;
                renderThreads.emplace_back([&, startX, startY, width, height] {
                    do {
                        // std::cout << "Rendering..." << std::endl;
                        for (int x = startX; x < startX + width; x++) {
                            for (int y = startY; y < startY + height; y++) {
                                Ray ray = camera.castRay({ x, y });
                                // TODO: switch between render modes
                                auto color = blinnPhong(ray, scene, maxBounces);
                                gammaCorrect(color);
                                color *= 100.f * camera.getExposure();
                                const auto toneMapped = toneFilmic(color) / toneFilmic(Vector3f(std::pow(100.f * camera.getExposure(), 2.2f)));
                                output.splatPixel(x, y, toneMapped);
                            }
                        }
                        // std::cout << "Frame rendered!" << std::endl;
                    } while (keepRendering);
                });
            }
        }
    }

    // Update GUI
    while (!viewport->shouldClose()) {
        viewport->updateInput();
        ImGui::Begin("Rendering");
        ImGui::Text("Rendering");
        ImGui::End();
        // for (int x = 0; x < output.getResolution().x; x++) {
        //     for (int y = 0; y < output.getResolution().y; y++) {
        //         Ray ray = camera.castRay({ x, y });
        //         // TODO: switch between render modes
        //         auto color = blinnPhong(ray, scene, maxBounces);
        //         gammaCorrect(color);
        //         color *= 100.f * camera.getExposure();
        //         const auto toneMapped = toneMapColor(color) / toneMapColor(Vector3f(std::pow(100.f * camera.getExposure(), 2.2f)));
        //         output.splatPixel(x, y, toneMapped);
        //     }
        // }
        viewport->refresh();
    }

    delete viewport;
    keepRendering = false;

    // Wait for the render threads to finish after the GUI was closed
    for (auto& thread : renderThreads) {
        thread.join();
    }

    output.writeToFile("../output.ppm");

    delete viewport;
    return 0;
}