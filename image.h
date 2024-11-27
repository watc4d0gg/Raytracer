//
// Created by Matej on 11/26/2024.
//
#pragma once
#include <vector>

#include "geometry.h"

constexpr float MAX_COLOR = 65535.f;


class Image {
public:
    explicit Image(const Vector2i& resolution);
    explicit Image(const std::string& filename);
    [[nodiscard]] Vector2i getResolution() const;
    void fill(const Vector3f& color);
    [[nodiscard]] Vector3f texel(const Point2f& uv) const;
    void splatPixel(int pixelX, int pixelY, const Vector3f& color);
    void writeToFile(const std::string& filename) const;
    [[nodiscard]] const std::vector<Vector3f>& getPixelData() const;
    std::vector<Vector3f>& getPixelData();
private:
    Vector2i resolution;
    std::vector<Vector3f> pixelData;
};
