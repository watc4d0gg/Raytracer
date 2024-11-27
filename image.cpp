//
// Created by Matej on 11/26/2024.
//

#include "image.h"

#include <fstream>
#include <iostream>
#include <bits/algorithmfwd.h>


Image::Image(const Vector2i& resolution):
    resolution(resolution),
    pixelData(static_cast<size_t>(resolution.x * resolution.y), Vector3f()) {}

Image::Image(const std::string& filename) {
    std::ifstream imageFile(filename);

    if (!imageFile.is_open()) {
        std::cerr << "Could not open file " << filename << std::endl;
        exit(1);
    }

    std::string magicNumber;
    int width, height, maxColor;

    imageFile >> magicNumber >> width >> height >> maxColor;

    if (magicNumber != "P3") {
        std::cerr << "Not a PPM (P3) image!" << std::endl;
        exit(1);
    }

    resolution = Vector2i(width, height);
    pixelData = std::vector(static_cast<size_t>(resolution.x * resolution.y), Vector3f());

    for (int y = 0; y < resolution.y; y++) {
        for (int x = 0; x < resolution.x; x++) {
            int r, g, b;
            imageFile >> r >> g >> b;
            splatPixel(x, y, Vector3f(static_cast<float>(r), static_cast<float>(g), static_cast<float>(b)) / static_cast<float>(maxColor));
        }
    }

    imageFile.close();
}

Vector2i Image::getResolution() const {
    return resolution;
}

void Image::fill(const Vector3f& color) {
    std::fill(std::begin(pixelData), std::end(pixelData), color);
}

Vector3f Image::texel(const Point2f& uv) const {
    const auto pixelX = static_cast<float>(resolution.x - 1) * uv.x;
    const auto pixelY = static_cast<float>(resolution.y - 1) * uv.y;

    // Interpolate on rounding error
    const Vector2i p1 {
        static_cast<int>(std::min(pixelX, resolution.x - 1.f)),
        static_cast<int>(std::min(pixelY, resolution.y - 1.f)),
    };
    const Vector2i p2 {
        static_cast<int>(std::min(pixelX + 1, resolution.x - 1.f)),
        static_cast<int>(std::min(pixelY, resolution.y - 1.f)),
    };
    const Vector2i p3 {
        static_cast<int>(std::min(pixelX, resolution.x - 1.f)),
        static_cast<int>(std::min(pixelY + 1, resolution.y - 1.f)),
    };
    const Vector2i p4 {
        static_cast<int>(std::min(pixelX + 1, resolution.x - 1.f)),
        static_cast<int>(std::min(pixelY + 1, resolution.y - 1.f)),
    };

    const auto c1 = pixelData[p1.y * resolution.x + p1.x];
    const auto c2 = pixelData[p2.y * resolution.x + p2.x];
    const auto c3 = pixelData[p3.y * resolution.x + p3.x];
    const auto c4 = pixelData[p4.y * resolution.x + p4.x];

    const auto alpha = pixelX - static_cast<float>(p1.x);
    const auto beta = pixelY - static_cast<float>(p1.y);

    return c1 * Vector3f((1.f - alpha) * (1.f - beta))
        + c2 * Vector3f(alpha * (1.f - beta))
        + c3 * Vector3f((1.f - alpha) * beta)
        + c4 * Vector3f(alpha * beta);
}

void Image::splatPixel(const int pixelX, const int pixelY, const Vector3f& color) {
    pixelData[pixelY * resolution.x + pixelX] = color;
}

void Image::writeToFile(const std::string& filename) const {
    std::ofstream imageFile(filename);

    if (!imageFile.is_open()) {
        std::cerr << "Could not open file " << filename << ", cannot save the image to disk" << std::endl;
        return;
    }

    imageFile << "P3\n" << resolution.x << " " << resolution.y << "\n" << static_cast<int>(MAX_COLOR) << "\n";

    for (int i = 0; i < resolution.x * resolution.y; i++) {
        const int r = static_cast<int>(MAX_COLOR * pixelData[i].x);
        const int g = static_cast<int>(MAX_COLOR * pixelData[i].y);
        const int b = static_cast<int>(MAX_COLOR * pixelData[i].z);
        imageFile << r << " " << g << " " << b << "\n";
    }

    imageFile.close();
}

const std::vector<Vector3f>& Image::getPixelData() const {
    return pixelData;
}

std::vector<Vector3f>& Image::getPixelData() {
    return pixelData;
}