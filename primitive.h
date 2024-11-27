//
// Created by Matej on 11/25/2024.
//
#pragma once
#include <memory>
#include <numeric>
#include <optional>
#include <utility>
#include <variant>
#include <vector>
#include <algorithm>
#include <queue>

#include "geometry.h"
#include "material.h"
#include "ray.h"


class Primitive {
public:
    virtual ~Primitive() = default;

    [[nodiscard]] std::optional<RayInteraction> intersect(const Ray& ray) const {
        if (!ray.intersects(bounds)) return {};
        return intersectShape(ray);
    }

    Point3f center {};
    Bounds bounds {};
    Material material {};
protected:
    Primitive(const Point3f& center, Bounds bounds, const Material& material):
        center(center), bounds(std::move(bounds)), material(material) {}

    [[nodiscard]] virtual std::optional<RayInteraction> intersectShape(const Ray& ray) const = 0;
};

class Sphere final : public Primitive {
public:
    Sphere(const Point3f& center, const float radius, const Material& material):
        Primitive(center, Bounds(center - radius, center + radius), material), radius(radius) {}

    [[nodiscard]] std::optional<RayInteraction> intersectShape(const Ray& ray) const override {
        const auto l = ray.origin - center;
        const auto angle = ray.direction.dot(l);

        const auto delta = angle * angle - std::pow(l.length(), 2.f) + radius * radius;

        if (delta < 0) {
            return {};
        }

        const auto t1 = -angle + std::sqrt(delta);
        const auto t2 = -angle - std::sqrt(delta);
        float tHit = 0.f;

        if (t1 > 0 && t2 > 0) {
            tHit = std::min(t1, t2);
        } else if (t1 > 0) {
            tHit = t1;
        } else if (t2 > 0) {
            tHit = t2;
        } else {
            return {};
        }

        RayInteraction result {};
        result.intersection = ray(tHit);
        result.tHit = tHit;
        result.normal = Normal3f(result.intersection - center).faceForward(-ray.direction);
        result.material = material;

        // Compute UV
        auto u = .5f * std::atan2(result.normal.y, result.normal.x) / pi;
        if (u < 0.f) u += 1.f;
        const auto v = .5f * (result.normal.z - radius) / radius;
        result.uv = { u, v };

        return result;
    }

    float radius;
};

class Cylinder final : public Primitive {
public:
    Cylinder(const Point3f& center, const Vector3f& axis, const float radius, const float height, const Material& material):
        Primitive(center, Bounds(), material),
        axis(axis),
        axisTransform(Quaternion::align({ 0.f, 0.f, 1.f}, axis)),
        inverseAxisTransform(axisTransform.inverse()),
        radius(radius),
        height(height) {
        // Shrink the bounds as tightly as possible to fit an axis-aligned box around the rotated cylinder
        std::vector<Point3f> points {};
        for (int i = 0; i < 6; i++) {
            const Vector3f offset { (i & 1 << 1 ? 1.f : -1.f) * radius, (i & 1 ? 1.f : -1.f) * radius, (i & 1 << 2 ? 1.f : -1.f) * height };
            points.push_back(center + axisTransform * offset);
        }
        this->bounds = std::accumulate(points.begin(), points.end(), bounds);
    }

    [[nodiscard]] std::optional<RayInteraction> intersectShape(const Ray& ray) const override {
        // Transform coordinates such that the cylinder's axis becomes the local z-axis
        const auto transformedOrigin = Point3f(inverseAxisTransform * (ray.origin - center));
        const auto transformedDirection = inverseAxisTransform * ray.direction;

        // Solve quadratic for t: (x_o + t*x_dir)^2 + (y_o + t*y_dir)^2 = r^2
        const auto a = transformedDirection.x * transformedDirection.x + transformedDirection.y * transformedDirection.y;
        const auto b = 2 * transformedOrigin.x * transformedDirection.x + 2 * transformedOrigin.y * transformedDirection.y;
        const auto c = transformedOrigin.x * transformedOrigin.x + transformedOrigin.y * transformedOrigin.y - radius * radius;
        const auto discriminant = b * b - 4 * a * c;

        if (discriminant < 0) return {};

        auto t1 = -.5f * (b - std::sqrt(discriminant)) / a;
        auto t2 = -.5f * (b + std::sqrt(discriminant)) / a;
        // Keep the smallest intersection at t1
        if (t1 > t2) {
            std::swap(t1, t2);
        }
        const auto z1 = transformedOrigin.z + t1 * transformedDirection.z;
        const auto z2 = transformedOrigin.z + t2 * transformedDirection.z;
        const auto zMin = -height;
        const auto zMax = height;
        float tHit = 0.f;

        // Both intersections are in front of the camera (the closest one is t1)
        if (t1 > 0 && t2 > 0) {
            // Either one is within z-axis bounds
            if (z1 > zMin && z1 < zMax || z2 > zMin && z2 < zMax) {
                // If the closest one is outside, set the intersection with the appropriate cap
                if (z1 <= zMin) {
                    tHit = (zMin - transformedOrigin.z) / transformedDirection.z;
                } else if (z1 >= zMax) {
                    tHit = (zMax - transformedOrigin.z) / transformedDirection.z;
                } else {
                    // Otherwise proceed with the closest one
                    tHit = t1;
                }
            } else {
                // None are within z-axis bounds, checking orderings z(1/2) < zMin <= z(2/1) < zMax
                if (z1 < zMin && z2 >= zMin || z2 < zMin && z1 >= zMin) {
                    // If the closest one is outside, set the intersection with the appropriate cap
                    if (z1 < zMin) {
                        tHit = (zMin - transformedOrigin.z) / transformedDirection.z;
                    } else {
                        tHit = (zMax - transformedOrigin.z) / transformedDirection.z;
                    }
                // Checking orderings zMin < z(2/1) <= zMax < z(1/2)
                } else if (z1 > zMax && z2 <= zMax || z2 > zMax && z1 <= zMax) {
                    // If the closest one is outside, set the intersection with the appropriate cap
                    if (z1 > zMax) {
                        tHit = (zMax - transformedOrigin.z) / transformedDirection.z;
                    } else {
                        tHit = (zMin - transformedOrigin.z) / transformedDirection.z;
                    }
                // Discard all other intersections
                } else {
                    return {};
                }
            }
        // If the camera is inside the infinite cylinder - t1 is behind the origin
        } else if (t1 > 0 || t2 > 0) {
            // The negative intersection is before/at zMin
            if (z2 > zMin && z1 <= zMin) {
                tHit = (zMin - transformedOrigin.z) / transformedDirection.z;
            // The negative intersection is at/after zMax
            } else if (z2 < zMax && z1 >= zMax) {
                tHit = (zMax - transformedOrigin.z) / transformedDirection.z;
            // Discard the rest
            } else {
                return {};
            }
        // Kept for correctness, however, this case should not occur at all
        } else {
            return {};
        }

        RayInteraction result {};
        result.intersection = ray(tHit);
        result.tHit = tHit;
        result.material = material;

        const auto relative = transformedOrigin + tHit * transformedDirection;

        // Normal pointing away from the bottom cap
        if (std::abs(relative.z - zMin) <= .0001f) {
            result.normal = axisTransform * Normal3f(0.f, 0.f, -1.f);
        // Normal pointing away from the top cap
        } else if (std::abs(relative.z - zMax) <= .0001f) {
            result.normal = axisTransform * Normal3f(0.f, 0.f, 1.f);
        // Normal pointing away from the side
        } else {
            // Not really necessary to check this (correctness)
            if (relative.z > zMin && relative.z < zMax) {
                result.normal = axisTransform * Normal3f(relative.x, relative.y , 0.f);
            }
        }
        result.normal = result.normal.faceForward(-ray.direction);

        // Compute UV
        auto u = .5f * std::atan2f(relative.y, relative.x) / pi;
        if (u < 0.f) u += 1.f;
        float v;
        if (std::abs(relative.z - zMin) <= .0001f || std::abs(relative.z - zMax) <= .0001f) {
            v = std::sqrt(relative.x * relative.x + relative.y * relative.y) / radius;
        } else {
            v = .5f * (zMax - relative.z) / height;
        }
        result.uv = { u, v };

        return result;
    }

    Vector3f axis;
    Quaternion axisTransform;
    Quaternion inverseAxisTransform;
    float radius;
    float height;
};

class Triangle final : public Primitive {
public:
    // TODO: meshes
    Triangle(
        const Point3f& a, const Point3f& b, const Point3f& c,
        const Normal3f& n1, const Normal3f& n2, const Normal3f& n3,
        const Point2f& uv1, const Point2f& uv2, const Point2f& uv3,
        const Material& material):
        Primitive((a + b + c) / 3.f, Bounds(a, b) + c, material),
        a(Point3f(a - center)), b(Point3f(b - center)), c(Point3f(c - center)),
        n1(n1), n2(n2), n3(n3),
        uv1(uv1), uv2(uv2), uv3(uv3) {}

    Triangle(const Point3f& a, const Point3f& b, const Point3f& c, const Material& material):
        Triangle(
            a, b, c,
            Normal3f((b - a).cross(c - a)), Normal3f((b - a).cross(c - a)), Normal3f((b - a).cross(c - a)),
            Point2f(0.f), Point2f(0.f, 1.f), Point2f(1.f, 1.f),
            material) {}

    // https://en.wikipedia.org/wiki/M%C3%B6ller%E2%80%93Trumbore_intersection_algorithm#C++_implementation
    [[nodiscard]] std::optional<RayInteraction> intersectShape(const Ray& ray) const override {
        const auto e1 = b - a;
        const auto e2 = c - a;
        const auto dxe2 = ray.direction.cross(e2);
        const auto det = e1.dot(dxe2);

        if (std::abs(det) <= epsilon) return {};

        const auto invDet = 1.f / det;
        const auto s = ray.origin - center - a;
        const auto u = invDet * s.dot(dxe2);

        if ((u < 0.f && std::abs(u) > epsilon) || (u > 1.f && std::abs(u - 1.f) > epsilon)) return {};

        const auto sxe1 = s.cross(e1);
        const auto v = invDet * ray.direction.dot(sxe1);

        if ((v < 0.f && std::abs(v) > epsilon) || (u + v > 1.f && std::abs(u + v - 1.f) > epsilon)) return {};

        // At this stage we can compute t to find out where the intersection point is on the line.
        const auto tHit = invDet * e2.dot(sxe1);

        if (tHit < epsilon) return {};

        // Compute barycentric coordinates
        const auto b1 = .5f * (b - a).cross(ray(tHit) - center - a).length();
        const auto b2 = .5f * (a - c).cross(ray(tHit) - center - c).length();
        const auto b3 = 1.f - b2 - b1;

        RayInteraction result {};
        result.intersection = ray(tHit);
        result.tHit = tHit;
        result.normal = Normal3f(b1 * n1 + b2 * n2 + b3 * n3).faceForward(-ray.direction);
        result.uv = b1 * uv1 + b2 * uv2 + b3 * uv3;
        result.material = material;

        return result;
    }

    Point3f a;
    Point3f b;
    Point3f c;
    Normal3f n1;
    Normal3f n2;
    Normal3f n3;
    Point2f uv1;
    Point2f uv2;
    Point2f uv3;
};

struct BVHNode final {
    [[nodiscard]] bool isLeaf() const {
        return left == nullptr && right == nullptr;
    }

    Bounds bounds {};
    BVHNode* left = nullptr;
    BVHNode* right = nullptr;
    std::vector<int> primitives {};
};

struct SAHBin final {
    Bounds bounds;
    std::vector<int> primitives {};
};

constexpr int N_BINS = 12;

class BoundingVolumeHierarchy {
public:
    explicit BoundingVolumeHierarchy(std::vector<std::shared_ptr<Primitive>>& primitives) : primitives(primitives) {
        auto* root = new BVHNode;
        this->root = root;

        auto bounds = Bounds();
        auto centroids = std::vector<Point3f>();

        auto primitiveIndex = 0;
        for (const auto& primitive : primitives) {
            const auto centroid = .5f * (primitive->bounds.minPoint + primitive->bounds.maxPoint);
            centroids.push_back(centroid);
            bounds += primitive->bounds;
            this->root->primitives.push_back(primitiveIndex++);
        }

        this->root->bounds = bounds;
        subdivideNode(this->root, centroids);
    }

    [[nodiscard]] std::optional<RayInteraction> intersect(Ray& ray) const {
        // The queue element used in BVH traversal
        struct NodeIntersection {
            bool operator>(const NodeIntersection& other) const {
                return tHit > other.tHit;
            }
            BVHNode* node;
            float tHit;
        };

        std::priority_queue<NodeIntersection, std::vector<NodeIntersection>, std::greater<>> queue {}; // Create priority queue
        std::optional<RayInteraction> minInteraction {};

        const auto rootBounds = root->bounds;

        const auto inside = rootBounds.contains(ray.origin);
        if (const auto intersect = ray.intersects(rootBounds); intersect || inside) {

            if (inside) {
                // Push root to queue with key 0 (matters in case both origin and end are inside root AABB)
                queue.push({ root, 0.f });
            } else {
                // Push root to queue
                queue.push({ root, *intersect });
            }

            // Go through the tree until no better options are available
            while (!queue.empty()) {
                const auto [node, tHit] = queue.top();
                queue.pop();

                // The current node has potential for a better intersection
                if (tHit <= ray.time) {
                    if (node->isLeaf()) {
                        // Find minimal intersection in all primitives at the leaf
                        for (const auto& index : node->primitives) {
                            if (const auto& interaction = primitives[index]->intersect(ray)) {
                                if (ray.time > interaction->tHit) {
                                    ray.time = interaction->tHit;
                                    minInteraction = interaction;
                                }
                            }
                        }
                    } else {
                        // Node has left child
                        if (const auto& left = node->left) {
                            if (left->bounds.contains(ray.origin)) {
                                queue.push({ left, 0 }); // Insert the node as a primitive can be intersected at any time
                            } else if (const auto childIntersect = ray.intersects(left->bounds)) {
                                // Insert into queue only if the node has potential for a better intersection
                                if (childIntersect <= ray.time) {
                                    queue.push({ left, *childIntersect });
                                }
                            }
                        }
                        // Node has right child
                        if (const auto& right = node->right) {
                            if (right->bounds.contains(ray.origin)) {
                                queue.push({ right, 0 }); // Insert the node as a primitive can be intersected at any time
                            } else if (const auto childIntersect = ray.intersects(right->bounds)) {
                                // Insert into queue only if the node has potential for a better intersection
                                if (childIntersect <= ray.time) {
                                    queue.push({ right, *childIntersect });
                                }
                            }
                        }
                    }
                }
            }
        }

        return minInteraction;
        // std::optional<RayInteraction> minInteraction {};
        // for (const auto& index : primitives) {
        //     if (const auto& interaction = index->intersect(ray)) {
        //         if (ray.time > interaction->tHit) {
        //             ray.time = interaction->tHit;
        //             minInteraction = interaction;
        //         }
        //     }
        // }
        // return minInteraction;
    }

    void subdivideNode(BVHNode* node, std::vector<Point3f>& centroids) {

        if (node->primitives.size() <= 1) {
            nodes.push_back(*node);
            return;
        }

        const Bounds centroidBounds = std::accumulate(centroids.begin(), centroids.end(), Bounds());
        const auto splitDimension = centroidBounds.maxDimension();

        // Allocate splitting bins
        std::vector<SAHBin> bins(N_BINS);

        // Add primitives to their appropriate bin and extend their bounds
        for (int i = 0; i < node->primitives.size(); i++) {
            const float offset = (centroids[node->primitives[i]][splitDimension] - centroidBounds.minPoint[splitDimension]) / (centroidBounds.maxPoint - centroidBounds.minPoint)[splitDimension];
            const int index = std::clamp(N_BINS * static_cast<int>(offset), 0, N_BINS - 1);
            bins[index].bounds += primitives[i]->bounds;
            bins[index].primitives.push_back(i);
        }

        // Initialize approximate costs in linear time
        const int nSplits = N_BINS - 1;
        float costs[nSplits];

        // Forwards scans over splits
        int countBelow = 0;
        Bounds boundBelow;
        for (int i = 0; i < nSplits; i++) {
            boundBelow += bins[i].bounds;
            countBelow += static_cast<int>(bins[i].primitives.size());
            costs[i] += static_cast<float>(countBelow) * boundBelow.surfaceArea();
        }

        // Backwards scans over splits
        int countAbove = 0;
        Bounds boundAbove;
        for (int i = nSplits; i >= 1; i--) {
            boundAbove += bins[i].bounds;
            countAbove += static_cast<int>(bins[i].primitives.size());
            costs[i - 1] += static_cast<float>(countAbove) * boundAbove.surfaceArea();
        }

        auto splitIndex = -1;
        auto minCost = std::numeric_limits<float>::max();
        for (int i = 0; i < nSplits; i++) {
           if (costs[i] < minCost) {
               minCost = costs[i];
               splitIndex = i;
           }
        }

        if (static_cast<float>(node->primitives.size()) < .5f + minCost / node->bounds.surfaceArea()) {
            return;
        }

        const auto partition = std::partition(node->primitives.begin(), node->primitives.end(), [=](const int& child) {
            const auto offset = (centroids[child][splitDimension] - centroidBounds.minPoint[splitDimension]) / (centroidBounds.maxPoint - centroidBounds.minPoint)[splitDimension];
            const auto index = std::clamp(N_BINS * static_cast<int>(offset), 0, N_BINS - 1);
            return index <= splitIndex;
        });

        const auto middle = partition - node->primitives.begin();

        auto* left = new BVHNode;
        auto* right = new BVHNode;
        left->primitives = { node->primitives.begin(), node->primitives.begin() + middle };
        left->bounds = std::accumulate(left->primitives.begin(), left->primitives.end(), Bounds(), [&](Bounds bounds, const int child) {
            return bounds + primitives[child]->bounds;
        });
        right->primitives = { node->primitives.begin() + middle, node->primitives.end() };
        right->bounds = std::accumulate(right->primitives.begin(), right->primitives.end(), Bounds(), [&](Bounds bounds, const int child) {
            return bounds + primitives[child]->bounds;
        });

        node->left = left;
        node->right = right;
        node->primitives.clear();
        nodes.push_back(*node);

        subdivideNode(left, centroids);
        subdivideNode(right, centroids);
    }

    std::vector<std::shared_ptr<Primitive>>& primitives;
private:
    BVHNode* root;
    std::vector<BVHNode> nodes;
};


