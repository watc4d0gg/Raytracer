//
// Created by Matej on 11/21/2024.
//

#include "geometry.h"


Quaternion Quaternion::fromAxisAngle(const Vector3f& axis, const float angle) {
    return { axis.normalized() * std::sin(angle / 2.f), std::cos(angle / 2.f) };
}

// https://en.wikipedia.org/wiki/Conversion_between_quaternions_and_Euler_angles#Source_code
Quaternion Quaternion::fromEulerAngles(const Vector3f& eulerAngles) {
    Quaternion result;
    const auto cosX = std::cos(eulerAngles.x * .5f);
    const auto sinX = std::sin(eulerAngles.x * .5f);
    const auto cosY = std::cos(eulerAngles.y * .5f);
    const auto sinY = std::sin(eulerAngles.y * .5f);
    const auto cosZ = std::cos(eulerAngles.z * .5f);
    const auto sinZ = std::sin(eulerAngles.z * .5f);
    result.u = {
        sinX * cosY * cosZ - cosX * sinY * sinZ,
        cosX * sinY * cosZ + sinX * cosZ * sinZ,
        cosX * cosY * sinZ - sinX * sinY * cosZ,
    };
    result.w = cosX * cosY * cosZ + sinX * sinY * sinZ;
    return result;
}

// https://stackoverflow.com/a/52551983
Quaternion Quaternion::lookAt(const Point3f& from, const Point3f& to, const Vector3f& up) {
    const auto forward = (to - from).normalized();
    const auto side = forward.cross(up).normalized();
    const auto rotatedUp = side.cross(forward);
    const auto trace = side.x + rotatedUp.y + forward.z;

    if (trace > 0.f) {
        const float s = .5f * std::sqrt(trace + 1.f);
        return {
            (rotatedUp.z - forward.y) * s,
            (forward.x - side.z) * s,
            (side.y - rotatedUp.x) * s,
            .25f * s,
        };
    }

    if (side.x > rotatedUp.y && side.x > rotatedUp.z) {
        const float s = 2.f * std::sqrt(1.f + side.x - rotatedUp.y - forward.z);
        return {
            (rotatedUp.z - forward.y) / s,
            .25f * s,
            (rotatedUp.x + side.y) / s,
            (forward.x + side.z) / s,
        };
    }

    if (rotatedUp.y > forward.z) {
        const float s = 2.f * std::sqrt(1.f + rotatedUp.y - side.x - forward.z);
        return {
            (forward.x - side.z) / s,
            (rotatedUp.x + side.y) / s,
            .25f * s,
            (forward.y + rotatedUp.z) / s,
        };
    }

    const float s = 2.f * std::sqrt(1.f + forward.z - side.x - rotatedUp.y);
    return {
        (side.y - rotatedUp.x) / s,
        (forward.x + side.z) / s,
        (forward.y + rotatedUp.z) / s,
        .25f * s,
    };
}

Quaternion Quaternion::align(const Vector3f& v1, const Vector3f& v2) {
    const Quaternion result { v1.cross(v2), std::sqrt(std::pow(v1.length(), 2.f) * std::pow(v2.length(), 2.f)) + v1.dot(v2) };
    return result.normalized();
}
