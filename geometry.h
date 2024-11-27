//
// Created by Matej on 11/21/2024.
//
#pragma once

#include <cmath>
#include <string>


constexpr float pi = std::numbers::pi;
constexpr float epsilon = std::numeric_limits<float>::epsilon();

inline float radians(const float degrees) {
    return static_cast<float>(degrees * std::numbers::pi / 180.f);
}

template<template <typename> class T, typename V>
class Tuple2 {
public:
    Tuple2() = default;
    Tuple2(V x, V y) {
        this->x = x;
        this->y = y;
    }
    explicit Tuple2(V value) {
        this->x = value;
        this->y = value;
    }
    explicit Tuple2(const T<V>& other) {
        this->x = other.x;
        this->y = other.y;
    }

    bool operator==(const T<V>& other) const {
        return x == other.x && y == other.y;
    }
    bool operator!=(const T<V>& other) const {
        return !(*this == other);
    }

    T<V>& operator=(const T<V>& other) {
        this->x = other.x;
        this->y = other.y;
        return static_cast<T<V>&>(*this);
    };
    T<V> operator+(const T<V>& other) const {
        return T<V>(x + other.x, y + other.y);
    }
    T<V>& operator+=(const T<V>& other) {
        x += other.x;
        y += other.y;
        return static_cast<T<V>&>(*this);
    }
    T<V> operator-(const T<V>& other) const {
        return T<V>(x - other.x, y - other.y);
    }
    T<V>& operator-=(const T<V>& other) {
        x -= other.x;
        y -= other.y;
        return static_cast<T<V>&>(*this);
    }
    T<V> operator*(const T<V>& other) const {
        return T<V>(x * other.x, y * other.y);
    }
    T<V>& operator*=(const T<V>& other) {
        x *= other.x;
        y *= other.y;
        return static_cast<T<V>&>(*this);
    }
    T<V> operator/(const T<V>& other) const {
        return T<V>(x / other.x, y / other.y);
    }
    T<V>& operator/=(const T<V>& other) {
        x /= other.x;
        y /= other.y;
        return static_cast<T<V>&>(*this);
    }

    T<V> operator-() const {
        return T<V>(-x, -y);
    }
    T<V> operator+(V value) const {
        return T<V>(x + value, y + value);
    }
    T<V>& operator+=(V value) {
        x += value;
        y += value;
        return static_cast<T<V>&>(*this);
    }
    T<V> operator-(V value) const {
        return T<V>(x - value, y - value);
    }
    T<V>& operator-=(V value) {
        x -= value;
        y -= value;
        return static_cast<T<V>&>(*this);
    }
    T<V> operator*(V value) const {
        return T<V>(x * value, y * value);
    }
    T<V>& operator*=(V value) {
        x *= value;
        y *= value;
        return static_cast<T<V>&>(*this);
    }
    T<V> operator/(V value) const {
        return T<V>(x / value, y / value);
    }
    T<V>& operator/=(V value) {
        x /= value;
        y /= value;
        return static_cast<T<V>&>(*this);
    }

    V operator[](const int i) const {
        return i == 0 ? x : y;
    }
    V& operator[](const int i) {
        return i == 0 ? x : y;
    }

    [[nodiscard]] std::string toString() const {
        return "(" + std::to_string(x) + ", " + std::to_string(y) + ")";
    }

    V x, y;
};

template <template <class> class T, typename V>
T<V> operator+(const V value, Tuple2<T, V> tuple) {
    return tuple + value;
}

template <template <class> class T, typename V>
T<V> operator-(const V value, Tuple2<T, V> tuple) {
    return tuple - value;
}

template <template <class> class T, typename V>
T<V> operator*(const V value, Tuple2<T, V> tuple) {
    return tuple * value;
}

template <template <class> class T, typename V>
T<V> lerp(V delta, const Tuple2<T, V>& start, const Tuple2<T, V>& end) {
    return start + delta * (end - start);
}

template <template <class> class T, typename V>
T<V> min(const Tuple2<T, V>& first, const Tuple2<T, V>& second) {
    return T<V>(std::min(first.x, second.x), std::min(first.y, second.y));
}

template <template <class> class T, typename V>
T<V> max(const Tuple2<T, V>& first, const Tuple2<T, V>& second) {
    return T<V>(std::max(first.x, second.x), std::max(first.y, second.y));
}

template <typename T>
class Vector2 : public Tuple2<Vector2, T> {
public:
    using Tuple2<Vector2, T>::x;
    using Tuple2<Vector2, T>::y;

    Vector2() = default;
    Vector2(T x, T y) : Tuple2<Vector2, T>(x, y) {}
    explicit Vector2(T value) : Tuple2<Vector2, T>(value) {}
    Vector2(const Vector2& other) : Tuple2<Vector2, T>(other) {}

    [[nodiscard]] T lengthSquared() const {
        return x * x + y * y;
    }
    [[nodiscard]] T length() const {
        return std::sqrt(lengthSquared());
    }
    [[nodiscard]] T dot(const Vector2& other) const {
        return x * other.x + y * other.y;
    }
    [[nodiscard]] Vector2 normalized() const {
        return *this / length();
    }
};

using Vector2i = Vector2<int>;
using Vector2f = Vector2<float>;

template <typename T>
class Point2 : public Tuple2<Point2, T> {
public:
    using Tuple2<Point2, T>::x;
    using Tuple2<Point2, T>::y;

    Point2() = default;
    Point2(T x, T y) : Tuple2<Point2, T>(x, y) {}
    explicit Point2(T value) : Tuple2<Point2, T>(value) {}
    Point2(const Point2& other) : Tuple2<Point2, T>(other.x, other.y) {}

    Vector2<T> operator-(const Point2& other) const  {
        return Vector2<T>(x - other.x, y - other.y);
    }

    [[nodiscard]] T distanceSquared(const Point2& other) const {
        return (*this - other).lengthSquared();
    }

    [[nodiscard]] T distance(const Point2& other) const {
        return (*this - other).length();
    }
};

using Point2i = Point2<int>;
using Point2f = Point2<float>;

template<template <typename> class T, typename V>
class Tuple3 {
public:
    Tuple3() = default;
    Tuple3(V x, V y, V z) {
        this->x = x;
        this->y = y;
        this->z = z;
    }
    explicit Tuple3(V value) {
        this->x = value;
        this->y = value;
        this->z = value;
    }
    explicit Tuple3(const T<V>& other) {
        this->x = other.x;
        this->y = other.y;
        this->z = other.z;
    }

    bool operator==(const T<V>& other) const {
        return x == other.x && y == other.y;
    }
    bool operator!=(const T<V>& other) const {
        return !(*this == other);
    }

    T<V>& operator=(const T<V>& other) {
        this->x = other.x;
        this->y = other.y;
        this->z = other.z;
        return static_cast<T<V>&>(*this);
    };
    T<V> operator+(const T<V> &other) const {
        return T<V>(x + other.x, y + other.y, z + other.z);
    }
    T<V>& operator+=(const T<V>& other) {
        x += other.x;
        y += other.y;
        z += other.z;
        return static_cast<T<V>&>(*this);
    }
    T<V> operator-(const T<V>& other) const {
        return T<V>(x - other.x, y - other.y, z - other.z);
    }
    T<V>& operator-=(const T<V>& other) {
        x -= other.x;
        y -= other.y;
        z -= other.z;
        return static_cast<T<V>&>(*this);
    }
    T<V> operator*(const T<V>& other) const {
        return T<V>(x * other.x, y * other.y, z * other.z);
    }
    T<V>& operator*=(const T<V>& other) {
        x *= other.x;
        y *= other.y;
        z *= other.z;
        return static_cast<T<V>&>(*this);
    }
    T<V> operator/(const T<V>& other) const {
        return T<V>(x / other.x, y / other.y, z / other.z);
    }
    T<V>& operator/=(const T<V>& other) {
        x /= other.x;
        y /= other.y;
        z /= other.z;
        return static_cast<T<V>&>(*this);
    }

    T<V> operator-() const {
        return T<V>(-x, -y, -z);
    }
    T<V> operator+(V value) const {
        return T<V>(x + value, y + value, z + value);
    }
    T<V> operator-(V value) const {
        return T<V>(x - value, y - value, z - value);
    }
    T<V> operator*(V value) const {
        return T<V>(x * value, y * value, z * value);
    }
    T<V> operator/(V value) const {
        return T<V>(x / value, y / value, z / value);
    }
    T<V>& operator+=(V value) {
        x += value;
        y += value;
        z += value;
        return static_cast<T<V>&>(*this);
    }
    T<V>& operator-=(V value) {
        x -= value;
        y -= value;
        z -= value;
        return static_cast<T<V>&>(*this);
    }
    T<V>& operator*=(V value) {
        x *= value;
        y *= value;
        z *= value;
        return static_cast<T<V>&>(*this);
    }
    T<V>& operator/=(V value) {
        x /= value;
        y /= value;
        z /= value;
        return static_cast<T<V>&>(*this);
    }

    V operator[](const int i) const {
        return i == 0 ? x : i == 1 ? y : z;
    }
    V& operator[](const int i) {
        return i == 0 ? x : i == 1 ? y : z;
    }

    [[nodiscard]] std::string toString() const {
        return "(" + std::to_string(x) + ", " + std::to_string(y) + std::string(", ") + std::to_string(z) + ")";
    }

    V x, y, z;
};

template <template <class> class T, typename V>
T<V> operator+(const V value, Tuple3<T, V> tuple) {
    return tuple + value;
}

template <template <class> class T, typename V>
T<V> operator-(const V value, Tuple3<T, V> tuple) {
    return tuple - value;
}

template <template <class> class T, typename V>
T<V> operator*(const V value, Tuple3<T, V> tuple) {
    return tuple * value;
}

template <template <class> class T, typename V>
T<V> lerp(V delta, const Tuple3<T, V>& start, const Tuple3<T, V>& end) {
    return start + delta * (end - start);
}

template <template <class> class T, typename V>
T<V> min(const Tuple3<T, V>& first, const Tuple3<T, V>& second) {
    return T<V>(std::min(first.x, second.x), std::min(first.y, second.y), std::min(first.z, second.z));
}

template <template <class> class T, typename V>
T<V> max(const Tuple3<T, V>& first, const Tuple3<T, V>& second) {
    return T<V>(std::max(first.x, second.x), std::max(first.y, second.y), std::max(first.z, second.z));
}

template <typename T>
class Vector3 : public Tuple3<Vector3, T> {
public:
    using Tuple3<Vector3, T>::x;
    using Tuple3<Vector3, T>::y;
    using Tuple3<Vector3, T>::z;

    Vector3() = default;
    Vector3(T x, T y, T z) : Tuple3<Vector3, T>(x, y, z) {}
    explicit Vector3(T value) : Tuple3<Vector3, T>(value) {}
    Vector3(const Vector3& other) : Tuple3<Vector3, T>(other) {}

    [[nodiscard]] T lengthSquared() const {
        return x * x + y * y + z * z;
    }
    [[nodiscard]] T length() const {
        return std::sqrt(lengthSquared());
    }
    [[nodiscard]] T dot(const Vector3& other) const {
        return x * other.x + y * other.y + z * other.z;
    }
    [[nodiscard]] Vector3 cross(const Vector3& other) const {
        Vector3 result;
        result.x = y * other.z - z * other.y;
        result.y = z * other.x - x * other.z;
        result.z = x * other.y - y * other.x;
        return result;
    }
    [[nodiscard]] Vector3 normalized() const {
        return *this / length();
    }
};

using Vector3i = Vector3<int>;
using Vector3f = Vector3<float>;

template <typename T>
class Point3 : public Tuple3<Point3, T> {
public:
    using Tuple3<Point3, T>::x;
    using Tuple3<Point3, T>::y;
    using Tuple3<Point3, T>::z;

    Point3() = default;
    Point3(T x, T y, T z) : Tuple3<Point3, T>(x, y, z) {}
    explicit Point3(T value) : Tuple3<Point3, T>(value) {}
    Point3(const Point3& other) : Tuple3<Point3, T>(other.x, other.y, other.z) {}
    explicit Point3(const Vector3<T>& other) {
        this->x = other.x;
        this->y = other.y;
        this->z = other.z;
    }

    [[nodiscard]] T distanceSquared(const Point3& other) const {
        return (*this - other).lengthSquared();
    }
    [[nodiscard]] T distance(const Point3& other) const {
        return (*this - other).length();
    }
};

template <typename T>
Point3<T> operator+(const Point3<T>& point, const Vector3<T>& other) {
    return Point3<T>(point.x + other.x, point.y + other.y, point.z + other.z);
}

template <typename T>
Point3<T>& operator+=(Point3<T>& point, const Vector3<T>& other) {
    point.x += other.x;
    point.y += other.y;
    point.z += other.z;
    return point;
}

template <typename T>
Point3<T> operator-(const Point3<T>& point, const Vector3<T>& vector) {
    return Point3(point.x - vector.x, point.y - vector.y, point.z - vector.z);
}

template <typename T>
Point3<T>& operator-=(const Point3<T>& point, const Vector3<T>& vector) {
    point.x -= vector.x;
    point.y -= vector.y;
    point.z -= vector.z;
    return *point;
}

using Point3i = Point3<int>;
using Point3f = Point3<float>;

template <typename T>
Vector3<T> operator-(const Point3<T>& point, const Point3<T>& other) {
    return Vector3<T>(point.x - other.x, point.y - other.y, point.z - other.z);
}

template <typename T>
Vector3<T> operator-(const Vector3<T>& vector, const Point3<T>& other) {
    return Vector3<T>(vector.x - other.x, vector.y - other.y, vector.z - other.z);
}

template <typename T>
class Normal3 : public Tuple3<Normal3, T> {
public:
    using Tuple3<Normal3, T>::x;
    using Tuple3<Normal3, T>::y;
    using Tuple3<Normal3, T>::z;

    Normal3() = default;
    Normal3(T x, T y, T z) : Tuple3<Normal3, T>(x, y, z) {
        const auto length = std::sqrt(x * x + y * y + z * z);
        this->x /= length;
        this->y /= length;
        this->z /= length;
    }
    explicit Normal3(T value) : Tuple3<Normal3, T>(value) {
        const auto length = std::sqrt(3) * value;
        this->x /= length;
        this->y /= length;
        this->z /= length;
    }
    Normal3(const Normal3& other) : Tuple3<Normal3, T>(other.x, other.y, other.z) {}
    explicit Normal3(const Vector3<T>& vector) {
        const auto normal = vector.normalized();
        x = normal.x;
        y = normal.y;
        z = normal.z;
    }

    Vector3<T> operator*(T value) const {
        return Vector3<T>(x * value, y * value, z * value);
    }

    [[nodiscard]] T dot(const Normal3& other) const {
        return x * other.x + y * other.y + z * other.z;
    }

    [[nodiscard]] T dot(const Vector3<T>& other) const {
        return x * other.x + y * other.y + z * other.z;
    }

    [[nodiscard]] Normal3 faceForward(const Vector3<T>& vector) const {
        return dot(vector) < 0.f ? -*this : *this;
    }
};

template <typename T>
Vector3<T> operator*(T value, const Normal3<T>& normal) {
    return Vector3<T>(normal.x * value, normal.y * value, normal.z * value);
}

template <typename T>
Vector3<T> operator+(const Normal3<T>& normal, const Vector3<T>& other) {
    return Vector3<T>(normal.x + other.x, normal.y + other.y, normal.z + other.z);
}

template <typename T>
Vector3<T> operator+(const Vector3<T>& vector, const Normal3<T>& other) {
    return Vector3<T>(vector.x + other.x, vector.y + other.y, vector.z + other.z);
}

template <typename T>
Vector3<T> operator-(const Normal3<T>& normal, const Vector3<T>& other) {
    return Vector3<T>(normal.x - other.x, normal.y - other.y, normal.z - other.z);
}

template <typename T>
Vector3<T> operator-(const Vector3<T>& vector, const Normal3<T>& other) {
    return Vector3<T>(vector.x - other.x, vector.y - other.y, vector.z - other.z);
}

using Normal3f = Normal3<float>;

template <typename T>
class AxisAlignedBoundingBox {
public:
    AxisAlignedBoundingBox() {
        T min = std::numeric_limits<T>::lowest();
        T max = std::numeric_limits<T>::max();
        this->minPoint = Point3<T>(max, max, max);
        this->maxPoint = Point3<T>(min, min, min);
    }
    AxisAlignedBoundingBox(const Point3<T>& minPoint, const Point3<T>& maxPoint) {
        this->minPoint = min(minPoint, maxPoint);
        this->maxPoint = max(minPoint, maxPoint);
    }
    explicit AxisAlignedBoundingBox(const Point3<T>& point) {
        this->minPoint = point;
        this->maxPoint = point;
    }

    bool isEmpty() {
        return minPoint.x >= maxPoint.x || minPoint.y >= maxPoint.y || minPoint.z >= maxPoint.z;
    }
    bool operator==(const AxisAlignedBoundingBox& other) const {
        return minPoint == other.minPoint && maxPoint == other.maxPoint;
    }
    bool operator!=(const AxisAlignedBoundingBox& other) const {
        return minPoint != other.minPoint || maxPoint != other.maxPoint;
    }

    AxisAlignedBoundingBox operator+(const AxisAlignedBoundingBox& other) {
        AxisAlignedBoundingBox result;
        result.minPoint = min<Point3, float>(minPoint, other.minPoint);
        result.maxPoint = max<Point3, float>(maxPoint, other.maxPoint);
        return result;
    }
    AxisAlignedBoundingBox& operator+=(const AxisAlignedBoundingBox& other) {
        this->minPoint = min<Point3, float>(minPoint, other.minPoint);
        this->maxPoint = max<Point3, float>(maxPoint, other.maxPoint);
        return *this;
    }
    AxisAlignedBoundingBox operator+(const Point3<T>& other) {
        AxisAlignedBoundingBox result;
        result.minPoint = min<Point3, float>(minPoint, other);
        result.maxPoint = max<Point3, float>(maxPoint, other);
        return result;
    }
    AxisAlignedBoundingBox& operator+=(const Point3<T>& other) {
        this->minPoint = min<Point3, float>(minPoint, other);
        this->maxPoint = max<Point3, float>(maxPoint, other);
        return *this;
    }

    Point3<T> operator[](const int index) const {
        return index == 0 ? minPoint : maxPoint;
    }
    Point3<T>& operator[](const int index) {
        return index == 0 ? minPoint : maxPoint;
    }
    Point3<T> corner(const int corner) {
        return Point3<T>(
            (*this)[corner & 1].x,
            (*this)[corner & 2 ? 1 : 0].y,
            (*this)[corner & 4 ? 1 : 0].z);
    }

    bool contains(const Point3<T>& point) const {
        return point.x >= minPoint.x && point.x <= maxPoint.x
            && point.y >= minPoint.y && point.y <= maxPoint.y
            && point.z >= minPoint.z && point.z <= maxPoint.z;
    }
    [[nodiscard]] int maxDimension() const {
        const auto diagonal = this->maxPoint - this->minPoint;
        if (diagonal.x > diagonal.y && diagonal.x > diagonal.z) {
            return 0;
        }
        if (diagonal.y > diagonal.z) {
            return 1;
        }
        return 2;
    }
    [[nodiscard]] float volume() const {
        const auto diagonal = this->maxPoint - this->minPoint;
        return diagonal.x * diagonal.y * diagonal.z;
    }
    [[nodiscard]] float surfaceArea() const {
        const auto diagonal = this->maxPoint - this->minPoint;
        return 2.f * (diagonal.x * diagonal.y + diagonal.x * diagonal.z + diagonal.y * diagonal.z);
    }

    Point3<T> minPoint, maxPoint;
};

using Bounds = AxisAlignedBoundingBox<float>;

class Quaternion {
public:
    Quaternion() = default;
    Quaternion(const float x, const float y, const float z, const float w) {
        this->u = Vector3f(x, y, z);
        this->w = w;
    }
    Quaternion(const Vector3f& u, const float w) {
        this->u = u;
        this->w = w;
    }

    bool operator==(const Quaternion& other) const {
        return u == other.u && w == other.w;
    }
    bool operator!=(const Quaternion& other) const {
        return u != other.u || w != other.w;
    }

    Quaternion operator+(const Quaternion& other) const {
        return { u + other.u, w + other.w };
    }
    Quaternion& operator+=(const Quaternion& other) {
        this->u += other.u;
        this->w += other.w;
        return *this;
    }
    Quaternion operator-(const Quaternion& other) const {
        return { u - other.u, w - other.w };
    }
    Quaternion& operator-=(const Quaternion& other) {
        this->u -= other.u;
        this->w -= other.w;
        return *this;
    }
    Quaternion operator+(const float value) const {
        return { u + value, w + value };
    }
    Quaternion& operator+=(const float value) {
        this->u += value;
        this->w += value;
        return *this;
    }
    Quaternion operator-(const float value) const {
        return { u - value, w - value };
    }
    Quaternion& operator-=(const float value) {
        this->u -= value;
        this->w -= value;
        return *this;
    }
    Quaternion operator*(const float value) const {
        return { u * value, w * value };
    }
    Quaternion& operator*=(const float value) {
        this->u *= value;
        this->w *= value;
        return *this;
    }
    Quaternion operator/(const float value) const {
        return { u / value, w / value };
    }
    Quaternion& operator/=(const float value) {
        this->u /= value;
        this->w /= value;
        return *this;
    }
    Quaternion operator-() const {
        return { -u, -w };
    }

    [[nodiscard]] float dot(const Quaternion& other) const {
        return u.dot(other.u) + w * other.w;
    }
    [[nodiscard]] float lengthSquared() const {
        return dot(*this);
    }
    [[nodiscard]] float length() const {
        return std::sqrt(lengthSquared());
    }
    [[nodiscard]] Quaternion normalized() const {
        return *this / length();
    }
    [[nodiscard]] Quaternion conjugate() const {
        return { -u.x, -u.y, -u.z, w };
    }
    [[nodiscard]] Quaternion inverse() const {
        return conjugate() / lengthSquared();
    }

    Vector3f operator*(const Vector3f& other) const {
        // GLM implementation
        const auto uv = u.cross(other);
        return other + (uv * w + u.cross(uv)) * 2.f;
    }
    Point3f operator*(const Point3f& other) const {
        return Point3f(*this * Vector3f(other.x, other.y, other.z));
    }
    Normal3f operator*(const Normal3f& other) const {
        return Normal3f(*this * Vector3f(other.x, other.y, other.z));
    }

    // https://en.wikipedia.org/wiki/Conversion_between_quaternions_and_Euler_angles#Source_code_2
    [[nodiscard]] Vector3f toEulerAngles() const {
        Vector3f result;

        const auto sinX_cosY = 2.f * (w * u.x + u.y * u.z);
        const auto cosX_cosY = 1.f - 2.f * (u.x * u.x + u.y * u.y);
        result.x = std::atan2(sinX_cosY, cosX_cosY);

        const auto sinY = std::sqrt(1.f + 2.f * (w * u.y - u.x * u.z));
        const auto cosY = std::sqrt(1.f - 2.f * (w * u.y - u.x * u.z));
        result.y = 2.f * std::atan2(sinY, cosY) - static_cast<float>(std::numbers::pi) / 2.f;

        const auto sinZ_cosY = 2.f * (w * u.z + u.x * u.y);
        const auto cosZ_cosY = 1.f - 2.f * (u.y * u.y + u.z * u.z);
        result.z = std::atan2(sinZ_cosY, cosZ_cosY);

        return result;
    }

    [[nodiscard]] std::string toString() const {
        return "(" + std::to_string(u.x) + ", " + std::to_string(u.y) + ", " + std::to_string(u.z) + ", " + std::to_string(w) + ")";
    }

    static Quaternion fromAxisAngle(const Vector3f& axis, float angle);
    static Quaternion fromEulerAngles(const Vector3f& eulerAngles);
    static Quaternion lookAt(const Point3f& from, const Point3f& to, const Vector3f& up = { 0.f, 1.f, 0.f });
    static Quaternion align(const Vector3f& v1, const Vector3f& v2);

    Vector3f u {};
    float w = 1.f;
};

inline Quaternion operator*(const float value, const Quaternion& quaternion) {
    return quaternion * value;
}