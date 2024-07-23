#ifndef __MATHLIB_H__
#define __MATHLIB_H__

#pragma once

#include <cmath>
#include <math.h>
#include <random>
#include <cassert>
#include <iostream>
#include <algorithm>
#include <array>
#include <tuple>
#include <immintrin.h>

static const float INF = std::numeric_limits<float>::infinity();
static const float EPS = std::numeric_limits<float>::epsilon();
static const float PI  = 3.141592653589793115997963468544185161590576171875;

/********************************************************************************
*                               Vector Defination                               *
********************************************************************************/
// template parameter used here
template<int n, typename T> struct vec {
    T data[n] = {0};

    vec() = default;
    T & operator[] (const int i)       { assert(i >= 0 && i < n); return data[i]; }
    T   operator[] (const int i) const { assert(i >= 0 && i < n); return data[i]; }
    vec operator-() const { 
        for (int i = 0; i < n; data[i] = -data[i])
        return *this;
    }
    vec& operator+=(const vec<n,T> &v) {
        for (int i = 0; i < n; data[i] += v[i])
        return *this;
    }
    vec& operator*=(T t) {
        for (int i = 0; i < n; data[i] *= t)
        return *this;
    }
    vec& operator/=(T t) { return *this *= 1/t; }

    T norm2() const { return *this * *this; }
    T norm()  const { return std::sqrt(norm2()); }
    vec normalize() const { return (*this) * (1 / norm()); }
};

// overload operator <<, vector print
template<int n,typename T> std::ostream& operator<< (std::ostream& out, const vec<n,T>& v) {
    std::cout << "[ ";
    for (int i = 0; i < n; i++) std::cout << v[i] << " ";
    std::cout << ']';
    return out;
};

// overload operator +, vector - vector addition
template<int n, typename T>
inline vec<n,T> operator+ (const vec<n,T>& u, const vec<n,T>& v) {
    vec<n,T> res = u;
    for (int i = n; i--; res[i] += v[i]) ;
    return res;
};

// overload operator +, vector - scalar addition
template<int n, typename T>
inline vec<n,T> operator+ (const vec<n,T>& u, const T& k) {
    vec<n,T> res = u;
    for (int i = n; i--; res[i] += k) ;
    return res;
};

// overload operator -, vector - vector subtraction
template<int n, typename T>
inline vec<n,T> operator- (const vec<n,T>& u, const vec<n,T>& v) {
    vec<n,T> res = u;
    for (int i = n; i--; res[i] -= v[i]) ;
    return res;
};

// overload operator -, vector - scalar subtraction
template<int n, typename T>
inline vec<n,T> operator- (const vec<n,T>& u, const T x) {
    vec<n,T> res = u;
    for (int i = n; i--; res[i] -= x) ;
    return res;
};

// overload operator *, vector dot product
template<int n, typename T>
inline T operator* (const vec<n,T>& u, const vec<n,T>& v) {
    T res = 0;
    for (int i = n; i--; res += u[i] * v[i]) ;
    return res;
};

// overload operator *, vector - scalar multipication
template<int n, typename T>
inline vec<n,T> operator* (const T& k, const vec<n,T>& u) {
    vec<n,T> res = u;
    for (int i = n; i--; res[i] *= k) ;
    return res;
};

// overload operator *, vector - scalar multipication
template<int n, typename T>
inline vec<n,T> operator* (const vec<n,T>& u, const T& k) {
    return k * u;
};

// overload operator /, vector - scalar division
template<int n, typename T> 
inline vec<n,T> operator/ (const vec<n,T>& u, const T& k) {
    return (1 / k) * u;
};

// overload operator ==, vector == vector
template<int n, typename T>
inline bool operator== (const vec<n,T>& u, const vec<n,T>& v) {
    for (int i = 0; i < n; ++ i)
        if (u[i] != v[i]) return false;
    return true;
};

// vector element-wise product
template<int n, typename T> 
inline vec<n,T> mult(const vec<n,T> &v1, const vec<n,T> &v2) {
    vec<n,T> res;
    for (int i = n; i--; res[i] = v1[i] * v2[i]) ;
    return res;
}

template<int n1, int n2, typename T>
inline vec<n1,T> embed(const vec<n2,T>& v, T fill = 0) {
    vec<n1,T> res;
    for (int i = n1; i--; res[i] = (i < n2) ? v[i] : fill) ;
    return res;
}

template<int n1, int n2, typename T>
inline vec<n1,T> slice(const vec<n2,T> &v, const int m=0, const int n=n1-1) {
    assert(n2 > (n - m + 1) && n >= m);
    vec<n1,T> res;
    for (int i = 0, j = m; j <= n; i++, j++) res[i] = v[j];
    return res;
}

// compute the projection of u onto v
template<int n, typename T> 
inline vec<n,T> proj(const vec<n,T>& u, const vec<n,T>& v) {
    vec<n,T> res = u;
    vec<n,T> vNorm = (1 / v.norm()) * v;
    return (res *  vNorm) * vNorm;
}

// normalize the vector
template<int n, typename T> 
inline vec<n,T> normalize(vec<n,T> &&v) {
    float k = 1.f / v.norm();
    v *= k;
    return std::forward<vec<n,T>>(v);
}

namespace std {
    template<int n, typename T> struct hash<vec<n,T>> {
        size_t operator()(const vec<n,T> &v) const {
            size_t seed = 0;
            for (int i = 0; i < n; ++ i)
                hashCombine(seed, v[i]);
            return seed;
        }
    };

    template<typename T>
    void hashCombine(size_t &seed, const T &v) {
        std::hash<T> hasher;
        seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    }
}

// Specialization for vec<2,float>
template<> struct vec<2,float> {
    float data[2] ={0};
    float& x = data[0];
    float& y = data[1];

    vec() = default;
    vec(float x_, float y_): data{x_, y_} {}
    vec(const vec &other) {
        data[0] = other.x;
        data[1] = other.y;
    }
    float & operator[] (const int i)       { assert(i >= 0 && i < 2); return data[i]; }
    float   operator[] (const int i) const { assert(i >= 0 && i < 2); return data[i]; }
    vec& operator=(const vec& other) {
        data[0] = other.data[0];
        data[1] = other.data[1];
        return *this;
    }
    vec operator-() const { return {-x, -y}; }
    vec& operator+=(const vec<2,float> &v) {
        x += v.x;
        y += v.y;
        return *this;
    }
    vec& operator*=(float t) {
        x *= t;
        y *= t;
        return *this;
    }
    vec& operator/=(float t) { return *this *= 1/t; }

    float norm2() const { return x*x + y*y; }
    float norm()  const { return sqrt(norm2()); }
    vec normalize() const { return (*this) * (1 / norm()); }
};

// template<> struct vec<2,float> {
//     // float data[2] ={0};
//     float x{}, y{};

//     vec() = default;
//     vec(float x_, float y_): x(x_), y(y_) {}

//     float & operator[] (const int i)       { assert(i >= 0 && i < 2); return i == 0 ? x : y; }
//     float   operator[] (const int i) const { assert(i >= 0 && i < 2); return i == 0 ? x : y; }

//     vec operator-() const { return {-x, -y}; }
//     vec& operator+=(const vec<2,float> &v) {
//         x += v.x;
//         y += v.y;
//         return *this;
//     }
//     vec& operator*=(float t) {
//         x *= t;
//         y *= t;
//         return *this;
//     }
//     vec& operator/=(float t) { return *this *= 1/t; }

//     float norm2() const { return x*x + y*y; }
//     float norm()  const { return sqrt(norm2()); }
//     vec normalize() const { return (*this) * (1 / norm()); }
// };

// overload operator +, vec2 + vec2
inline vec<2,float> operator+(const vec<2,float> &u, const vec<2,float> &v) {
    return {u.x + v.x, u.y + v.y};
}

// overload operator +, vec2 + float
inline vec<2,float> operator+(const vec<2,float> &u, const float &k) {
    return {u.x + k, u.y + k};
}

// overload operator -, vec2 - vec2
inline vec<2,float> operator-(const vec<2,float> &u, const vec<2,float> &v) {
    return {u.x - v.x, u.y - v.y};
}

// overload operator -, vec2 - float
inline vec<2,float> operator-(const vec<2,float> &u, const float &k) {
    return {u.x - k, u.y - k};
}

// overload operator *, float * vec2
inline vec<2,float> operator*(const float &k, const vec<2,float> &u) {
    return {u.x * k, u.y * k};
}

// overload operator *, vec2 * float addition
inline vec<2,float> operator*(const vec<2,float> &u, const float &k) {
    return k * u;
}

// overload operator *, vec2 dot vec2
inline float operator*(const vec<2,float> &u, const vec<2,float> &v) {
    return u.x * v.x + u.y * v.y;
}

// overload operator /, vec2 / float
inline vec<2,float> operator/(const vec<2,float> &v, const float &k) {
    return (1/k) * v;
}

// overload operator ==, vec2 == vec2
inline bool operator==(const vec<2,float> &u, const vec<2,float> &v) {
    return u.x == v.x && u.y == v.y;
}

// vector compnent-wise product
inline vec<2,float> mult(const vec<2,float> &v1, const vec<2,float> &v2) {
    return {v1.x*v2.x, v1.y*v2.y};
}

// normalize the vector
inline vec<2,float> normalize(vec<2,float> &&v) {
    float k = 1.f / v.norm();
    v *= k;
    return std::forward<vec<2,float>>(v);
}

namespace std {
    template<> struct hash<vec<2,float>> {
        size_t operator()(const vec<2,float> &v) const {
            size_t seed = 0;
            hashCombine(seed, v.x);
            hashCombine(seed, v.y);
            return seed;
        }
    };
}

// // Specialization for vec<3,float>
// template<> struct vec<3,float> {
//     float data[3] = {0};
//     float& x = data[0];
//     float& y = data[1];
//     float& z = data[2];

//     vec() = default;
//     vec(float x_, float y_, float z_ = 0.0): data{x_, y_, z_} {}
//     vec(const vec &other) {
//         data[0] = other.data[0];
//         data[1] = other.data[1];
//         data[2] = other.data[2];
//     }
//     vec(const vec<2,float> &v, float z_ = 0.0) {
//         data[0] = v.data[0];
//         data[1] = v.data[1];
//         data[2] = z_;
//     }
//     float & operator[] (const int i)       { assert(i >= 0 && i < 3); return data[i]; }
//     float   operator[] (const int i) const { assert(i >= 0 && i < 3); return data[i]; }
//     vec& operator=(const vec &other) {
//         data[0] = other.data[0];
//         data[1] = other.data[1];
//         data[2] = other.data[2];
//         return *this;
//     }
//     vec  operator-() const { return {-x, -y, -z}; }
//     vec& operator+=(const vec &v) {
//         x += v.x;
//         y += v.y;
//         z += v.z;
//         return *this;
//     }
//     vec& operator*=(float t) {
//         x *= t;
//         y *= t;
//         z *= t;
//         return *this;
//     }
//     vec& operator/=(float t) { return *this *= 1/t; }
//     bool operator<(const vec<3,float> &v) const {
//         return ((x < v.x) && (y < v.y) && (z < v.z));
//     }

//     float norm2() const { return x*x + y*y + z*z; }
//     float norm()  const { return sqrt(norm2()); }
//     vec normalize() const { return (*this) * (1 / norm()); }
//     // vector compnent-wise product
//     vec mult(const vec<3,float> &v) const {
//         return {x*v.x, y*v.y, z*v.z};
//     }
//     // vector cross product
//     vec cross(const vec<3,float> &v) const {
//         return {y*v.z - z*v.y, z*v.x - x*v.z, x*v.y - y*v.x};
//     }
// };

// Specialization for vec<3,float>
template<> struct vec<3,float> {
    float x{}, y{}, z{};

    vec() = default;
    vec(float x_, float y_, float z_ = 0.0): x(x_), y(y_), z(z_) {}
    vec(const vec<2,float> &v, float z = 0.0): x(v.x), y(v.y), z(z) {}
    float & operator[] (const int i)       { assert(i >= 0 && i < 3); return i ? (i == 2 ? z : y) : x; }
    float   operator[] (const int i) const { assert(i >= 0 && i < 3); return i ? (i == 2 ? z : y) : x; }
    vec  operator-() const { return {-x, -y, -z}; }
    vec& operator+=(const vec &v) {
        x += v.x;
        y += v.y;
        z += v.z;
        return *this;
    }
    vec& operator*=(float t) {
        x *= t;
        y *= t;
        z *= t;
        return *this;
    }
    vec& operator/=(float t) { return *this *= 1/t; }
    bool operator<(const vec<3,float> &v) const {
        return ((x < v.x) && (y < v.y) && (z < v.z));
    }

    float norm2() const { return x*x + y*y + z*z; }
    float norm()  const { return sqrt(norm2()); }
    vec normalize() const { return (*this) * (1 / norm()); }
    // vector compnent-wise product
    vec mult(const vec<3,float> &v) const {
        return {x*v.x, y*v.y, z*v.z};
    }
    // vector cross product
    vec cross(const vec<3,float> &v) const {
        return {y*v.z - z*v.y, z*v.x - x*v.z, x*v.y - y*v.x};
    }
};

// overload operator +, vec3 + vec3
inline vec<3,float> operator+(const vec<3,float> &u, const vec<3,float> &v) {
    return {u.x + v.x, u.y + v.y, u.z + v.z};
}

// overload operator +, vec3 + float
inline vec<3,float> operator+(const vec<3,float> &u, const float &k) {
    return {u.x + k, u.y + k, u.z + k};
}

// overload operator -, vec3 - vec3
inline vec<3,float> operator-(const vec<3,float> &u, const vec<3,float> &v) {
    return {u.x - v.x, u.y - v.y, u.z - v.z};
}

// overload operator -, vec3 - float
inline vec<3,float> operator-(const vec<3,float> &u, const float &k) {
    return {u.x - k, u.y - k, u.z - k};
}

// overload operator *, float * vec3
inline vec<3,float> operator*(const float &k, const vec<3,float> &u) {
    return {u.x * k, u.y * k, u.z * k};
}

// overload operator *, vec3 * float addition
inline vec<3,float> operator*(const vec<3,float> &u, const float &k) {
    return k * u;
}

// overload operator *, vec3 dot vec3
inline float operator*(const vec<3,float> &u, const vec<3,float> &v) {
    return u.x * v.x + u.y * v.y + u.z * v.z;
}

// overload operator /, vec3 / float
inline vec<3,float> operator/(const vec<3,float> &v, const float &k) {
    return (1/k) * v;
}

// overload operator ==, vec3 == vec3
inline bool operator==(const vec<3,float> &u, const vec<3,float> &v) {
    return u.x == v.x && u.y == v.y && u.z == u.z;
}

// vector compnent-wise product
inline vec<3,float> mult(const vec<3,float> &v1, const vec<3,float> &v2) {
    return {v1.x*v2.x, v1.y*v2.y, v1.z*v2.z};
}

// vector cross product
inline vec<3,float> cross(const vec<3,float> &v1, const vec<3,float> &v2) {
    return {v1.y*v2.z - v1.z*v2.y, v1.z*v2.x - v1.x*v2.z, v1.x*v2.y - v1.y*v2.x};
}

// normalize the vector
inline vec<3,float> normalize(vec<3,float> &&v) {
    float k = 1.f / v.norm();
    v *= k;
    return std::forward<vec<3,float>>(v);
}

namespace std {
    template<> struct hash<vec<3,float>> {
        size_t operator()(const vec<3,float> &v) const {
            size_t seed = 0;
            hashCombine(seed, v.x);
            hashCombine(seed, v.y);
            hashCombine(seed, v.z);
            return seed;
        }
    };
}

// // Specialization for vec<4,float>
// template<> struct vec<4,float> {
//     float data[4] = {0};
//     float& x = data[0];
//     float& y = data[1];
//     float& z = data[2];
//     float& w = data[3];

//     vec() = default;
//     vec(float x_, float y_, float z_, float w_ = 0.0): data{x_, y_, z_, w_} {}
//     vec(const vec<3,float> &v, float w=0.0) {
//         data[0] = v.x;
//         data[1] = v.y;
//         data[2] = v.z;
//         data[3] = w;
//     }
//     vec(const vec& other) {
//         data[0] = other.data[0];
//         data[1] = other.data[1];
//         data[2] = other.data[2];
//         data[3] = other.data[3];
//     }
//     float & operator[] (const int i)       { assert(i >= 0 && i < 4); return data[i]; }
//     float   operator[] (const int i) const { assert(i >= 0 && i < 4); return data[i]; }

//     vec& operator=(const vec& other) {
//         data[0] = other.data[0];
//         data[1] = other.data[1];
//         data[2] = other.data[2];
//         data[3] = other.data[3];
//         return *this;
//     }
//     vec  operator-() const { return {-x, -y, -z, -w}; }
//     vec& operator+=(const vec &v) {
//         x += v.x;
//         y += v.y;
//         z += v.z;
//         w += v.w;
//         return *this;
//     }
//     vec& operator*=(float t) {
//         x *= t;
//         y *= t;
//         z *= t;
//         w *= t;
//         return *this;
//     }
//     vec& operator/=(float t) { return *this *= 1.f/t; }
//     bool operator<(const vec<4,float> &v) const {
//         return ((x < v.x) && (y < v.y) && (z < v.z)) && (w < v.w);
//     }

//     vec<3,float> xyz() const { return {x, y, z}; }
//     float norm2() const { return x*x + y*y + z*z + w*w; }
//     float norm()  const { return sqrt(norm2()); }
//     vec normalize() const { return (*this) * (1 / norm()); }
// };

// Specialization for vec<4,float>
template<> struct vec<4,float> {
    float x{}, y{}, z{}, w{};

    vec() = default;
    vec(float x_, float y_, float z_, float w_ = 0.0): x(x_), y(y_), z(z_), w(w_) {}
    vec(const vec<3,float> &v, float w=0.0): x(v.x), y(v.y), z(v.z), w(w) {}
    float & operator[] (const int i)       { assert(i >= 0 && i < 4); return i ? (i == 3 ? w : (i == 2 ? z : y)) : x; }
    float   operator[] (const int i) const { assert(i >= 0 && i < 4); return i ? (i == 3 ? w : (i == 2 ? z : y)) : x; }

    vec  operator-() const { return {-x, -y, -z, -w}; }
    vec& operator+=(const vec &v) {
        x += v.x;
        y += v.y;
        z += v.z;
        w += v.w;
        return *this;
    }
    vec& operator*=(float t) {
        x *= t;
        y *= t;
        z *= t;
        w *= t;
        return *this;
    }
    vec& operator/=(float t) { return *this *= 1.f/t; }
    bool operator<(const vec<4,float> &v) const {
        return ((x < v.x) && (y < v.y) && (z < v.z)) && (w < v.w);
    }

    vec<3,float> xyz() const { return {x, y, z}; }
    float norm2() const { return x*x + y*y + z*z + w*w; }
    float norm()  const { return sqrt(norm2()); }
    vec normalize() const { return (*this) * (1 / norm()); }
};

// overload operator +, vec3 + vec3
inline vec<4,float> operator+(const vec<4,float> &u, const vec<4,float> &v) {
    return {u.x + v.x, u.y + v.y, u.z + v.z, u.w + v.w};
}

// overload operator +, vec3 + float
inline vec<4,float> operator+(const vec<4,float> &u, const float &k) {
    return {u.x + k, u.y + k, u.z + k, u.w + k};
}

// overload operator -, vec3 - vec3
inline vec<4,float> operator-(const vec<4,float> &u, const vec<4,float> &v) {
    return {u.x - v.x, u.y - v.y, u.z - v.z, u.w - v.w};
}

// overload operator -, vec3 - float
inline vec<4,float> operator-(const vec<4,float> &u, const float &k) {
    return {u.x - k, u.y - k, u.z - k, u.w - k};
}

// overload operator *, float * vec3
inline vec<4,float> operator*(const float &k, const vec<4,float> &u) {
    return {u.x * k, u.y * k, u.z * k, u.w * k};
}

// overload operator *, vec3 * float addition
inline vec<4,float> operator*(const vec<4,float> &u, const float &k) {
    return k * u;
}

// overload operator *, vec3 dot vec3
inline float operator*(const vec<4,float> &u, const vec<4,float> &v) {
    return u.x * v.x + u.y * v.y + u.z * v.z + u.w * v.w;
}

// overload operator /, vec3 / float
inline vec<4,float> operator/(const vec<4,float> &v, const float &k) {
    return (1/k) * v;
}

inline bool operator==(const vec<4,float> &u, const vec<4,float> &v) {
    return u.x == v.x && u.y == v.y && u.z == u.z && u.w == v.w;
}

// vector compnent-wise product
inline vec<4,float> mult(const vec<4,float> &v1, const vec<4,float> &v2) {
    return {v1.x*v2.x, v1.y*v2.y, v1.z*v2.z, v1.w*v2.w};
}

// normalize the vector
inline vec<4,float> normalize(vec<4,float> &&v) {
    float k = 1.f / v.norm();
    v *= k;
    return std::forward<vec<4,float>>(v);
}

namespace std {
    template<> struct hash<vec<4,float>> {
        size_t operator()(const vec<4,float> &v) const {
            size_t seed = 0;
            hashCombine(seed, v.x);
            hashCombine(seed, v.y);
            hashCombine(seed, v.z);
            hashCombine(seed, v.w);
            return seed;
        }
    };
}

typedef vec<2,float> vec2;
typedef vec<3,float> vec3;
typedef vec<4,float> vec4;

/********************************************************************************
*                               Matrix Defination                               *
********************************************************************************/
template<int nrows, int ncols, typename T> struct mat;

template<typename T> T detRecursive(mat<2,2,T> m) {
    return m[0][0] * m[1][1] - m[0][1] * m[1][0];
}

template<int N, typename T> T detRecursive(mat<N,N,T> m) {
    mat<N-1,N-1,T> submatrix;
    vec<N-1,T> row;
    T det = 0.0;

    for (int j = 0; j < N; ++j) {
        for (int i = 1; i < N; ++i) {
            for (int k = 0; k < N; ++k) {
                if (k != j) {
                    row[i-1] = m[i][k];
                }
            }
            submatrix[i-1] = row;
        }

        det += m[0][j] * detRecursive(submatrix) * ((j % 2 == 0) ? 1 : -1);
    }

    return det;
}

template<int N, typename T> T detGauss(mat<N,N,T> m) {
    int swaps = 0;

    for (int i = 0; i < N - 1; ++i) {
        int pivot_row = i;
        for (int j = i + 1; j < N; ++j) {
            if (std::abs(m[j][i]) > std::abs(m[pivot_row][i])) {
                pivot_row = j;
            }
        }

        if (pivot_row != i) {
            std::swap(m[i], m[pivot_row]);
            swaps++;
        }

        for (int j = i + 1; j < N; ++j) {
            T factor = m[j][i] / m[i][i];
            for (int k = i; k < N; ++k) {
                m[j][k] -= factor * m[i][k];
            }
        }
    }

    T determinant = 1;
    for (int i = 0; i < N; ++i) {
        determinant *= m[i][i];
    }

    return (swaps % 2 == 0) ? determinant : -determinant;
}

template<int N, typename T> inline T det(mat<N,N,T> m) {
    if (N < 6)  return detRecursive(m);
    else        return detGauss(m);
}

// matrix define
template<int nrows, int ncols, typename T> struct mat {
    vec<ncols,T> rows[nrows] = {{}};

    vec<ncols,T>& operator[] (const int idx)       { assert(idx >= 0 && idx < nrows); return rows[idx]; }
    const vec<ncols,T>& operator[] (const int idx) const { assert(idx >= 0 && idx < nrows); return rows[idx]; }

    vec<nrows,T> col(const int idx) const {
        assert(idx >= 0 && idx < ncols);
        vec<nrows,T> res;
        for (int i = nrows; i--; res[i] = rows[i][idx]);
        return res;
    }

    mat<nrows,ncols,T>& setCol(const int idx, const vec<nrows,T>& v) {
        assert(idx >= 0 && idx < ncols);
        for (int i = nrows; i--; rows[i][idx] = v[i]);
        return *this;
    }

    // transpose the matrix
    mat<ncols,nrows,T> transpose() const {
        mat<ncols,nrows,T> res;
        for (int i = ncols; i--; res[i] = this->col(i));
        return res;
    }

    // get the complementary submatrix
    mat<nrows-1,ncols-1,T> getMinor(const int row, const int col) const {
        mat<nrows-1,ncols-1,T> res;
        for (int i = nrows - 1; i-- ;)
            for (int j = ncols - 1; j--; res[i][j] = rows[i < row ? i : i+1][j < col ? j : j+1]);
        return res;
    }

    // compute the cofactor
    T cofactor(const int row, const int col) const {
        return det(getMinor(row, col)) * ((row + col)%2 ? -1 : 1);
    }

    // compute the adjugate matrix
    mat<nrows,ncols,T> adjugate() const {
        mat<nrows,ncols,T> res;
        for (int i = nrows; i--; )
            for (int j = ncols; j--; res[j][i] = cofactor(i, j));
        return res;
    }

    // compute the transposed invert matrix
    mat<nrows,ncols,T> invertTranspose() const {
        return invert().transpose();
    }

    // compute the invert matrix using Gauss-Jordan elimination
    mat<nrows,ncols,T> invert() const {
        if (nrows != ncols) return zero();
        mat<nrows, nrows*2, T> augmentedMatrix;

        // Construct an augmented matrix [matrix | I]
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nrows; ++j) {
                augmentedMatrix[i][j] = rows[i][j];
                augmentedMatrix[i][j + nrows] = (i == j) ? 1 : 0;
            }
        }

        // Use Gauss-Jordan elimination to change the left part into an identity matrix
        for (int i = 0; i < nrows; ++i) {
            // Set the diagonal element to 1
            T pivot = augmentedMatrix[i][i];
            for (int j = 0; j < nrows * 2; ++j)
                augmentedMatrix[i][j] /= pivot;

            // elimination
            for (int k = 0; k < nrows; ++k) {
                if (k != i) {
                    T factor = augmentedMatrix[k][i];
                    for (int j = 0; j < nrows * 2; ++j)
                        augmentedMatrix[k][j] -= factor * augmentedMatrix[i][j];
                }
            }
        }

        // Extract the right part as the inverse matrix
        mat<nrows, nrows, T> inverseMatrix;
        for (int i = 0; i < nrows; ++i)
            for (int j = 0; j < nrows; ++j)
                inverseMatrix[i][j] = augmentedMatrix[i][j + nrows];

        return inverseMatrix;
    }

    // get zero matrix
    static mat<nrows,ncols,T> zero() {
        mat<nrows,ncols,T> res;
        for (int i = nrows; i--; )
            for (int j = ncols; j--; res[i][j] = 0);
        return res;
    }

    // get identity matrix
    static mat<nrows,ncols,T> identity() {
        mat<nrows,ncols,T> res;
        for (int i = nrows; i--; )
            for (int j = ncols; j--; res[i][j] = (i == j));
        return res;
    }
};

// override operator <<, print matrix
template<int nrows, int ncols, typename T>
std::ostream& operator<< (std::ostream& out, const mat<nrows,ncols,T>& m) {
    for (int i = 0; i < nrows; i++) std::cout << m[i] << std::endl;
    return out;
}

// override operator +, matrix - matrix addition
template<int nrows, int ncols, typename T> mat<nrows,ncols,T>
inline operator+ (const mat<nrows,ncols,T>& m, const mat<nrows,ncols,T>& n) {
    mat<nrows,ncols,T> res;
    for (int i = nrows; i--; res[i] = m[i] + n[i]) ;
    return res;
}

// override operator +, matrix - scalar addition
template<int nrows, int ncols, typename T> mat<nrows,ncols,T>
inline operator+ (const mat<nrows,ncols,T>& m, const T& k) {
    mat<nrows,ncols,T> res;
    for (int i = nrows; i--; res[i] = m[i] + k);
    return res;
}

// override operator -, matrix - matrix subtraction
template<int nrows, int ncols, typename T> mat<nrows,ncols,T>
inline operator- (const mat<nrows,ncols,T>& m, const mat<nrows,ncols,T>& n) {
    mat<nrows,ncols,T> res;
    for (int i = nrows; i--; res[i] = m[i] - n[i]);
    return res;
}

// override operator -, matrix - scalar addition
template<int nrows, int ncols, typename T> mat<nrows,ncols,T>
inline operator- (const mat<nrows,ncols,T>& m, const T& k) {
    mat<nrows,ncols,T> res;
    for (int i = nrows; i--; res[i] = m[i] - k) ;
    return res;
}

// override operator *, matrix - vector multiplication
template<int nrows, int ncols, typename T>
inline vec<nrows,T> operator* (const mat<nrows,ncols,T>& m, const vec<ncols,T>& v) {
    vec<nrows,T> res;
    for (int i = nrows; i-- ; res[i] = m[i] * v) ;
    return res;
}

// override operator *, matrix - matrix multiplication
template<int M, int N, int K, typename T>
inline mat<M,K,T>operator* (const mat<M,N,T>& lhs, const mat<N,K,T>& rhs) {
    mat<M,K,T> res;
    for (int i = M; i--; )
        for (int j = K; j--; res[i][j] = lhs[i] * rhs.col(j)) ;
    return res;
}

// override operator *, matrix - scalar multiplicat number
template<int nrows, int ncols, typename T>
inline mat<nrows,ncols,T> operator* (const T& k, const mat<nrows,ncols,T>& m) {
    mat<nrows,ncols,T> res;
    for (int i = nrows; i--; res[i] = m[i] * k) ;
    return res;
}

// override operator *, matrix - scalar multiplicat number
template<int nrows, int ncols, typename T>
inline mat<nrows,ncols,T> operator* (const mat<nrows,ncols,T>& m, const T& k) {
    return k * m;
}

// override operator /, matrix - scalar division
template<int nrows, int ncols, typename T>
inline mat<nrows,ncols,T> operator/ (const mat<nrows,ncols,T>& m, const T& k) {
    mat<nrows,ncols,T> res;
    for (int i = nrows; i--; res[i] = m[i] / k) ;
    return res;
}

// // Specialization for nrows == ncols == 3
// template<typename T> struct mat<3,3,T> {
//     vec<3,T> rows[3] = {{}};

//     vec<3,T>& operator[] (const int idx)       { assert(idx >= 0 && idx < 3); return rows[idx]; }
//     const vec<3,T>& operator[] (const int idx) const { assert(idx >= 0 && idx < 3); return rows[idx]; }

//     vec<3,T> col(const int idx) const {
//         assert(idx >= 0 && idx < 3);
//         return { rows[0][idx], rows[1][idx], rows[2][idx] };
//     }

//     void setCol(const int idx, const vec<3,T>& v) {
//         assert(idx >= 0 && idx < 3);
//         rows[0][idx] = v[0];
//         rows[1][idx] = v[1];
//         rows[2][idx] = v[2];
//     }

//     mat<3,3,T> transpose() const {
//         mat<3, 3, T> res;
//         res[0] = this->col(0);
//         res[1] = this->col(1);
//         res[2] = this->col(2);
//         return res;
//     }
// };

// // Specialization for nrows == ncols == 4
// template<typename T> struct mat<4,4,T> {
//     vec<4,T> rows[4] = {{}};

//     vec<4,T>& operator[] (const int idx)       { assert(idx >= 0 && idx < 4); return rows[idx]; }
//     const vec<4,T>& operator[] (const int idx) const { assert(idx >= 0 && idx < 4); return rows[idx]; }

//     vec<4,T> col(const int idx) const {
//         assert(idx >= 0 && idx < 4);
//         return { rows[0][idx], rows[1][idx], rows[2][idx], rows[3][idx] };
//     }

//     void setCol(const int idx, const vec<4,T>& v) {
//         assert(idx >= 0 && idx < 4);
//         rows[0][idx] = v[0];
//         rows[1][idx] = v[1];
//         rows[2][idx] = v[2];
//         rows[3][idx] = v[3];
//     }

//     mat<4,4,T> transpose() const {
//         mat<4,4,T> res;
//         res[0] = this->col(0);
//         res[1] = this->col(1);
//         res[2] = this->col(2);
//         res[3] = this->col(3);
//         return res;
//     }
// };

typedef mat<3,3,float> mat3;
typedef mat<4,4,float> mat4;

inline vec3 operator* (const mat3& m, const vec3& v) {
    return {m[0]*v, m[1]*v, m[2]*v};
}

inline vec4 operator* (const mat4& m, const vec4& v) {
    return {m[0]*v, m[1]*v, m[2]*v, m[3]*v};
}

// Const value defination
static const vec3 XA = {1, 0, 0};
static const vec3 YA = {0, 1, 0};
static const vec3 ZA = {0, 0, 1};
static const vec3 ONE_VEC3 = {1, 1, 1};
static const vec3 ZERO_VEC3 = {0, 0, 0};
static const mat3 E33 = mat3::identity();
static const mat4 E44 = mat4::identity();

/********************************************************************************
*                                Utility Methods                                *
********************************************************************************/
// Convert degree to radian
inline float deg2rad(float degrees) {
    return degrees * PI / 180.0;
}

// Clamp x to [min, max], default to [0, 1]
inline float clamp(float x, float min=0, float max=1) {
    return (x < min) ? min : ((x > max) ? max : x);
}

// Clamp vec2 to [min, max], default to [0, 1]
inline vec2 clamp(const vec2& v, float min=0, float max=1) {
    return {clamp(v.x, min, max), clamp(v.y, min, max)};
}

// Clamp vec3 to [min, max], default to [0, 1]
inline vec3 clamp(const vec3& v, float min=0, float max=1) {
    return {clamp(v.x, min, max), clamp(v.y, min, max), clamp(v.z, min, max)};
}

// Clamp vec4 to [min, max], default to [0, 1]
inline vec4 clamp(const vec4& v, float min=0, float max=1) {
    return {clamp(v.x, min, max), clamp(v.y, min, max), clamp(v.z, min, max), clamp(v.w, min, max)};
}

// Clamp vecn to [min, max], default to [0, 1]
template <int n, typename T>
inline vec<n,T> clamp(const vec<n,T>& v, float min=0, float max=1) {
    vec<n,T> res;
    for (int i = 0; i < n; -- i)
        res[i] = clamp(v[i], min, max);
    return res;
}

// Clamp vec3 to [min, max], default to [0, 1]
inline vec3 clampVec3(const vec3& v, float min=0, float max=1) {
    return {clamp(v.x, min, max), clamp(v.y, min, max), clamp(v.z, min, max)};
}

// Gamma correction with gamma level
inline void gammaCorrection(vec3 &color, float gamma = 2.2) {
    float gamma_recip = 1.0 / gamma;
    color.x = pow(color.x, gamma_recip);
    color.y = pow(color.y, gamma_recip);
    color.z = pow(color.z, gamma_recip);
}

// Generate a random real in [0,1].
inline float randDouble() {
    static std::uniform_real_distribution<float> distribution(0.0, 1.0);
    static std::mt19937 generator;
    return distribution(generator);
}

// Generate a random real in [min,max].
inline float randDouble(float min, float max) {
    return min + (max-min)*randDouble();
}

// Generate a random int in [min,max].
inline float randInt(float min, float max) {
    return static_cast<int>(randDouble(min, max+1));
}

// Sample in [0, 1]^2
static vec2 sampleUnitSpuare() {
    float x = randDouble();
    float y = randDouble();
    return vec2(x, y);
}

// Sample in uniu disk
static vec2 sampleUnitDisk() {
    float r = randDouble();
    float phi = 2*PI*randDouble();
    return vec2(r*cos(phi), r*sin(phi));
}

static vec3 randVec3() {
    return vec3(randDouble(), randDouble(), randDouble());
}

static vec3 randVec3(float min, float max) {
    return vec3(randDouble(min, max), randDouble(min, max), randDouble(min, max));
}

// Generate random vector in sphere
static vec3 randVecSphere() {
    // Assume that n is normalized
    // float phi = 2*PI*randDouble(), r = randDouble(-1, 1);
    // float a, b;
    // if (r > 0) { a = sqrt(1 - r), b = sqrt(r); }
    // else {a = sqrt(1 + r), b = sqrt(-r); }
    // return vec3(cos(phi)*a, sin(phi)*a, b).normalize();
    float phi = 2*PI*randDouble(), theta = PI*randDouble(-1, 1);
    float a = cos(theta), b = sin(theta);
    return vec3(cos(phi)*a, sin(phi)*a, b).normalize();
}

// Generate random vector in hemisphere
static vec3 randVecSemisphere(const vec3& n) {
    // Assume that n is normalized
    float phi = 2*PI*randDouble(), r = randDouble();
    float a = sqrt(r), b = sqrt(1 - r);
    vec3 u = YA.cross(n);
    u = (u.norm() > EPS ? u : XA.cross(n)).normalize();
    // u = (u.norm() > EPS ? u : XA.cross(n));
    vec3 v = n.cross(u);
    return (u*cos(phi)*a + v*sin(phi)*a + n*b).normalize();
}

// Compute Reflect Ray
inline vec3 reflect(const vec3& v, const vec3& n) {
    return (v - 2*(v*n)*n).normalize();
}

inline vec3 refract(const vec3& rin, const vec3& n, float etai_over_etat) {
    float cos_theta1 = -rin * n;
    float cos2_theta2 = 1 - etai_over_etat*etai_over_etat*(1 - cos_theta1*cos_theta1);
    // Use Schlick's approximation for reflectance.
    float r0 = (1 - etai_over_etat) / (1 + etai_over_etat);
    r0 = r0 * r0;
    float rate = r0 + (1 - r0)*pow(1 - cos_theta1, 5);
    if (cos2_theta2 < 0 || rate > randDouble()) {
        // std::cout << "Total reflection" << std::endl;
        return reflect(rin, n);    // Total reflection
    }
    vec3 rout_perp =  etai_over_etat * (rin + cos_theta1*n);
    vec3 rout_para = -n * sqrt(cos2_theta2);
    return (rout_perp + rout_para).normalize();
}

inline mat4 scale(const vec3 &s) {
    mat4 S = E44;
    for (int i = 0; i < 3; i ++) S[i][i] = s[i];
    return S;
}

inline mat4 translate(const vec3 &t) {
    mat4 T = E44;
    for (int i = 0; i < 3; i ++) T[i][3] = t[i];
    return T;
}

inline mat4 rotateX(float degree) {
    float radian = deg2rad(degree);
    mat4 R;
    R = {{{1.f, 0.f, 0.f, 0.f},
                {0.f, cos(radian), -sin(radian), 0.f},
                {0.f, sin(radian), cos(radian), 0.f},
                {0.f, 0.f, 0.f, 1.f}}};
    return R;
}

inline mat4 rotateY(float degree) {
    float radian = deg2rad(degree);
    mat4 R;
    R = {{{cos(radian), 0.f, sin(radian), 0.f},
                {0.f, 1.f, 0.f, 0.f},
                {-sin(radian), 0.f, cos(radian), 0.f},
                {0.f, 0.f, 0.f, 1.f}}};
    return R;
}

inline mat4 rotateZ(float degree) {
    float radian = deg2rad(degree);
    mat4 R;
    R = {{{ cos(radian), -sin(radian), 0.f, 0.f},
                { sin(radian),  cos(radian), 0.f, 0.f},
                { 0.f, 0.f, 1.f, 0.f},
                { 0.f, 0.f, 0.f, 1.f}}};
    return R;
}

inline mat4 rotate(const vec3 &n, float degree) {
    mat4 R = E44;
    return R;
}


// View transform
inline mat4 lookAt(const vec3 &position, const vec3 &direction, const vec3 &up) {
    vec3 z = -direction.normalize();
    vec3 x = cross(up, z).normalize();
    vec3 y = cross(z, x).normalize();
    mat4 Rview = E44;
    mat4 Tview = E44;
    Tview[0][3] = -position.x;
    Tview[1][3] = -position.y;
    Tview[2][3] = -position.z;
    Rview[0] = vec4(x);
    Rview[1] = vec4(y);
    Rview[2] = vec4(z);
    return Rview * Tview;
}

// Projection transform
inline mat4 projection(float fov, float aspectRatio, float near, float far) {
    // fear and near are always grateter than 0
    float h = 2*near*tan(deg2rad(fov / 2));
    float w = aspectRatio * h;
    float n = near;
    float f = far;
    mat4 ortho, persp2ortho;
    ortho = {{{2.f/w, 0.f, 0.f, 0.f},
                    {0.f, 2.f/h, 0.f, 0.f},
                    {0.f, 0.f, 2.f/(n-f), (n+f)/(n-f)},
                    {0.f, 0.f, 0.f, 1.f}}};
    persp2ortho = {{{n, 0.f, 0.f, 0.f},
                        {0.f, n, 0.f, 0.f},
                        {0.f, 0.f, n+f, n*f},
                        {0.f, 0.f, -1.f, 0.f}}};
    return ortho * persp2ortho;
}

// Viewport transform
inline mat4 viewport(int screen_w, int screen_h) {
    mat4 m = E44;
    m[0][0] = screen_w / 2.f;
    m[1][1] = screen_h / 2.f;
    m[2][2] = 0.5f;

    m[0][3] = screen_w / 2.f;
    m[1][3] = screen_h / 2.f;
    m[2][3] = 0.5f;

    return m;
}

// Linear Interpolation
inline float mylerp(float a, float b, float t) {
    return a + t * (b - a);
}

// Bilinear Interpolation
inline float bilerp(float q11, float q12, float q21, float q22, float s, float t) {
    float r1 = mylerp(q11, q21, s);
    float r2 = mylerp(q12, q22, s);
    return mylerp(r1, r2, t);
}

// Interpolation with perspective correction
inline float interpLine(float a1, float a2, float z1, float z2, float s) {
    // z1, z2 is linear depth, s is interpolation factor
    // a1, a2 is attribute to be interpolated
    return ((1-s)*a1*z1 + s*a2*z1) / (s*z1 + (1-s)*z2);
}

// Depth Interpolation with perspective correction
inline float interpZ(const vec3 &bc, float depth[3]) {
    // x, y is the corrdination in screen space
    // bc is the barycentric coordination of point (x, y) 
    // depth is the depth of vertics in view space
    // 1/Z = alpha/Za + beta/Zb + gamma/Zc
    float z_interp = 1.0 / (bc.x / depth[0] + bc.y / depth[1] + bc.z / depth[2]);
    return z_interp;
}

// // Attritube Interpolation with perspective correction
// inline float interp(const vec3 &I, const vec3 &bc, float depth[3], float z_interp) {
//     // x, y is the corrdination in screen space
//     // bc is the barycentric coordination of point (x, y) 
//     // depth is the depth of vertics in view space
//     // I is the attritube to be interpolated
//     float I_interp = bc.x * I[0] / depth[0] + bc.y * I[1] / depth[1] + bc.z * I[2] / depth[2];
//     I_interp *= z_interp;
//     return I_interp;
// }

// Attritube Interpolation with perspective correction
template<typename T>
inline T interp(const T &i1, const T &i2, const T &i3, float alpha, float beta, float gamma) {
    return i1 * alpha + i2 * beta + i3 * gamma;
}

// Compute barycentric coordinate in screen space
static std::tuple<float, float, float> computeBarycentric(float x, float y, const std::array<vec4,3> &v) {
    float c1 = (x*(v[1].y - v[2].y) + y*(v[2].x - v[1].x) + v[1].x*v[2].y - v[2].x*v[1].y) / (v[0].x*(v[1].y - v[2].y) + v[0].y*(v[2].x - v[1].x) + v[1].x*v[2].y - v[2].x*v[1].y);
    float c2 = (x*(v[2].y - v[0].y) + y*(v[0].x - v[2].x) + v[2].x*v[0].y - v[0].x*v[2].y) / (v[1].x*(v[2].y - v[0].y) + v[1].y*(v[0].x - v[2].x) + v[2].x*v[0].y - v[0].x*v[2].y);
    float c3 = (x*(v[0].y - v[1].y) + y*(v[1].x - v[0].x) + v[0].x*v[1].y - v[1].x*v[0].y) / (v[2].x*(v[0].y - v[1].y) + v[2].y*(v[1].x - v[0].x) + v[0].x*v[1].y - v[1].x*v[0].y);
    return {c1, c2, c3};
}

#endif