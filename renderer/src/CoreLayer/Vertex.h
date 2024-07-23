#ifndef __VERTEX_H__
#define __VERTEX_H__

#include "MathLib.h"

struct Vertex {
    vec4 pos;
    vec4 normal;
    vec4 tangent;
    vec4 bitangent;
    vec3 color;
    vec2 texcoord;

    bool operator==(const Vertex &other) const {
        return pos == other.pos && texcoord == other.texcoord && normal == other.normal && color == other.color;
    }
};

namespace std {
    template<> struct hash<Vertex> {
        size_t operator()(const Vertex& v) const {
            size_t hashValue = 0;
            hashCombine(hashValue, std::hash<vec4>()(v.pos));
            hashCombine(hashValue, std::hash<vec4>()(v.normal));
            hashCombine(hashValue, std::hash<vec3>()(v.color));
            hashCombine(hashValue, std::hash<vec2>()(v.texcoord));
            return hashValue;
        }
    };
};

#endif