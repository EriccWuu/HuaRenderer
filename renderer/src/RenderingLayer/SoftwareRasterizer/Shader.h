#ifndef __SHADER_H__
#define __SHADER_H__

#include "CoreLayer/MathLib.h"
#include "CoreLayer/Vertex.h"
#include "CoreLayer/Texture.h"
#include <unordered_map>
#include <variant>
#include <memory>

using ShaderContext = std::unordered_map<int, std::variant<int, float, vec2, vec3, vec4>>;

class IShader {
public:
    ShaderContext context;
    mat<2,3,float> varying_uv;
    mat<3,3,float> varying_norm;
    mat<3,3,float> varying_tang;
    mat<3,3,float> varying_bitang;
    mat<3,3,float> varying_pos;

    // Basic Shader uniform resources
    mat4 uniform_MODEL;
    mat4 uniform_VIEW;
    mat4 uniform_PROJECTION;
    mat4 uniform_VIEWPORT;
    mat4 uniform_NORMAL;

    IShader() {};
    virtual void vertex(const Vertex &inVertex, ShaderContext &vertexOut) = 0;
    virtual bool fragment(vec3 &gl_Fragcolor) = 0;
};

class Shader: public IShader {
public:
    struct Binding {
        const int pos = 0;
        const int uv = 1;
        const int norm = 2;
        const int tang = 3;
        const int bitang = 4;
        const int fragPos = 5;
    } binding;

    // Shader uniform resources
    vec3 uniform_lightDir;
    vec3 uniform_lightPos;
    mat4 uniform_Mshadow;
    std::shared_ptr<Texture> uniform_shadowmapPtr;
    std::shared_ptr<Texture> uniform_diffuseTexPtr;
    std::shared_ptr<Texture> uniform_specularTexPtr;
    std::shared_ptr<Texture> uniform_normalTexPtr;

    Shader() {}
    void vertex(const Vertex &inVertex, ShaderContext &vertexOut) override;
    bool fragment(vec3 &gl_Fragcolor) override;
};

class DepthShader: public IShader {
public:
    struct Binding {
        const int pos = 0;
    } binding;
    vec3 varying_depth;

    DepthShader() {}
    void vertex(const Vertex &inVertex, ShaderContext &vertexOut) override;
    bool fragment(vec3 &gl_Fragcolor) override;
};

#endif