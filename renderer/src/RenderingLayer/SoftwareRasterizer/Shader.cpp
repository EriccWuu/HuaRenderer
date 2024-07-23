#include "Shader.h"
#include <algorithm>

/***********************************************************************************
********                                Shader                              ********
***********************************************************************************/
void Shader::vertex(const Vertex &inVertex, ShaderContext &vertexOut) {
    vec3 n = (uniform_NORMAL * inVertex.normal).xyz();
    vec3 t = (uniform_NORMAL * inVertex.tangent).xyz();
    t = normalize(t - (t*n)*n);
    vec3 b = cross(n, t);
    vec4 fragPos = uniform_MODEL * inVertex.pos;
    vec4 viewPos = uniform_VIEW * fragPos;

    vertexOut[binding.fragPos] = fragPos;
    vertexOut[binding.norm] = n;
    vertexOut[binding.tang] = t;
    vertexOut[binding.bitang] = b;
    vertexOut[binding.uv] = inVertex.texcoord;
    vertexOut[binding.pos] = uniform_PROJECTION * viewPos;
}

bool Shader::fragment(vec3 &gl_Fragcolor) {
    vec2 uv = std::get<vec2>(context[binding.uv]);
    vec3 fragPos = std::get<vec4>(context[binding.fragPos]).xyz();
    vec3 bn = std::get<vec3>(context[binding.norm]).normalize();
    vec3 t = std::get<vec3>(context[binding.tang]).normalize();
    vec3 b = std::get<vec3>(context[binding.bitang]).normalize();

    vec3 n, l, r;
    mat3 TBN;
    TBN.setCol(0, t).setCol(1, b).setCol(2, bn);
    n = normalize(TBN * (2 * uniform_normalTexPtr->value(uv).xyz() - 1.f));

    l = -uniform_lightDir;
    r = normalize(n*(n*l*2.f) - l);

    // Texture Map + Normal Map (Tangent Space) + bling-phong + shadow map:
    vec4 currPixPos = uniform_Mshadow * vec4(fragPos, 1.0);  // Now in light space
    float rhw = 1 / currPixPos.w;
    currPixPos = currPixPos * rhw;
    // float depth = 0.5 + 0.5 * currPixPos.z;
    float depth = rhw;
    float shadow;
    float u = (currPixPos.x + 1) / 2;
    float v = (currPixPos.y + 1) / 2;
    if (u >= 1.f || v >= 1.f || u < 0 || v < 0)
        shadow = 1;
    else {
        float depthval = uniform_shadowmapPtr->value({u, v}).x;
        float dotProduct = n * normalize(uniform_lightPos - fragPos);
        float bias = 0.05 * std::max(0.f, 1.f - dotProduct);
        shadow = (depth + bias< depthval) ? 0.5f : 1.f;
    }

    auto c = uniform_diffuseTexPtr->value(uv).xyz() * 255;
    float ambi = 0.3; // Ambient term
    float diff = std::max(n*l, 0.f);    // diffuse term
    diff = 0.6;
    float f = uniform_specularTexPtr->value(uv).z * 255;
    float spec;
    if (f < 1e-6) spec = 0;
    else spec = 0.5*pow(std::max(r.z, 0.f), f);
    float intensity = diff + spec;

    gl_Fragcolor = clamp(c * (ambi + intensity*shadow), 0, 255);

    return false;
}

/***********************************************************************************
********                          Depth Shader                              ********
***********************************************************************************/
void DepthShader::vertex(const Vertex &inVertex, ShaderContext &vertexOut) {
    vec4 fragPos = uniform_MODEL * inVertex.pos;
    vec4 viewPos = uniform_VIEW * fragPos;
    vec4 gl_Position = uniform_PROJECTION * viewPos;

    vertexOut[binding.pos] = gl_Position;
}

bool DepthShader::fragment(vec3 &gl_Fragcolor) {
    float depth = std::get<vec4>(context[binding.pos]).w;
    gl_Fragcolor = ONE_VEC3 * (1 / depth);
    return false;
}
