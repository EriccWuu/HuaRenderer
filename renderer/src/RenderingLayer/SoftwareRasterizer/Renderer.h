#ifndef __RENDERER_H__
#define __RENDERER_H__

#include "CoreLayer/MathLib.h"
#include "Camera.h"
#include "Shader.h"
#include <memory>

class Renderer {
private:
    enum ReginCode {
        eCenter = 0b000000,
        eLeft = 0b000001,
        eRight = 0b000010,
        eButtom = 0b000100,
        eTop = 0b001000,
        eFront = 0b010000,
        eBack = 0b100000,
    };

    int width_, height_;
    std::shared_ptr<std::vector<Vertex>> vertexBufferPtr_;
    std::shared_ptr<std::vector<size_t>> indexBufferPtr_;
    std::vector<ShaderContext> shaderContextBuffer_;
    std::vector<float> zbuffer_;
    mat4 viewport_;
    bool drawMesh_ = false;

    // Find bounding box of triangle
    std::vector<vec3> bounding_box(const std::array<vec4,3> &vertices);

    // draw line
    void drawline(int x1, int y1, int x2, int y2, std::vector<vec4> &image, const vec3 &color);

    // draw triangle
    void drawMesh(const std::array<size_t,3> &indices, std::vector<vec4> &image, const vec3 &color);
    void drawPrimitive(const std::array<size_t,3> &indices, IShader &s, std::vector<vec4> &image);

    int encode(float x, float y);
    bool lineClip(float ix1, float iy1, float ix2, float iy2, float &ox1, float &oy1, float &ox2, float &oy2);

    void interpolation(const std::array<size_t,3> &idx, const vec3 &bar, ShaderContext &interpolated);

public:
    Renderer(int width, int height);
    ~Renderer();

    void setViewport(const mat4& viewport);
    void initZbuffer();
    void bindVertexBuffer(std::shared_ptr<std::vector<Vertex>> ptr);
    void bindIndexBuffer(std::shared_ptr<std::vector<size_t>> ptr);
    void draw( IShader &shader, std::vector<vec4> &image);
    void meshMode(bool drawMesh);

};

#endif