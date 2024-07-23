#include "Camera.h"

Camera::Camera() {
    position = vec4(0.f, 0.f, 0.f, 1.f);
    direction = vec4(0.f, 0.f, -1.f, 0.f);
    up = vec4(0.f, 1.f, 0.f, 0.f);
    right = vec4(cross(direction.xyz(), up.xyz()).normalize(), 0.f);
    up = vec4(cross(right.xyz(), direction.xyz()).normalize(), 0.f);
    fov = 45;
    aspectRatio = 1;
    near = -0.5;
    far = -50;
    height = 100;
    width = 100;
}

Camera::Camera(const vec3 &pos, const vec3 &dir, const vec3 &up, double w, double h, double fov, double r, double n, double f) {
    init(pos, dir, up, w, h, fov, r, n, f);
}

void Camera::init(const vec3 &pos, const vec3 &dir, const vec3 &up, double w, double h, double fov, double r, double n, double f) {
    position = vec4(pos, 1.f);
    direction = vec4(dir.normalize(), 0.f);
    right = vec4(cross(direction.xyz(), up).normalize(), 0.f);
    this->up = vec4(cross(right.xyz(), direction.xyz()).normalize(), 0.f);
    this->fov = fov;
    aspectRatio = r;
    near = n;
    far = f;
    height = w;
    width = h;
}