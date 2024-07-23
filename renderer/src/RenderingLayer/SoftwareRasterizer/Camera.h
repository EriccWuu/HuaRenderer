#ifndef __CAMERA_H__
#define __CAMERA_H__

#include "CoreLayer/MathLib.h"

class Camera {
public:
    vec4 position;
    vec4 direction;
    vec4 up;
    vec4 right;
    double fov;
    double aspectRatio;
    double width, height;
    double near, far;

    Camera();
    Camera(const vec3 &pos, const vec3 &dir={0,0,-1}, const vec3 &up={0,1,0}, double w=100, double h=100, double fov=45, double r=1, double n=-0.5, double f=-50);

private:
    void init(const vec3 &pos, const vec3 &dir, const vec3 &up, double w=100, double h=100, double fov=45, double r=1, double n=-0.5, double f=-50);
};

#endif