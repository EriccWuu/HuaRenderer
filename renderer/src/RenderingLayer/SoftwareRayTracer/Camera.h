#ifndef __CAMERA_H__
#define __CAMERA_H__

#include "CoreLayer/MathLib.h"
#include "CoreLayer/TGAImage.h"
#include "Utility.h"

class Camera {
public:
    vec4 position;
    vec4 direction;
    vec4 up;
    vec4 right;
    double fov;
    double aspectRatio;
    double near, far;
    double defocusAngle = 0;
    int screen_w, screen_h;
    int spp = 100, maxDepth = 5;
    int SSAA = 1;
    vec3 backgroundLight = ZERO_VEC3;

    TGAImage frame;
    TGAImage zbuffer;

    Camera() = default;
    Camera(vec3 pos, vec3 dir, vec3 u, int h=100, double r=1, double fov=60, double n=-0.5, double f=-50);

    void set(vec3 pos, vec3 dir, vec3 u, double h=100, double r=1, double fov=60, double n=-0.5, double f=-50);
    mat4 view();
    mat4 projection();
    mat4 viewport();
    void init();
    void render(const BVHNode &objects, TGAImage &image);

private:
    vec3 viewO;
    vec3 viewportO, viewportU, viewportV;
    vec3 defocusDiskU, defocusDiskV; 

    inline double width();
    inline double height();
    Ray  getRay(const int i, const int j, const int sx, const int sy);
    vec3 radiance(const Ray &r, int depth, const BVHNode &obj);
    inline vec2 pixSampleSquare();
    inline vec2 pixSampleDisk(double radius = 1.0);
};

#endif