#ifndef __TRIANGLE_H__
#define __TRIANGLE_H__

#include "Vertex.h"
#include <array>

struct Triangle {
    std::array<size_t,3> indices;

    Triangle(size_t v0, size_t v1, size_t v2): indices({v0, v1, v2}) {}
    Triangle(std::array<size_t, 3> indices_): indices(indices_) {}
};

#endif