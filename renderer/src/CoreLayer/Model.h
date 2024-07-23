#ifndef __MODEL_H__
#define __MODEL_H__

#include <vector>
#include <string>
#include "MathLib.h"
#include "TGAImage.h"
#include "Vertex.h"
#include "Texture.h"
#include "Triangle.h"
#include <unordered_map>
#include <memory>

struct Index {
    int vertex_index;
    int texcoord_index;
    int normal_index;
};

struct Mesh {
    std::vector<Index> indices;
};

struct Attribute {
    std::vector<vec3> vertices;
    std::vector<vec3> normals;
    std::vector<vec2> texcoords;
    std::vector<vec3> colors;

    Attribute() {}
};

class Model {
    // Texture texDiffuse;
    // Texture texSpecular;
    // Texture texNormal;
    std::shared_ptr<std::vector<Vertex>> uniqueVerticesPtr;
    std::shared_ptr<std::vector<size_t>> uniqueIndicesPtr;

    std::string texture_path(const std::string &filename, const std::string &suffix);
    void load_texture(const std::string filename, const std::string suffix, TGAImage &img);
    void makeVerticesUnique();

public:
    Attribute attribute;
    Mesh mesh;
    std::shared_ptr<Texture> diffuseTexPtr;
    std::shared_ptr<Texture> specularTexPtr;
    std::shared_ptr<Texture> normalTexPtr;

    Model(const std::string filename);
    std::shared_ptr<std::vector<Vertex>> uniqueVertices();
    std::shared_ptr<std::vector<size_t>> uniqueIndices();
    std::shared_ptr<Texture> diffuseTexture();
    std::shared_ptr<Texture> specularTexture();
    std::shared_ptr<Texture> normalTexture();
    int nverts() const;
    int nfaces() const;
    vec3 normal(const int iface, const int nthvert) const;
    vec3 vert(const int i) const;
    vec3 vert(const int iface, const int nthvert) const;
    Index face(const int i) const;
    vec2 uv(const int iface, const int nthvert) const;
    vec4 diffuse(vec2 uv) const;
    vec4 specular(vec2 uv) const;
    vec4 normal(const vec2 &uv) const;
};

#endif