#ifndef __TEXTURE_H__
#define __TEXTURE_H__

#include "TGAImage.h"

struct Texture {
    unsigned int id;
    std::string type;
    std::string path;

    Texture();
    Texture(const TGAImage &textureImage, const std::string type = "undifined");
    Texture(std::vector<vec4> &textureData, const int width, const int height, const std::string type = "undifined");
    Texture(const std::string path, const std::string type = "undifined");
    Texture(const int width, const int height, const std::string type = "undifined");

    void swap(std::vector<vec4> &textureData, const int width, const int height);
    void set(int x, int y, const vec4 &c);
    vec4 get(int x, int y);
    vec4 value(const vec2 &uv) const;
    size_t width();
    size_t height();

private:
    int w;
    int h;
    std::vector<vec4> data;

    void load_texture(const std::string path);
    void load_from_image(const TGAImage &image);
    void load_from_vector(const std::vector<vec4> &vec);
};

#endif