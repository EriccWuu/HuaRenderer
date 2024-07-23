#include "Texture.h"

Texture::Texture() {
    this->path = "None";
    this->type = "None";
}

Texture::Texture(const TGAImage &textureImage, const std::string type): type(type) {
    load_from_image(textureImage);
};

Texture::Texture(std::vector<vec4> &textureData, int width, int height, const std::string type) {
    this->type = type;
    w = width;
    h = height;
    load_from_vector(textureData);
}

Texture::Texture(const std::string path, const std::string type) {
    this->path = path;
    this->type = type;
    load_texture(this->path);
}

Texture::Texture(const int width, const int height, const std::string type) {
    data = std::vector<vec4>(width*height);
    w = width;
    h = height;
    this->type = type;
}

void Texture::swap(std::vector<vec4> &textureData, const int width, const int height) {
    w = width;
    h = height;
    data.swap(textureData);
}

void Texture::set(int x, int y, const vec4 &c) {
    data[y*w + x] = c;
}

vec4 Texture::get(int x, int y) {
    return data[x + y * w];
}

size_t Texture::width() { return w; }

size_t Texture::height() { return h; }

vec4 Texture::value(const vec2 &uv) const {
    if (data.size() == 0) return vec4{0.f, 0.f, 0.f, 0.f};
    float u = uv.x;
    float v = uv.y;
    u = u * (w - 1) + 0.5f;
    v = v * (h - 1) + 0.5f;
    if (u >= w || u < 0 || v >= h || v < 0) return ZERO_VEC3;
    unsigned int u0 = uint32_t(u);
    unsigned int v0 = uint32_t(v);
    unsigned int u1 = std::min(u0 + 1, (unsigned int)w - 1);
    unsigned int v1 = std::min(v0 + 1, (unsigned int)h - 1);
    float s = u - u0, t = v - v0;
    vec4 c00 = data[u0 + v0 * w];
    vec4 c01 = data[u0 + v1 * w];
    vec4 c10 = data[u1 + v0 * w];
    vec4 c11 = data[u1 + v1 * w];
    float r = bilerp(c00.x, c01.x, c10.x, c11.x, s, t);
    float g = bilerp(c00.y, c01.y, c10.y, c11.y, s, t);
    float b = bilerp(c00.z, c01.z, c10.z, c11.z, s, t);
    float a = bilerp(c00.w, c01.w, c10.w, c11.w, s, t);

    return {r, g, b, a};
}

void Texture::load_texture(const std::string path) {
    TGAImage image;  
    std::cerr << "Loading texture file \"" << path << "\"...\n";
    std::cerr << (image.read_tga_file(path.c_str()) ? "- Load successfully" : "- Fail to load.") << std::endl;
    load_from_image(image);
}

void Texture::load_from_image(const TGAImage &image) {
    w = image.width();
    h = image.height();
    data = std::vector<vec4>(w*h);
    for (int i = 0; i < w; ++ i) {
        for (int j = 0; j < h; ++ j) {
            auto c = image.get(i, j);
            data[j*w + i].x = c.bgra[2] / 255.f;
            data[j*w + i].y = c.bgra[1] / 255.f;
            data[j*w + i].z = c.bgra[0] / 255.f;
            data[j*w + i].w = c.bgra[3] / 255.f;
        }
    }
}

void Texture::load_from_vector(const std::vector<vec4> &vec) {
    data = std::vector<vec4>(w*h, -1 * ONE_VEC3);
    for (int i = 0; i < w; ++ i) {
        for (int j = 0; j < h; ++ j) {
            data[j*w + i] = vec[j*w + i];
        }
    }
}