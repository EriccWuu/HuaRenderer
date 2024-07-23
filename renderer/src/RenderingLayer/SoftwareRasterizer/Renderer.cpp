#include "Renderer.h"

Renderer::Renderer(int width, int height): width_(width), height_(height) {
    initZbuffer();
}

Renderer::~Renderer() {
}

void Renderer::setViewport(const mat4& viewport) {
    viewport_ = viewport;
}

void Renderer::meshMode(bool drawMesh) {
    drawMesh_ = drawMesh;
}

void Renderer::initZbuffer() {
    int size = width_ * height_;
    zbuffer_ = std::vector<float>(size, 0);
}

void Renderer::bindVertexBuffer(std::shared_ptr<std::vector<Vertex>> ptr) {
    vertexBufferPtr_ = ptr;
    shaderContextBuffer_.clear();
    shaderContextBuffer_.resize(vertexBufferPtr_->size());
}

void Renderer::bindIndexBuffer(std::shared_ptr<std::vector<size_t>> ptr) {
    indexBufferPtr_ = ptr;
}

void Renderer::draw(IShader &shader, std::vector<vec4> &image) {
    if (vertexBufferPtr_ == nullptr || indexBufferPtr_ == nullptr) {
        std::cerr << "No Model can be draw." << '\n';
        return ;
    }

    // Vertex shader
    int numVert = vertexBufferPtr_->size();
    for (int i = 0; i < numVert; ++ i) {
        shader.vertex(vertexBufferPtr_->at(i), shaderContextBuffer_[i]);
    }

    int nface = indexBufferPtr_->size() / 3;
    for (int i = 0; i < nface; i ++) {
        std::array<size_t, 3> idx = {
            indexBufferPtr_->at(3*i), 
            indexBufferPtr_->at(3*i + 1), 
            indexBufferPtr_->at(3*i + 2)
        };

        // Render
        if (drawMesh_)
            drawMesh(idx, image, {255, 255, 255});
        else
            drawPrimitive(idx, shader, image);
    }
}

// Find bounding box of triangle
std::vector<vec3> Renderer::bounding_box(const std::array<vec4,3> &vertices) {
    auto v0 = vertices[0];
    auto v1 = vertices[1];
    auto v2 = vertices[2];
    vec3 bottomleft, upperright;
    bottomleft.x = std::floor(std::min(v0.x, std::min(v1.x, v2.x)));
    bottomleft.y = std::floor(std::min(v0.y, std::min(v1.y, v2.y)));
    bottomleft.z = std::floor(std::min(v0.z, std::min(v1.z, v2.z)));
    upperright.x = std::ceil(std::max(v0.x, std::max(v1.x, v2.x)));
    upperright.y = std::ceil(std::max(v0.y, std::max(v1.y, v2.y)));
    upperright.z = std::ceil(std::max(v0.z, std::max(v1.z, v2.z)));
    std::vector<vec3> bbox = {bottomleft, upperright};
    return bbox;
}

// draw line
void Renderer::drawline(int x1, int y1, int x2, int y2, std::vector<vec4> &image, const vec3 &color) {
    int dx = x2 - x1, dy = y2 - y1;
    bool isSteep = abs(dx) < abs(dy);
    if (isSteep) {
        std::swap(x1, y1);
        std::swap(x2, y2);
    }

    if (x1 > x2) {
        std::swap(x1, x2);
        std::swap(y1, y2);
    }
    dx = x2 - x1, dy = y2 - y1;
    bool isNegSlope = dy < 0;

    int eps = 0;
    for (int x = x1, y = y1; x <= x2; x ++) {

        if (!isSteep)
            image[y * width_ + x] = color;
        else
            image[x * width_ + y] = color;

        eps += dy;
        if (!isNegSlope) {
            if ((eps << 1) >= dx) {
                y ++;
                eps -= dx;
            }
        }
        else {
            if ((eps << 1) <= -dx) {
                y --;
                eps += dx;
            }
        }
    }
}

int Renderer::encode(float x, float y) {
    int l = 0;
    int r = width_ - 1;
    int b = 0;
    int t = height_ - 1;
    int rc = eCenter;
    if (x < l) rc |= ReginCode::eLeft;
    if (x > r) rc |= ReginCode::eRight;
    if (y < b) rc |= ReginCode::eButtom;
    if (y > t) rc |= ReginCode::eTop;
    return rc;
}

bool Renderer::lineClip(float ix1, float iy1, float ix2, float iy2, float &ox1, float &oy1, float &ox2, float &oy2) {
    // Cohen-Sutherland algorithm to clip the line
    int l = 0;
    int r = width_ - 1;
    int b = 0;
    int t = height_ - 1;

    int rc1 = encode(ix1, iy1);
    int rc2 = encode(ix2, iy2);
    bool accept = false;

    while (1) {
        if (rc1 == 0 && rc2 == 0) {
            accept = true;
            break;
        } else if (rc1 & rc2) {
            accept = false;
            break;
        } else {
            int outcode = (rc1 == 0) ? rc2 : rc1;
            float x, y;

            if (outcode & ReginCode::eTop) {
                x = ix1 + (ix2 - ix1) * (t - iy1) / float(iy2 - iy1);
                y = t;
            } else if (outcode & ReginCode::eButtom) {
                x = ix1 + (ix2 - ix1) * (b - iy1) / float(iy2 - iy1);
                y = b;
            } else if (outcode & ReginCode::eRight) {
                y = iy1 + (iy2 - iy1) * (r - ix1) / float(ix2 - ix1);
                x = r;
            } else if (outcode & ReginCode::eLeft) {
                y = iy1 + (iy2 - iy1) * (l - ix1) / float(ix2 - ix1);
                x = l;
            }

            if (outcode == rc1) {
                ix1 = x;
                iy1 = y;
                rc1 = encode(ix1, iy1);
            } else {
                ix2 = x;
                iy2 = y;
                rc2 = encode(ix2, iy2);
            }
        }
    }

    ox1 = ix1;
    ox2 = ix2;
    oy1 = iy1;
    oy2 = iy2;

    return accept;
}

void Renderer::drawMesh(const std::array<size_t,3> &indices, std::vector<vec4> &image, const vec3 &color) {
    vec4 v0 = std::get<vec4>(shaderContextBuffer_[indices[0]][0]);
    vec4 v1 = std::get<vec4>(shaderContextBuffer_[indices[1]][0]);
    vec4 v2 = std::get<vec4>(shaderContextBuffer_[indices[2]][0]);

    // Reciprocal of the Homogeneous W
    float rhws[3] = {1.f / v0.w, 1.f / v1.w, 1.f / v2.w};

    // Perspective division
    v0 = v0 * rhws[0];
    v1 = v1 * rhws[1];
    v2 = v2 * rhws[2];

    // Face culling
    vec3 faceFront = cross((v1 - v0).xyz(), (v2 - v0).xyz());
    if (faceFront.z < 0) return;

    // Viewport transform
    v0 = viewport_ * v0;
    v1 = viewport_ * v1;
    v2 = viewport_ * v2;

    // Rasterization
    float x1, x2, y1, y2;
    if (lineClip(v0.x, v0.y, v1.x, v1.y, x1, y1, x2, y2)) {
        drawline(x1, y1, x2, y2, image, color);
    }
    if (lineClip(v1.x, v1.y, v2.x, v2.y, x1, y1, x2, y2)) {
        drawline(x1, y1, x2, y2, image, color);
    }
    if (lineClip(v2.x, v2.y, v0.x, v0.y, x1, y1, x2, y2)) {
        drawline(x1, y1, x2, y2, image, color);
    }
}

void Renderer::drawPrimitive(
        const std::array<size_t,3> &indices, IShader &shader, 
        std::vector<vec4> &image
    ) {
    vec4 v0 = std::get<vec4>(shaderContextBuffer_[indices[0]][0]);
    vec4 v1 = std::get<vec4>(shaderContextBuffer_[indices[1]][0]);
    vec4 v2 = std::get<vec4>(shaderContextBuffer_[indices[2]][0]);

    // Reciprocal of the Homogeneous W
    double rhws[3] = {1 / v0.w, 1 / v1.w, 1 / v2.w};

    // Perspective division
    v0 = v0 * rhws[0];
    v1 = v1 * rhws[1];
    v2 = v2 * rhws[2];

    // Face culling
    vec3 faceFront = cross((v1 - v0).xyz(), (v2 - v0).xyz());
    if (faceFront.z < 0) return;

    // Viewport transform
    v0 = viewport_ * v0;
    v1 = viewport_ * v1;
    v2 = viewport_ * v2;

    // Rasterization
    std::array<vec4,3> fragCoords = {v0, v1, v2};
    auto bbox = bounding_box(fragCoords);
    bbox[0].x = clamp(int(bbox[0].x + 0.5f), 0, width_-1);
    bbox[0].y = clamp(int(bbox[0].y + 0.5f), 0, height_-1);
    bbox[1].x = clamp(int(bbox[1].x + 0.5f), 0, width_-1);
    bbox[1].y = clamp(int(bbox[1].y + 0.5f), 0, height_-1);

    for (int x = bbox[0].x; x <= bbox[1].x; ++ x) {
        for (int y = bbox[0].y; y <= bbox[1].y; ++ y) {
            if (x >= width_ || y >= height_ || x < 0 || y < 0) continue;
            int idx = y * width_ + x;
            int px = x + 0.5;
            int py = y + 0.5;

            auto [a, b, c] = computeBarycentric(px, py, fragCoords);
            if (a < 0 || b < 0 || c < 0) continue;   // Determine if the pixel is inside triangle

            float rhw = (rhws[0] * a + rhws[1] * b + rhws[2] * c);
            float w = 1.f / rhw;
            float c0 = a * rhws[0] * w;
            float c1 = b * rhws[1] * w;
            float c2 = c * rhws[2] * w;
            float z_interp = rhw; // depth value

            if (zbuffer_.at(idx) > z_interp) continue; // Depth test
            zbuffer_.at(idx) = z_interp;  // Updata z-buffer

            interpolation(indices, {c0, c1, c2}, shader.context);

            // Fragment shader
            vec3 frag_color;
            if (shader.fragment(frag_color)) continue; // Fragment shader process
            image[idx] = frag_color;
        }
    }
}

void Renderer::interpolation(const std::array<size_t,3> &idx, const vec3 &bar, ShaderContext &interpolated) {
    auto& ctx0 = shaderContextBuffer_[idx[0]];
    auto& ctx1 = shaderContextBuffer_[idx[1]];
    auto& ctx2 = shaderContextBuffer_[idx[2]];

    // Iterate through all keys in the first context
    for (const auto& [key, value] : ctx0) {
        const auto &k = key;
        std::visit([&](const auto& val0) {
            using T = std::decay_t<decltype(val0)>;
            const auto& val1 = std::get<T>(ctx1.at(k));
            const auto& val2 = std::get<T>(ctx2.at(k));
            interpolated[k] = interp(val0, val1, val2, bar.x, bar.y, bar.z);
        }, value);
    }
}