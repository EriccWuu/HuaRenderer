#include <iostream>
#include "CoreLayer/MathLib.h"
#include "CoreLayer/Model.h"
#include "RenderingLayer/SoftwareRasterizer/Camera.h"
#include "RenderingLayer/SoftwareRasterizer/Renderer.h"
#include <chrono>
#include "SDL.h"

// This porgram uses right-hand coordinate

const std::string ROOT_PATH = "E:/Doc/ProjectFiles/CPP/HuaRenderer/renderer/";

const int width = 1200, height = 900;
// const int width = 900, height = 900;
const double fov = 60, aspectRatio = float(width) / height;
const double near = 0.5, far = 100;
// vec3 camPos = 1.2*vec3(1, 1, 3);
// vec3 up(0, 1, 0);
vec3 lightPos = 3*vec3(1, 2, 1);
vec3 center(0, 0, 0);
vec3 lightdir = (center - lightPos).normalize();

vec3 camPos = 2*vec3(0, 1, 2);
vec3 camDir = vec3{0, 1, 0} - camPos;
vec3 up(0, 1, 0);

// vec3 camPos = 1.5*vec3(0, 1, 0);
// vec3 up(0, 0, -1);

double degSpeed = 60;

std::vector<vec4> img_shadowmap(width * height, ZERO_VEC3);
std::vector<vec4> img_output(width * height, ZERO_VEC3);
Texture shadowmap(width, height);

Model African_head(ROOT_PATH + "assets/models/african_head/african_head.obj");
Model African_head_eye_inner(ROOT_PATH + "assets/models/african_head/african_head_eye_inner.obj");

Model Diablo3_pose(ROOT_PATH + "assets/models/diablo3_pose/diablo3_pose.obj");

// Model boggie_body(ROOT_PATH + "assets/models/boggie/body.obj");
// Model boggie_eyes(ROOT_PATH + "assets/models/boggie/eyes.obj");
// Model boggie_head(ROOT_PATH + "assets/models/boggie/head.obj");

Model ground(ROOT_PATH + "assets/models/floor/floor.obj");

Model red_bird(ROOT_PATH + "assets/models/Red/Red.obj");
// Texture red_bird_texture(ROOT_PATH + "assets/models/Red/Red.tga");

std::vector<Model> scene;

Camera camera;
Renderer renderer(width, height);
DepthShader shadowShader;
Shader modelShader;

mat4 MODEL;
mat4 VIEW;
mat4 MV;
mat4 NORMAL;
mat4 PROJECTION;
mat4 VIEWPORT;

auto startTime = std::chrono::high_resolution_clock::now();

bool meshMode = false;

void init() {
    camera = Camera(camPos, camDir, up, width, height, fov, aspectRatio, near, far);

    MODEL = E44;
    VIEW = lookAt(lightPos, lightdir, up);
    MV = VIEW * MODEL;
    NORMAL = MV.invertTranspose();
    PROJECTION = projection(camera.fov, camera.aspectRatio, camera.near, camera.far);
    VIEWPORT = viewport(width, height);

    renderer.setViewport(VIEWPORT);

    shadowShader.uniform_MODEL = MODEL;
    shadowShader.uniform_VIEW = lookAt(lightPos, lightdir, up);
    shadowShader.uniform_NORMAL = shadowShader.uniform_MODEL.invertTranspose();
    shadowShader.uniform_PROJECTION = PROJECTION;
    shadowShader.uniform_VIEWPORT = VIEWPORT;

    modelShader.uniform_MODEL = MODEL;
    modelShader.uniform_VIEW = lookAt(camera.position.xyz(), camera.direction.xyz(), camera.up.xyz());
    modelShader.uniform_NORMAL = modelShader.uniform_MODEL.invertTranspose();
    modelShader.uniform_PROJECTION = projection(camera.fov, camera.aspectRatio, camera.near, camera.far);
    modelShader.uniform_VIEWPORT = VIEWPORT;
    modelShader.uniform_Mshadow = shadowShader.uniform_PROJECTION * shadowShader.uniform_VIEW;
    modelShader.uniform_lightDir = (vec4(lightdir)).xyz().normalize();
    modelShader.uniform_lightPos = lightPos;

}

void loadScene() {
    // scene.push_back(African_head);
    // scene.push_back(African_head_eye_inner);
    // scene.push_back(Diablo3_pose);
    // scene.push_back(boggie_head);
    // scene.push_back(boggie_eyes);
    // scene.push_back(boggie_body);
    scene.push_back(red_bird);
}

void testShadowMap() {
    auto start_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(start_time - startTime);
    float rotateAngle = degSpeed * duration.count() / 1000;

    img_shadowmap = std::vector<vec4>(width * height, ONE_VEC3);
    std::fill(img_output.begin(), img_output.end(), ZERO_VEC3);

    auto Mold = translate({0, 0, -2})* rotateY(0) * scale({4, 1, 4});
    auto Mnew = translate({0, 1, 0}) * rotateY(rotateAngle) * scale({1, 1, 1});

    /***************************************************************
    *            First pass: rendering the shadow map              *
    ***************************************************************/
    renderer.initZbuffer();

    shadowShader.uniform_MODEL = Mold;
    shadowShader.uniform_NORMAL = shadowShader.uniform_MODEL.invertTranspose();

    renderer.bindVertexBuffer(ground.uniqueVertices());
    renderer.bindIndexBuffer(ground.uniqueIndices());
    renderer.draw(shadowShader, img_shadowmap);

    shadowShader.uniform_MODEL = Mnew;
    shadowShader.uniform_NORMAL = shadowShader.uniform_MODEL.invertTranspose();
    for (auto &model : scene) {
        renderer.bindVertexBuffer(model.uniqueVertices());
        renderer.bindIndexBuffer(model.uniqueIndices());
        renderer.draw(shadowShader, img_shadowmap);
    }

    shadowmap.swap(img_shadowmap, width, height);
    modelShader.uniform_shadowmapPtr = std::make_shared<Texture>(shadowmap);

    /***************************************************************
    *       Second pass: rendering objects with shadow map         *
    ***************************************************************/
    renderer.initZbuffer();

    renderer.meshMode(meshMode);

    modelShader.uniform_MODEL = Mold;
    modelShader.uniform_NORMAL = modelShader.uniform_MODEL.invertTranspose();

    modelShader.uniform_diffuseTexPtr = ground.diffuseTexture();
    modelShader.uniform_specularTexPtr = ground.specularTexture();
    modelShader.uniform_normalTexPtr = ground.normalTexture();
    renderer.bindVertexBuffer(ground.uniqueVertices());
    renderer.bindIndexBuffer(ground.uniqueIndices());
    renderer.draw(modelShader, img_output);

    modelShader.uniform_MODEL = Mnew;
    modelShader.uniform_NORMAL = modelShader.uniform_MODEL.invertTranspose();
    for (auto &model :  scene) {
        modelShader.uniform_diffuseTexPtr = model.diffuseTexture();
        modelShader.uniform_specularTexPtr = model.specularTexture();
        modelShader.uniform_normalTexPtr = model.normalTexture();
        renderer.bindVertexBuffer(model.uniqueVertices());
        renderer.bindIndexBuffer(model.uniqueIndices());
        renderer.draw(modelShader, img_output);
    }

}

void present() {
    SDL_Init(SDL_INIT_EVERYTHING);
    SDL_Window* window = SDL_CreateWindow("HuaRenderer",
        SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, 
        width, height, 
        SDL_WINDOW_SHOWN | SDL_WINDOW_RESIZABLE);
    SDL_Renderer* renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED);
    SDL_Texture* texture = SDL_CreateTexture(renderer, SDL_PIXELFORMAT_RGB888, SDL_TEXTUREACCESS_STREAMING, width, height);

    bool isQuit = false;
    SDL_Event event;

    while (!isQuit) {
        // Dealing with events
        while (SDL_PollEvent(&event)) {
            if (event.type == SDL_QUIT) {
                isQuit = true;
            }
            if (event.type == SDL_KEYDOWN) {
                if (event.key.keysym.sym == SDLK_SPACE) {
                    meshMode = !meshMode;
                }
            }
        }

        auto start_time = std::chrono::high_resolution_clock::now();

        testShadowMap();

        void* pixels;
        int pitch;
        SDL_LockTexture(texture, NULL, &pixels, &pitch);

        for (int y = 0; y < height; ++ y) {
            for (int x = 0; x < width; ++ x) {
                int xx = x, yy = height - 1 - y;
                auto color = img_output[xx + width * yy];
                uint32_t* pixelPtr = (uint32_t*)((uint8_t*)pixels + y * pitch) + x;
                *pixelPtr = (uint8_t(color.x) << 16) | (uint8_t(color.y) << 8) | uint8_t(color.z);
            }
        }

        SDL_UnlockTexture(texture);

        // Clear renderer
        SDL_RenderClear(renderer);

        // 
        SDL_RenderCopy(renderer, texture, NULL, NULL);

        // Render the image
        SDL_RenderPresent(renderer);

        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
        std::cout << "Time cost: " << duration.count() << " ms" << std::endl;

    }

    SDL_DestroyTexture(texture);
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();
}

int main(int argc, char** argv) {

    loadScene();

    init();

    present();

    TGAImage output_image(img_output, width, height);
    output_image.write_tga_file("output.tga");

    return 0;
}