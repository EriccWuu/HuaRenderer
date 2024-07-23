#include <iostream>
#include "SDL.h"

int main(int argc, char* argv[]) {

    std::cout << "Hello world!" << std::endl;
    SDL_Init(SDL_INIT_EVERYTHING);
    SDL_Window* window = SDL_CreateWindow("HuaRenderer", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, 800, 800, SDL_WINDOW_SHOWN);
    SDL_Renderer* renderer = SDL_CreateRenderer(window, -1, 0);

    bool isQuit = false;
    SDL_Event event;

    SDL_SetRenderDrawColor(renderer, 0, 255, 0, 255);
    while (!isQuit) {
        while (SDL_PollEvent(&event)) {
            if (event.type == SDL_QUIT) {
                isQuit = true;
            }
        }

        SDL_Rect rect;
        rect.x = 400;
        rect.y = 400;
        rect.w = 200;
        rect.h = 100;
        SDL_RenderDrawRect(renderer, &rect);
        SDL_RenderPresent(renderer);
    }

    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();
    return 0;
}