// avbd_viewer.cpp — minimal SDL2 frontend to visualise the World
// ---------------------------------------------------------------
// Build (Linux/Mac):
//   g++ avbd_viewer.cpp -o viewer -lSDL2 -std=c++17
// Build (MSVC / Windows):
//   cl /std:c++17 /EHsc avbd_viewer.cpp /Ipath\to\SDL2\include path\to\SDL2.lib
// ---------------------------------------------------------------
// The file expects your existing physics headers in the include path:
//   * World.h      (contains World, WorldSettings, etc.)
//   * BasicObjects.h and friends
// ---------------------------------------------------------------
// The program creates a window, steps the simulation (fixed‑timestep),
// and draws each Body as a filled circle.  Adjust SCALE, WINDOW_* or
// add extra rendering (constraints, contacts) as you wish.

#include <SDL2/SDL.h>
#include <chrono>
#include <cmath>
#include "physics/World.h"
#include "physics/BasicObjects.h"

// ---------------------------------------------------------------
// Tunables
// ---------------------------------------------------------------
static constexpr int   WINDOW_W = 960;
static constexpr int   WINDOW_H = 720;
static constexpr float SCALE     = 20.0f;   // world‑unit → pixel
static constexpr float DT        = 1.0f / 60.0f;

// Convert world coordinates (meters) to pixel space
static inline void WorldToScreen(const Vector2 &p, int &sx, int &sy)
{
    sx = static_cast<int>(p.X * SCALE);
    sy = static_cast<int>(WINDOW_H - p.Y * SCALE); // origin at ground
}

// Naïve filled‑circle helper (Bresenham variant)
static void FillCircle(SDL_Renderer *r, int cx, int cy, int radius)
{
    for (int dy = -radius; dy <= radius; ++dy)
    {
        int dx = static_cast<int>(std::sqrt(radius*radius - dy*dy));
        SDL_RenderDrawLine(r, cx - dx, cy + dy, cx + dx, cy + dy);
    }
}

static World *Initialize()
{
    World *world = new World();

    // Two falling discs stacked — tweak as needed
    for (int i = 0; i < 200; ++i)
    {
        static std::mt19937 rng{ std::random_device{}() };
        std::uniform_real_distribution<float> unit(0.0f, 1.0f);

        float radius = unit(rng) * 2.f + 0.2f;
        world->AddBody({ unit(rng) * 50, unit(rng) * 50 }, {unit(rng) * 10 - 5, unit(rng) * 10 - 5}, radius, 3.14f * radius * radius);
    }

    world->AddPlane(Vector2(0, 1.0), 2);
    world->AddPlane(Vector2(1.0, 0.0), 1);
    world->AddPlane(Vector2(-1.0, 0.0), -45);

    return world;
}

int main(int argc, char **argv)
{
    if (SDL_Init(SDL_INIT_VIDEO) != 0) {
        SDL_Log("SDL init failed: %s", SDL_GetError());
        return 1;
    }

    SDL_Window   *win = SDL_CreateWindow("AVBD Viewer",
                                         SDL_WINDOWPOS_CENTERED,
                                         SDL_WINDOWPOS_CENTERED,
                                         WINDOW_W, WINDOW_H,
                                         SDL_WINDOW_SHOWN);
    SDL_Renderer *ren = SDL_CreateRenderer(win, -1, SDL_RENDERER_ACCELERATED);

    // -----------------------------------------------------------
    // Main loop
    // -----------------------------------------------------------
    bool quit = false;
    auto last = std::chrono::high_resolution_clock::now();

    long frame = 0;

    World *world = Initialize();

    float acc = 0;
    while (!quit)
    {
        // --- handle events
        SDL_Event e;
        while (SDL_PollEvent(&e)) {
            if (e.type == SDL_QUIT)
            {
                quit = true;
            }
            else if (e.type == SDL_KEYDOWN)
            {
                if (e.key.keysym.sym == SDLK_SPACE)
                {
                    delete world;
                    world = Initialize();
                }
            }
        }

        // --- fixed‑timestep update
        auto now   = std::chrono::high_resolution_clock::now();
        acc += std::chrono::duration<float>(now - last).count();
        while (acc >= DT) {
            world->Step(DT);
            acc -= DT;
        }
        last = now;

        // --- render
        SDL_SetRenderDrawColor(ren, 30, 30, 30, 255);
        SDL_RenderClear(ren);

        for (BodyIndex bi = 0; bi < world->GetNumBodies(); ++bi)
        {
            const Body &b = *world->GetBody(bi);
            SDL_SetRenderDrawColor(ren, b.Color[0],b.Color[1],b.Color[2], 255);

            int sx, sy;
            WorldToScreen(b.Pos, sx, sy);
            FillCircle(ren, sx, sy, static_cast<int>(b.Radius * SCALE));
        }

        SDL_SetRenderDrawColor(ren, 80, 80, 80, 255);
        for (BodyIndex pi = 0; pi < world->GetNumPlanes(); ++pi)
        {
            const Plane &p = world->GetPlane(pi);

            Vector2 tangent(-p.Normal.Y, p.Normal.X);
            constexpr float len = 100.f;
            Vector2 p0 = tangent * len + p.Normal * p.Distance;
            Vector2 p1 = tangent * -1 * len + p.Normal * p.Distance;

            int startX, startY;
            int endX, endY;
            WorldToScreen(p0, startX, startY);
            WorldToScreen(p1, endX, endY);
            SDL_RenderDrawLine(ren, startX, startY, endX, endY);
        }

        SDL_RenderPresent(ren);
        SDL_Delay(1); // yield to OS
    }

    SDL_DestroyRenderer(ren);
    SDL_DestroyWindow(win);
    SDL_Quit();
    return 0;
}
