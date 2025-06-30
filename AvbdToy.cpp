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
#include "World.h"
#include "physics/BasicObjects.h"
#include "physics/World.h"

// ---------------------------------------------------------------
// Tunables
// ---------------------------------------------------------------
static constexpr int   WINDOW_W = 960;
static constexpr int   WINDOW_H = 720;
static constexpr float SCALE     = 60.0f;   // world‑unit → pixel
static constexpr float DT        = 1.0f / 60.0f;

// Convert world coordinates (meters) to pixel space
static inline void WorldToScreen(const Vector2 &p, int &sx, int &sy)
{
    sx = static_cast<int>(p.X * SCALE + WINDOW_W * 0.5f);
    sy = static_cast<int>(WINDOW_H * 0.75f - p.Y * SCALE); // origin at ground
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
    // Construct a simple demo world
    // -----------------------------------------------------------
    World world;

    // Two falling discs stacked — tweak as needed
    const float radius = 0.25f;
    world.AddBody({ 0.0f, 5.0f }, radius, 1.0f);
    world.AddBody({ 0.0f, 6.0f }, radius, 1.0f);

    // Ground contact represented by a constraint per body
    const Vector2 groundNormal(0.0f, 1.0f);
    for (BodyIndex bi = 1; bi <= 2; ++bi)
        world.AddConstraint(bi, INVALID, ConstraintType::Normal,
                            groundNormal,   /*rest=*/0.0f,
                            world.Settings.KStart);

    // -----------------------------------------------------------
    // Main loop
    // -----------------------------------------------------------
    bool quit = false;
    auto last = std::chrono::high_resolution_clock::now();

    while (!quit)
    {
        // --- handle events
        SDL_Event e;
        while (SDL_PollEvent(&e)) {
            if (e.type == SDL_QUIT) quit = true;
        }

        // --- fixed‑timestep update
        auto now   = std::chrono::high_resolution_clock::now();
        float acc  = std::chrono::duration<float>(now - last).count();
        while (acc >= DT) {
            world.Step(DT);
            acc -= DT;
        }
        last = now;

        // --- render
        SDL_SetRenderDrawColor(ren, 30, 30, 30, 255);
        SDL_RenderClear(ren);

        SDL_SetRenderDrawColor(ren, 255, 160, 40, 255);
        for (BodyIndex bi = 0; bi <= world.GetNumBodies(); ++bi)
        {
            const Body &b = *world.GetBody(bi);
            int sx, sy;
            WorldToScreen(b.Pos, sx, sy);
            FillCircle(ren, sx, sy, static_cast<int>(b.Radius * SCALE));
        }

        // ground line for reference
        SDL_SetRenderDrawColor(ren, 80, 80, 80, 255);
        SDL_RenderDrawLine(ren, 0, WINDOW_H * 0.75f, WINDOW_W, WINDOW_H * 0.75f);

        SDL_RenderPresent(ren);

        SDL_Delay(1); // yield to OS
    }

    SDL_DestroyRenderer(ren);
    SDL_DestroyWindow(win);
    SDL_Quit();
    return 0;
}
