#include "SDLauxiliary.h"
#include "TestModel.h"
#include <SDL2/SDL.h>
#include <glm/glm.hpp>
#include <iostream>
#include <stdint.h>

using namespace std;
using glm::mat3;
using glm::vec3;

#define SCREEN_WIDTH 320
#define SCREEN_HEIGHT 256
#define FULLSCREEN_MODE false

/* ----------------------------------------------------------------------------*/
/* GLOBAL VARIABLES                                                            */
int t;
vector<vec3> stars(1000);

/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */
void Update();
void Draw(screen *screen);
void Interpolate(float a, float b, vector<float> &result);
void Interpolate(vec3 a, vec3 b, vector<vec3> &result);

int main(int argc, char *argv[]) {

  screen *screen = InitializeSDL(SCREEN_WIDTH, SCREEN_HEIGHT, FULLSCREEN_MODE);
  t = SDL_GetTicks(); /*Set start value for timer.*/
  for (size_t s = 0; s < stars.size(); ++s) {
    stars[s].x = float(rand()) * 2.0f / float(RAND_MAX) - 1.0f;
    stars[s].y = float(rand()) * 2.0f / float(RAND_MAX) - 1.0f;
    stars[s].z = float(rand()) / float(RAND_MAX);
  }

  while (NoQuitMessageSDL()) {
    Draw(screen);
    Update();
    SDL_Renderframe(screen);
  }

  SDL_SaveImage(screen, "screenshot.bmp");

  KillSDL(screen);
  return 0;
}

/*Place your drawing here*/
void Draw(screen *screen) {
  /* Clear buffer */
  memset(screen->buffer, 0, screen->height * screen->width * sizeof(uint32_t));
  for (size_t s = 0; s < stars.size(); ++s) {
    float f = SCREEN_HEIGHT / 2.0f;
    int x = f * stars[s].x / stars[s].z + SCREEN_WIDTH / 2;
    int y = f * stars[s].y / stars[s].z + SCREEN_HEIGHT / 2;
    if (x > 0 && x < SCREEN_WIDTH && y > 0 && y < SCREEN_HEIGHT) {
      vec3 colour = 0.2f * vec3(1, 1, 1) / (stars[s].z * stars[s].z);
      PutPixelSDL(screen, x, y, colour);
    }
  }
}

/*Place updates of parameters here*/
void Update() {
  /* Compute frame time */
  int t2 = SDL_GetTicks();
  float dt = float(t2 - t);
  t = t2;
  /*Good idea to remove this*/
  std::cout << "Render time: " << dt << " ms." << std::endl;
  /* Update variables*/
  float velocity = 0.5f;
  for (size_t s = 0; s < stars.size(); ++s) {
    stars[s].z += dt / 1000.0f * velocity;
    if (stars[s].z <= 0)
      stars[s].z += 1;
    if (stars[s].z > 1)
      stars[s].z -= 1;
  }
}

void Interpolate(float a, float b, vector<float> &result) {
  int size = result.size();
  if (size == 0)
    return;
  if (size == 1) {
    result[0] = (a + b) / 2;
    return;
  }
  float step = (b - a) / (size - 1);
  for (int i = 0; i < size; ++i) {
    result[i] = a;
    a += step;
  }
}

void Interpolate(vec3 a, vec3 b, vector<vec3> &result) {
  int size = result.size();
  if (size == 0)
    return;
  if (size == 1) {
    result[0] = (a + b) / 2.0f;
    return;
  }
  vec3 step = (b - a) / (float)(size - 1);
  for (int i = 0; i < size; ++i) {
    result[i] = a;
    a += step;
  }
}