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

/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */

void Update();
void Draw(screen *screen);
void Interpolate(float a, float b, vector<float> &result);
void Interpolate(vec3 a, vec3 b, vector<vec3> &result);

int main(int argc, char *argv[]) {

  screen *screen = InitializeSDL(SCREEN_WIDTH, SCREEN_HEIGHT, FULLSCREEN_MODE);
  t = SDL_GetTicks(); /*Set start value for timer.*/

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

  vec3 topLeft(1, 0, 0);     // red
  vec3 topRight(0, 0, 1);    // blue
  vec3 bottomRight(0, 1, 0); // green
  vec3 bottomLeft(1, 1, 0);  // yellow
  vector<vec3> leftSide(SCREEN_HEIGHT);
  vector<vec3> rightSide(SCREEN_HEIGHT);
  vector<vec3> line(SCREEN_WIDTH);
  Interpolate(topLeft, bottomLeft, leftSide);
  Interpolate(topRight, bottomRight, rightSide);
  for (int y = 0; y < SCREEN_HEIGHT; ++y) {
    Interpolate(leftSide[y], rightSide[y], line);
    for (int x = 0; x < SCREEN_WIDTH; ++x) {
      PutPixelSDL(screen, x, y, line[x]);
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
  float current = a;
  for (int i = 0; i < size; ++i) {
    result[i] = current;
    current += step;
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
  vec3 current = a;
  for (int i = 0; i < size; ++i) {
    result[i] = current;
    current += step;
  }
}