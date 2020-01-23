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
#define CAM_VELOCITY 0.4f
#define STAR_VELOCITY 0.5f
#define FULLSCREEN_MODE false

/* ----------------------------------------------------------------------------*/
/* GLOBAL VARIABLES                                                            */
int t;
vec3 camera_pos;
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
  camera_pos = vec3(0,0,0);

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
    float z = stars[s].z - camera_pos.z;
    if (z <= 0) continue;
    float x = stars[s].x - camera_pos.x;
    float y = stars[s].y - camera_pos.y;
    int px = f * x / z + SCREEN_WIDTH / 2;
    int py = f * y / z + SCREEN_HEIGHT / 2;
    if (px > 0 && px < SCREEN_WIDTH && py > 0 && py < SCREEN_HEIGHT) {
      vec3 colour = 0.2f * vec3(1, 1, 1) / (z * z);
      PutPixelSDL(screen, px, py, colour);
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
  const Uint8 *state = SDL_GetKeyboardState(NULL);
  float ds = dt / 1000.0f * CAM_VELOCITY;
  if (state[SDL_SCANCODE_A]) {
    camera_pos.x -= ds;
  }
  if (state[SDL_SCANCODE_D]) {
    camera_pos.x += ds;
  }
  if (state[SDL_SCANCODE_W]) {
    camera_pos.z += ds;
  }
  if (state[SDL_SCANCODE_S]) {
    camera_pos.z -= ds;
  }
  if (state[SDL_SCANCODE_LSHIFT]) {
    camera_pos.y += ds;
  }
  if (state[SDL_SCANCODE_SPACE]) {
    camera_pos.y -= ds;
  }
  for (size_t s = 0; s < stars.size(); ++s) {
    stars[s].z += dt / 1000.0f * STAR_VELOCITY;
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