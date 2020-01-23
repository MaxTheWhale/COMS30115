#include "SDLauxiliary.h"
#include "TestModel.h"
#include <SDL2/SDL.h>
#include <glm/glm.hpp>
#include <iostream>
#include <stdint.h>

using namespace std;
using glm::mat3;
using glm::vec3;

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define SCREEN_WIDTH 320
#define SCREEN_HEIGHT 256
#define CAM_VELOCITY 0.4f
#define STAR_VELOCITY 0.5f
#define MOUSE_SENSITIVITY 0.0015f;
#define FULLSCREEN_MODE false

/* ----------------------------------------------------------------------------*/
/* GLOBAL VARIABLES                                                            */
int t;
vec3 camera_pos;
vec3 camera_rot;
vec3 cam_forward;
vec3 cam_right;
mat3 R;
vector<vec3> stars(100000);

/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */
void Update();
void Draw(screen *screen);
void Interpolate(float a, float b, vector<float> &result);
void Interpolate(vec3 a, vec3 b, vector<vec3> &result);
void HandleKeyboard(const Uint8 *state, float dt);
void HandleMouse();
void UpdateStars(float dt);
void ClampCamera(vec3 &rot);
void CalcRMatrix(vec3 rot, mat3 &R);

int main(int argc, char *argv[]) {

  screen *screen = InitializeSDL(SCREEN_WIDTH, SCREEN_HEIGHT, FULLSCREEN_MODE);
  SDL_SetRelativeMouseMode(SDL_TRUE);
  t = SDL_GetTicks(); /*Set start value for timer.*/
  for (size_t s = 0; s < stars.size(); ++s) {
    stars[s].x = float(rand()) * 2.0f / float(RAND_MAX) - 1.0f;
    stars[s].y = float(rand()) * 2.0f / float(RAND_MAX) - 1.0f;
    stars[s].z = float(rand()) / float(RAND_MAX);
  }
  camera_pos = vec3(0, 0, 0);
  camera_rot = vec3(0, 0, 0);
  CalcRMatrix(camera_rot, R);
  cam_forward = glm::transpose(R) * vec3(0, 0, 1);
  cam_right = glm::transpose(R) * vec3(1, 0, 0);

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
    vec3 pos = R * (stars[s] - camera_pos);
    if (pos.z <= 0) continue;
    int px = f * pos.x / pos.z + SCREEN_WIDTH / 2;
    int py = f * pos.y / pos.z + SCREEN_HEIGHT / 2;
    if (px > 0 && px < SCREEN_WIDTH && py > 0 && py < SCREEN_HEIGHT) {
      vec3 colour = 0.2f * vec3(1, 1, 1) / (pos.z * pos.z);
      PutPixelSDL(screen, px, py, colour);
    }
  }
}

/*Place updates of parameters here*/
void Update() {
  /* Compute frame time */
  int t2 = SDL_GetTicks();
  float dt = float(t2 - t) / 1000.0f;
  t = t2;
  /*Good idea to remove this*/
  std::cout << "Render time: " << dt * 1000 << " ms." << std::endl;
  /* Update variables*/
  HandleKeyboard(SDL_GetKeyboardState(NULL), dt);
  HandleMouse();
  UpdateStars(dt);
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

void HandleKeyboard(const Uint8 *state, float dt) {
  float ds = dt * CAM_VELOCITY;
  if (state[SDL_SCANCODE_LSHIFT])
    ds *= 2;
  if (state[SDL_SCANCODE_A])
    camera_pos -= cam_right * ds;
  if (state[SDL_SCANCODE_D])
    camera_pos += cam_right * ds;
  if (state[SDL_SCANCODE_W])
    camera_pos += cam_forward * ds;
  if (state[SDL_SCANCODE_S])
    camera_pos -= cam_forward * ds;
  if (state[SDL_SCANCODE_LCTRL])
    camera_pos.y += ds;
  if (state[SDL_SCANCODE_SPACE])
    camera_pos.y -= ds;
}

void HandleMouse() {
  int motion_x = 0;
  int motion_y = 0;
  SDL_GetRelativeMouseState(&motion_x, &motion_y);
  if (motion_x + motion_y) {
    camera_rot.y += motion_x * MOUSE_SENSITIVITY;
    camera_rot.x -= motion_y * MOUSE_SENSITIVITY;
    ClampCamera(camera_rot);
    CalcRMatrix(camera_rot, R);
    cam_forward = glm::transpose(R) * vec3(0, 0, 1);
    cam_right = glm::transpose(R) * vec3(1, 0, 0);
  }
}

void UpdateStars(float dt) {
  for (size_t s = 0; s < stars.size(); ++s) {
    stars[s].z += dt * STAR_VELOCITY;
    if (stars[s].z <= 0) {
      stars[s].z += 1;
      stars[s].x = float(rand()) * 2.0f / float(RAND_MAX) - 1.0f;
      stars[s].y = float(rand()) * 2.0f / float(RAND_MAX) - 1.0f;
    }
    if (stars[s].z > 1) {
      stars[s].z -= 1;
      stars[s].x = float(rand()) * 2.0f / float(RAND_MAX) - 1.0f;
      stars[s].y = float(rand()) * 2.0f / float(RAND_MAX) - 1.0f;
    }
  }
}

void ClampCamera(vec3 &rot) {
  if (camera_rot.x > M_PI / 2)
    camera_rot.x = M_PI / 2;
  if (camera_rot.x < -M_PI / 2)
    camera_rot.x = -M_PI / 2;
  if (camera_rot.y > M_PI)
    camera_rot.y -= 2 * M_PI;
  if (camera_rot.y < M_PI)
    camera_rot.y += 2 * M_PI;
}

void CalcRMatrix(vec3 rot, mat3 &R) {
  float sin_x = sinf(rot.x);
  float cos_x = cosf(rot.x);
  float sin_y = sinf(rot.y);
  float cos_y = cosf(rot.y);
  mat3 Rx = mat3(1, 0, 0, 0, cos_x, -sin_x, 0, sin_x, cos_x);
  mat3 Ry = mat3(cos_y, 0, sin_y, 0, 1, 0, -sin_y, 0, cos_y);
  R = Rx * Ry;
};