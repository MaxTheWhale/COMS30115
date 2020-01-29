#include <ModelTriangle.h>
#include <CanvasTriangle.h>
#include <DrawingWindow.h>
#include <Utils.h>
#include <fstream>
#include <glm/glm.hpp>
#include <vector>

using namespace std;
using namespace glm;

#define WIDTH 320
#define HEIGHT 240

void draw();
void update();
void handleEvent(SDL_Event event);
vector<float> Interpolate(float a, float b, int n);
vector<vec3> Interpolate(vec3 a, vec3 b, int n);

DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);

int main(int argc, char *argv[]) {
  SDL_Event event;
  while (true) {
    // We MUST poll for events - otherwise the window will freeze !
    if (window.pollForInputEvents(&event))
      handleEvent(event);
    update();
    draw();
    // Need to render the frame at the end, or nothing actually gets shown on
    // the screen !
    window.renderFrame();
  }
}

void draw() {
  window.clearPixels();
  vec3 red = vec3(255, 0, 0);
  vec3 green = vec3(0, 255, 0);
  vec3 blue = vec3(0, 0, 255);
  vec3 yellow = vec3(255, 255, 0);
  vector<vec3> leftSide = Interpolate(red, yellow, HEIGHT);
  vector<vec3> rightSide = Interpolate(blue, green, HEIGHT);
  for (int y = 0; y < window.height; y++) {
    vector<vec3> line = Interpolate(leftSide[y], rightSide[y], WIDTH);
    for (int x = 0; x < window.width; x++) {
      uint32_t colour = (255 << 24) + (int(line[x].r) << 16) +
                        (int(line[x].g) << 8) + int(line[x].b);
      window.setPixelColour(x, y, colour);
    }
  }
}

void update() {
  // Function for performing animation (shifting artifacts or moving the camera)
}

void handleEvent(SDL_Event event) {
  if (event.type == SDL_KEYDOWN) {
    if (event.key.keysym.sym == SDLK_LEFT)
      cout << "LEFT" << endl;
    else if (event.key.keysym.sym == SDLK_RIGHT)
      cout << "RIGHT" << endl;
    else if (event.key.keysym.sym == SDLK_UP)
      cout << "UP" << endl;
    else if (event.key.keysym.sym == SDLK_DOWN)
      cout << "DOWN" << endl;
  } else if (event.type == SDL_MOUSEBUTTONDOWN)
    cout << "MOUSE CLICKED" << endl;
}

vector<float> Interpolate(float a, float b, int n) {
  vector<float> result;
  if (n == 1) {
    result.push_back((a + b) / 2);
    return result;
  } else if (n > 1) {
    float step = (b - a) / (n - 1);
    for (int i = 0; i < n; ++i) {
      result.push_back(a);
      a += step;
    }
  }
  return result;
}

vector<vec3> Interpolate(vec3 a, vec3 b, int n) {
  vector<vec3> result;
  if (n == 1) {
    result.push_back((a + b) / 2.0f);
    return result;
  } else if (n > 1) {
    vec3 step = (b - a) / (n - 1.0f);
    for (int i = 0; i < n; ++i) {
      result.push_back(a);
      a += step;
    }
  }
  return result;
}
