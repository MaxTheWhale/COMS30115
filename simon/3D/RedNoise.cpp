#include <ModelTriangle.h>
#include <CanvasTriangle.h>
#include <DrawingWindow.h>
#include <Utils.h>
#include <fstream>
#include <glm/glm.hpp>
#include <unordered_map>
#include <vector>
#include <sys/time.h>
#include "Camera.hpp"
#include "Model.hpp"
#include "Times.hpp"
#include "VectorOutput.hpp"
using namespace std;
using namespace glm;

#define WIDTH 320
#define HEIGHT 240

void draw();
void line(CanvasPoint p, CanvasPoint q, int colour);
void triangle(CanvasTriangle t, int colour, bool filled = false);
int *loadPPM(string fileName, int &width, int &height);
void savePPM(string fileName, DrawingWindow *window);
void skipHashWS(ifstream &f);
void update(Camera &cam, vector<Model> models);
void handleEvent(SDL_Event event, Camera &cam);
float depthBuffer[WIDTH * HEIGHT];
vector<float> Interpolate(float a, float b, int n);
vector<vec3> Interpolate(vec3 a, vec3 b, int n);

void drawTriangles(Model &model, Camera &cam)
{
  mat4 MVP = cam.projection * cam.worldToCamera() * model.transform;
  for (auto tri : model.tris)
  {
    vec4 camToV1 = MVP * tri.vertices[0];
    camToV1 /= camToV1.w;
    vec4 camToV2 = MVP * tri.vertices[1];
    camToV2 /= camToV2.w;
    vec4 camToV3 = MVP * tri.vertices[2];
    camToV3 /= camToV3.w;
    CanvasPoint v1 = CanvasPoint(
        (camToV1.x + 1) * 0.5f * WIDTH,
        (1 - (camToV1.y + 1) * 0.5f) * HEIGHT,
        (99.9f / 2) * camToV1.z + (100.1f / 2));
    CanvasPoint v2 = CanvasPoint(
        (camToV2.x + 1) * 0.5f * WIDTH,
        (1 - (camToV2.y + 1) * 0.5f) * HEIGHT,
        (99.9f / 2) * camToV2.z + (100.1f / 2));
    CanvasPoint v3 = CanvasPoint(
        (camToV3.x + 1) * 0.5f * WIDTH,
        (1 - (camToV3.y + 1) * 0.5f) * HEIGHT,
        (99.9f / 2) * camToV3.z + (100.1f / 2));
    if (v1.depth > 0 || v2.depth > 0 || v3.depth > 0)
      triangle(CanvasTriangle(v1, v2, v3), tri.colour.toPackedInt(), true);
  }
}

class Texture
{
public:
  int width, height;
  int *buff;
  Texture(string fileName) { buff = loadPPM(fileName, width, height); }

  void draw(DrawingWindow window)
  {
    for (int y = 0; y < height; y++)
    {
      for (int x = 0; x < width; x++)
      {
        window.setPixelColour(x, y, buff[y * width + x]);
      }
    }
  }
};

void texturedTriangle(CanvasTriangle screenTri, CanvasTriangle texTri,
                      Texture tex);

DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);

//vector<vec3> cameraPositions{ vec3(5.0f, 2.5f, 3.0f), vec3(5.0f, 0.0f, 3.0f), vec3(5.0f, 0.0f, 6.0f) };
vector<mat4> cameraTransforms = vector<mat4>();
int main(int argc, char *argv[])
{
  SDL_Event event;

  Model cornell = Model("cornell-box");

  Camera cam;
  cam.setProjection(90.0f, WIDTH / (float)HEIGHT, 0.1f, 100.0f);
  cam.lookAt(vec3(5.0f, 2.5f, 3.0f), vec3(0.0f, 2.5f, 0.0f));
  cam.moves.push(Movement(cam.transform));
  cam.moves.top().move(vec3(0.0f, -2.0f, 0.0f));
  //cam.moves.top().rotate(vec3(2.0f, 0.0f, 3.0f));

  Times::init();

  while (true)
  {
    //cout << "camera transform = " << cam.transform << endl;
    Times::update();
    //cout << "deltaTime: " << Times::deltaTime() << endl;
    // We MUST poll for events - otherwise the window will freeze !
    if (window.pollForInputEvents(&event))
      handleEvent(event, cam);
    //cout << "deltaTime = " << Times::deltaTime() << endl;
    update(cam, vector<Model>{cornell});
    draw();
    drawTriangles(cornell, cam);
    // Need to render the frame at the end, or nothing actually gets shown on
    // the screen !
    window.renderFrame();
  }
}

void draw()
{
  window.clearPixels();
  for (int i = 0; i < WIDTH * HEIGHT; i++)
  {
    depthBuffer[i] = std::numeric_limits<float>::infinity();
  }
}

int moveStage = 0;

void update(Camera &cam, vector<Model> models)
{
  // Function for performing animation (shifting artifacts or moving the camera)
  cam.update();
  for (unsigned int i = 0; i < models.size(); i++)
  {
    models[i].update();
  }
}

void handleEvent(SDL_Event event, Camera &cam)
{
  if (event.type == SDL_KEYDOWN)
  {
    if (event.key.keysym.sym == SDLK_LEFT)
    {
      cout << "LEFT" << endl;
      cam.move(-0.5f * cam.right);
    }
    else if (event.key.keysym.sym == SDLK_RIGHT)
    {
      cout << "RIGHT" << endl;
      cam.move(0.5f * cam.right);
    }
    else if (event.key.keysym.sym == SDLK_UP)
    {
      cout << "UP" << endl;
      cam.move(-0.5f * cam.forward);
    }
    else if (event.key.keysym.sym == SDLK_DOWN)
    {
      cout << "DOWN" << endl;
      cam.move(0.5f * cam.forward);
    }
    else if (event.key.keysym.sym == SDLK_LSHIFT)
    {
      cout << "LSHIFT" << endl;
      cam.move(-0.5f * cam.up);
    }
    else if (event.key.keysym.sym == SDLK_SPACE)
    {
      cout << "SPACE" << endl;
      cam.move(0.5f * cam.up);
    }
  }
}

vector<float> Interpolate(float a, float b, int n)
{
  vector<float> result;
  if (n == 1)
  {
    result.push_back((a + b) / 2);
    return result;
  }
  else if (n > 1)
  {
    float step = (b - a) / (n - 1);
    for (int i = 0; i < n; ++i)
    {
      result.push_back(a);
      a += step;
    }
  }
  return result;
}

vector<vec3> Interpolate(vec3 a, vec3 b, int n)
{
  vector<vec3> result;
  if (n == 1)
  {
    result.push_back((a + b) / 2.0f);
    return result;
  }
  else if (n > 1)
  {
    vec3 step = (b - a) / (n - 1.0f);
    for (int i = 0; i < n; ++i)
    {
      result.push_back(a);
      a += step;
    }
  }
  return result;
}

vector<CanvasPoint> Interpolate(CanvasPoint a, CanvasPoint b, int n)
{
  vector<CanvasPoint> result = vector<CanvasPoint>();
  if (n == 1)
  {
    result.push_back((a + b) / 2.0f);
    return result;
  }
  else if (n > 1)
  {
    CanvasPoint step = (b - a) / (n - 1.0f);
    for (int i = 0; i < n; ++i)
    {
      result.push_back(a);
      a += step;
    }
  }
  return result;
}

void line(CanvasPoint p, CanvasPoint q, int colour)
{
  float x_diff = p.x - q.x;
  float y_diff = p.y - q.y;
  float max_diff = std::max(abs(x_diff), abs(y_diff));
  vector<float> interpolate_x = Interpolate(p.x, q.x, max_diff + 1);
  vector<float> interpolate_y = Interpolate(p.y, q.y, max_diff + 1);
  for (int i = 0; i < max_diff; i++)
  {
    window.setPixelColour(round(interpolate_x[i]), round(interpolate_y[i]),
                          colour);
  }
}

void triangle(CanvasTriangle t, int colour, bool filled)
{
  if (filled)
  {
    sort(begin(t.vertices), end(t.vertices));
    // if (t.vertices[0].x < 0 || t.vertices[0].x >= WIDTH ||
    //     t.vertices[0].y < 0 || t.vertices[0].y >= HEIGHT ||
    //     t.vertices[1].x < 0 || t.vertices[1].x >= WIDTH ||
    //     t.vertices[1].y < 0 || t.vertices[1].y >= HEIGHT ||
    //     t.vertices[2].x < 0 || t.vertices[2].x >= WIDTH ||
    //     t.vertices[2].y < 0 || t.vertices[2].y >= HEIGHT) return;
    vector<CanvasPoint> line1 =
        Interpolate(t.vertices[0], t.vertices[1],
                    t.vertices[1].y - t.vertices[0].y + 1);
    vector<CanvasPoint> line2 =
        Interpolate(t.vertices[0], t.vertices[2],
                    t.vertices[2].y - t.vertices[0].y + 1);
    vector<CanvasPoint> line3 =
        Interpolate(t.vertices[1], t.vertices[2],
                    t.vertices[2].y - t.vertices[1].y + 1);

    for (int y = t.vertices[0].y; y < t.vertices[1].y; y++)
    {
      int xStart, xEnd;
      vector<CanvasPoint> currentLine;
      if (line1[y - t.vertices[0].y].x < line2[y - t.vertices[0].y].x)
      {
        xStart = line1[y - t.vertices[0].y].x;
        xEnd = line2[y - t.vertices[0].y].x;
        currentLine = Interpolate(line1[y - t.vertices[0].y], line2[y - t.vertices[0].y], xEnd - xStart);
      }
      else
      {
        xStart = line2[y - t.vertices[0].y].x;
        xEnd = line1[y - t.vertices[0].y].x;
        currentLine = Interpolate(line2[y - t.vertices[0].y], line1[y - t.vertices[0].y], xEnd - xStart);
      }
      for (int x = xStart; x < xEnd; x++)
      {
        if (x >= 0 && y >= 0 && x < WIDTH && y < HEIGHT)
        {
          if (currentLine[x - xStart].depth < depthBuffer[x + y * WIDTH])
          {
            window.setPixelColour(x, y, colour);
            depthBuffer[x + y * WIDTH] = currentLine[x - xStart].depth;
          }
        }
      }
    }
    for (int y = t.vertices[1].y; y < t.vertices[2].y; y++)
    {
      int xStart, xEnd;
      vector<CanvasPoint> currentLine;
      if (line3[y - t.vertices[1].y].x < line2[y - t.vertices[0].y].x)
      {
        xStart = line3[y - t.vertices[1].y].x;
        xEnd = line2[y - t.vertices[0].y].x;
        currentLine = Interpolate(line3[y - t.vertices[1].y], line2[y - t.vertices[0].y], xEnd - xStart);
      }
      else
      {
        xStart = line2[y - t.vertices[0].y].x;
        xEnd = line3[y - t.vertices[1].y].x;
        currentLine = Interpolate(line2[y - t.vertices[0].y], line3[y - t.vertices[1].y], xEnd - xStart);
      }
      for (int x = xStart; x < xEnd; x++)
      {
        if (x >= 0 && y >= 0 && x < WIDTH && y < HEIGHT)
        {
          if (currentLine[x - xStart].depth < depthBuffer[x + y * WIDTH])
          {
            window.setPixelColour(x, y, colour);
            depthBuffer[x + y * WIDTH] = currentLine[x - xStart].depth;
          }
        }
      }
    }
  }
  else
  {
    line(t.vertices[0], t.vertices[1], colour);
    line(t.vertices[1], t.vertices[2], colour);
    line(t.vertices[2], t.vertices[0], colour);
  }
}

void texturedTriangle(CanvasTriangle screenTri, CanvasTriangle texTri,
                      Texture tex)
{
  if (screenTri.vertices[0].y > screenTri.vertices[1].y)
  {
    swap(screenTri.vertices[0], screenTri.vertices[1]);
    swap(texTri.vertices[0], texTri.vertices[1]);
  }
  if (screenTri.vertices[1].y > screenTri.vertices[2].y)
  {
    swap(screenTri.vertices[1], screenTri.vertices[2]);
    swap(texTri.vertices[1], texTri.vertices[2]);
  }
  if (screenTri.vertices[0].y > screenTri.vertices[1].y)
  {
    swap(screenTri.vertices[0], screenTri.vertices[1]);
    swap(texTri.vertices[0], texTri.vertices[1]);
  }

  vector<float> sline1 =
      Interpolate(screenTri.vertices[0].x, screenTri.vertices[1].x,
                  abs(screenTri.vertices[0].y - screenTri.vertices[1].y) + 1);
  vector<float> sline2 =
      Interpolate(screenTri.vertices[0].x, screenTri.vertices[2].x,
                  abs(screenTri.vertices[0].y - screenTri.vertices[2].y) + 1);
  vector<float> sline3 =
      Interpolate(screenTri.vertices[1].x, screenTri.vertices[2].x,
                  abs(screenTri.vertices[1].y - screenTri.vertices[2].y) + 1);
  vector<float> txline1 =
      Interpolate(texTri.vertices[0].x, texTri.vertices[1].x,
                  abs(screenTri.vertices[0].y - screenTri.vertices[1].y) + 1);
  vector<float> txline2 =
      Interpolate(texTri.vertices[0].x, texTri.vertices[2].x,
                  abs(screenTri.vertices[0].y - screenTri.vertices[2].y) + 1);
  vector<float> txline3 =
      Interpolate(texTri.vertices[1].x, texTri.vertices[2].x,
                  abs(screenTri.vertices[1].y - screenTri.vertices[2].y) + 1);
  vector<float> tyline1 =
      Interpolate(texTri.vertices[0].y, texTri.vertices[1].y,
                  abs(screenTri.vertices[0].y - screenTri.vertices[1].y) + 1);
  vector<float> tyline2 =
      Interpolate(texTri.vertices[0].y, texTri.vertices[2].y,
                  abs(screenTri.vertices[0].y - screenTri.vertices[2].y) + 1);
  vector<float> tyline3 =
      Interpolate(texTri.vertices[1].y, texTri.vertices[2].y,
                  abs(screenTri.vertices[1].y - screenTri.vertices[2].y) + 1);

  for (int y = screenTri.vertices[0].y; y < screenTri.vertices[1].y; y++)
  {
    float tx1 = txline1[y - screenTri.vertices[0].y];
    float ty1 = tyline1[y - screenTri.vertices[0].y];
    float tx2 = txline2[y - screenTri.vertices[0].y];
    float ty2 = tyline2[y - screenTri.vertices[0].y];
    int xStart = std::min(sline1[y - screenTri.vertices[0].y],
                          sline2[y - screenTri.vertices[0].y]);
    int xEnd = std::max(sline1[y - screenTri.vertices[0].y],
                        sline2[y - screenTri.vertices[0].y]);
    vector<float> txs = Interpolate(tx1, tx2, xEnd - xStart);
    vector<float> tys = Interpolate(ty1, ty2, xEnd - xStart);
    for (int x = xStart; x < xEnd; x++)
    {
      window.setPixelColour(
          x, y,
          tex.buff[(int)tys[x - xStart] * tex.width + (int)txs[x - xStart]]);
    }
  }
  for (int y = screenTri.vertices[1].y; y < screenTri.vertices[2].y; y++)
  {
    float tx1 = txline2[y - screenTri.vertices[0].y];
    float ty1 = tyline2[y - screenTri.vertices[0].y];
    float tx2 = txline3[y - screenTri.vertices[1].y];
    float ty2 = tyline3[y - screenTri.vertices[1].y];
    if (tx1 > tx2)
    {
      swap(tx1, tx2);
      swap(ty1, ty2);
    }
    int xStart = std::min(sline3[y - screenTri.vertices[1].y],
                          sline2[y - screenTri.vertices[0].y]);
    int xEnd = std::max(sline3[y - screenTri.vertices[1].y],
                        sline2[y - screenTri.vertices[0].y]);
    vector<float> txs = Interpolate(tx1, tx2, xEnd - xStart);
    vector<float> tys = Interpolate(ty1, ty2, xEnd - xStart);
    for (int x = xStart; x < xEnd; x++)
    {
      window.setPixelColour(
          x, y,
          tex.buff[(int)tys[x - xStart] * tex.width + (int)txs[x - xStart]]);
    }
  }
}

int *loadPPM(string fileName, int &width, int &height)
{
  ifstream f;
  string s;
  f.open(fileName, ios::in | ios::binary);
  f >> s;
  skipHashWS(f);
  f >> s;
  width = stoi(s);
  skipHashWS(f);
  f >> s;
  height = stoi(s);
  skipHashWS(f);
  f >> s;
  f.seekg(1, f.cur);

  int *buff = new int[width * height];
  for (int i = 0; i < width * height; i++)
  {
    buff[i] = 0xff000000;
    f.read((char *)&buff[i], 3);
  }
  return buff;
}

void savePPM(string fileName, DrawingWindow *window)
{
  ofstream f;
  f.open(fileName, ios::out | ios::binary);
  f << "P6" << endl;
  f << window->width << " " << window->height << endl;
  f << 255 << endl;
  for (int y = 0; y < window->height; y++)
  {
    unsigned char *buffer = new unsigned char[window->width * 3];
    for (int x = 0; x < window->width; x++)
    {
      uint32_t value = window->getPixelColour(x, y);
      buffer[x * 3] = value & 0xff;
      buffer[x * 3 + 1] = (value & 0xff00) >> 8;
      buffer[x * 3 + 2] = (value & 0xff0000) >> 16;
    }
    f.write((char *)buffer, window->width * 3);
  }
  f.close();
}

void skipHashWS(ifstream &f)
{
  ws(f);
  char current;
  f.read(&current, 1);
  if (current == '#')
  {
    f.ignore(1000, '\n');
    ws(f);
  }
  else
  {
    f.seekg(-1, f.cur);
  }
}