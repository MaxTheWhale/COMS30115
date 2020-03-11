#include <ModelTriangle.h>
#include <CanvasTriangle.h>
#include <DrawingWindow.h>
#include <Utils.h>
#include <fstream>
#include <glm/glm.hpp>
#include <unordered_map>
#include <vector>
#include <sys/time.h>

using namespace std;
using namespace glm;

#define WIDTH 320
#define HEIGHT 240

unsigned long long getTime();
void draw();
void line(CanvasPoint p, CanvasPoint q, int colour);
void triangle(CanvasTriangle t, int colour, bool filled = false);
int *loadPPM(string fileName, int &width, int &height);
void savePPM(string fileName, DrawingWindow *window);
void skipHashWS(ifstream &f);
void update(mat4 &cam);
void handleEvent(SDL_Event event, mat4 &camToWorld);
float depthBuffer[WIDTH * HEIGHT];
vector<float> Interpolate(float a, float b, int n);
vector<vec3> Interpolate(vec3 a, vec3 b, int n);
mat4 buildProjection(float fov, float near, float far);
mat4 lookAt(const vec3 &from, const vec3 &to);
unordered_map<string, Colour> loadMTL(string fileName);
vector<ModelTriangle> loadOBJ(string fileName,
                              unordered_map<string, Colour> palette);
vector<ModelTriangle> loadOBJ(string fileName,
                              unordered_map<string, Colour> palette)
{
  ifstream f;
  string s;
  Colour colour;
  vector<vec3> vertices;
  vector<ModelTriangle> faces;
  f.open(fileName, ios::in);
  while (!f.eof())
  {
    if (f >> s)
    {
      if (s == "mtllib")
      {
        f >> s;
      }
      if (s == "o")
      {
        f >> s;
      }
      if (s == "usemtl")
      {
        f >> s;
        colour = palette[s];
      }
      if (s == "v")
      {
        float x, y, z;
        f >> x >> y >> z;
        vertices.push_back(vec3(x, y, z));
      }
      if (s == "f")
      {
        string a, b, c;
        f >> a >> b >> c;
        faces.push_back(ModelTriangle(vertices[stoi(split(a, '/')[0]) - 1],
                                      vertices[stoi(split(b, '/')[0]) - 1],
                                      vertices[stoi(split(c, '/')[0]) - 1],
                                      colour));
      }
    }
  }
  return faces;
}

unordered_map<string, Colour> loadMTL(string fileName)
{
  unordered_map<string, Colour> palette;

  ifstream f;
  string s;
  f.open(fileName, ios::in);
  while (!f.eof())
  {
    f >> s;
    if (s == "newmtl")
    {
      string key, r, g, b;
      f >> key;
      f >> s;
      f >> r;
      f >> g;
      f >> b;
      palette[key] = Colour(key, stof(r) * 255, stof(g) * 255, stof(b) * 255);
    }
  }
  return palette;
}

std::ostream &operator<<(std::ostream &os, const std::vector<float> vector)
{
  os << '[';
  for (int i = 0; i < (int)vector.size(); i++)
  {
    os << vector[i];
    if (i != (int)vector.size() - 1)
    {
      os << ',';
    }
  }
  os << ']';
  return os;
}
std::ostream &operator<<(std::ostream &os, const vec3 vec3)
{
  std::vector<float> vec{vec3.x, vec3.y, vec3.z};
  os << vec;
  return os;
}
std::ostream &operator<<(std::ostream &os, const vec4 vec4)
{
  std::vector<float> vec{vec4.x, vec4.y, vec4.z, vec4.w};
  os << vec;
  return os;
}
std::ostream &operator<<(std::ostream &os, const mat3 mat3)
{
  for (int i = 0; i < 3; i++)
  {
    os << mat3[i];
  }
  return os;
}
std::ostream &operator<<(std::ostream &os, const mat4 mat4)
{
  for (int i = 0; i < 4; i++)
  {
    os << mat4[i];
  }
  return os;
}

float mat4Dist(mat4 &a, mat4 &b)
{
  float result = 0;
  for (int i = 0; i < 4; i++)
  {
    result += distance(a[i], b[i]);
  }
  return result;
}

class Model
{
public:
  vector<ModelTriangle> tris;
  mat4 transform;
  unordered_map<string, Colour> palette;
  Model(string filename)
  {
    palette = loadMTL(filename + ".mtl");
    tris = loadOBJ(filename + ".obj", palette);
    transform = lookAt(vec3(0, 0, 0), vec3(0, 0, -1));
  }
};

long long lastFrameTime = getTime();

long long frameCount = 0;

//returns milliseconds since epoch
unsigned long long getTime()
{
  struct timeval tv;

  gettimeofday(&tv, NULL);

  unsigned long long millisecondsSinceEpoch =
      (unsigned long long)(tv.tv_sec) * 1000 +
      (unsigned long long)(tv.tv_usec) / 1000;
  return millisecondsSinceEpoch;
}

float deltaTime()
{
  long long ms = getTime() - lastFrameTime;
  return (float)ms / (float)1000;
}

// x = Xf/Z
// y = Yf/Z
// x, y = 2d screen coords
// X, Y, Z = 3d world coords
// f = focal length

void drawTriangles(Model model, mat4 camToWorld, mat4 projection)
{
  mat4 posMat = mat4(vec4(1, 0, 0, 0), vec4(0, 1, 0, 0), vec4(0, 0, 1, 0), camToWorld[3]);
  camToWorld[3] = vec4(0, 0, 0, 1);
  for (auto tri : model.tris)
  {
    vec4 camToV1 = posMat * model.transform * tri.vertices[0];
    camToV1 = camToWorld * camToV1;
    camToV1 = projection * camToV1;
    camToV1.x /= camToV1.w;
    camToV1.y /= camToV1.w;
    camToV1.z /= camToV1.w;
    vec4 camToV2 = posMat * model.transform * tri.vertices[1];
    camToV2 = camToWorld * camToV2;
    camToV2 = projection * camToV2;
    camToV2.x /= camToV2.w;
    camToV2.y /= camToV2.w;
    camToV2.z /= camToV2.w;
    vec4 camToV3 = posMat * model.transform * tri.vertices[2];
    camToV3 = camToWorld * camToV3;
    camToV3 = projection * camToV3;
    camToV3.x /= camToV3.w;
    camToV3.y /= camToV3.w;
    camToV3.z /= camToV3.w;
    CanvasPoint v1 = CanvasPoint(
        (camToV1.x + 1) * 0.5f * WIDTH,
        (1 - (camToV1.y + 1) * 0.5f) * HEIGHT,
        1.0 / camToV1.z);
    CanvasPoint v2 = CanvasPoint(
        (camToV2.x + 1) * 0.5f * WIDTH,
        (1 - (camToV2.y + 1) * 0.5f) * HEIGHT,
        1.0 / camToV2.z);
    CanvasPoint v3 = CanvasPoint(
        (camToV3.x + 1) * 0.5f * WIDTH,
        (1 - (camToV3.y + 1) * 0.5f) * HEIGHT,
        1.0 / camToV3.z);
    if (v1.depth > 0 || v2.depth > 0 || v3.depth > 0)
      triangle(CanvasTriangle(v1, v2, v3), tri.colour.toPackedInt(), true);
  }
}

mat4 lookAt(const vec3 &from, const vec3 &to)
{
  vec3 forward = normalize(from - to);
  vec3 right = normalize(cross(vec3(0, 1, 0), forward));
  vec3 up = normalize(cross(forward, right));
  cout << forward.x << ' ' << forward.y << ' ' << forward.z << '\n';
  cout << up.x << ' ' << up.y << ' ' << up.z << '\n';
  cout << right.x << ' ' << right.y << ' ' << right.z << '\n';
  return mat4(vec4(right.x, up.x, forward.x, 0),
              vec4(right.y, up.y, forward.y, 0),
              vec4(right.z, up.z, forward.z, 0),
              vec4(-from.x, -from.y, -from.z, 1));
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
  // Texture img = Texture("texture.ppm");
  // unordered_map<string, Colour> palette = loadMTL("cornell-box.mtl");
  // vector<ModelTriangle> tris = loadOBJ("cornell-box.obj", palette);

  Model cornell = Model("cornell-box");

  for (int i = 0; i < WIDTH * HEIGHT; i++)
  {
    depthBuffer[i] = 0;
  }

  mat4 camToWorld = lookAt(vec3(5.0f, 2.5f, 3.0f), vec3(0.0f, 2.5f, -3.0f));
  cameraTransforms.push_back(lookAt(vec3(5.0f, 2.5f, 3.0f), vec3(0.0f, 2.5f, -3.0f)));
  cameraTransforms.push_back(lookAt(vec3(5.0f, 2.5f, 3.0f), vec3(0.0f, 4.0f, -3.0f)));
  cameraTransforms.push_back(lookAt(vec3(5.0f, 2.5f, 0.0f), vec3(0.0f, 3.0f, 0.0f)));
  cameraTransforms.push_back(lookAt(vec3(0.0f, 2.5f, 3.0f), vec3(0.0f, 3.0f, -3.0f)));
  cameraTransforms.push_back(lookAt(vec3(0.0f, 2.5f, 4.0f), vec3(0.0f, 4.0f, 0.0f)));
  cameraTransforms.push_back(lookAt(vec3(0.0f, 2.5f, 4.0f), vec3(0.0f, 2.0f, 0.0f)));
  mat4 projection = buildProjection(90.0f, 0.1f, 100.0f);
  cout << camToWorld[3].x << ' ' << camToWorld[3].y << ' ' << camToWorld[3].z << ' ' << camToWorld[3].w << '\n';

  while (true)
  {
    // We MUST poll for events - otherwise the window will freeze !
    if (window.pollForInputEvents(&event))
      handleEvent(event, camToWorld);
    if (frameCount % 60 == 0)
    {
      cornell.transform[3].y++;
    }
    update(camToWorld);
    lastFrameTime = getTime();
    draw();
    drawTriangles(cornell, camToWorld, projection);
    // Need to render the frame at the end, or nothing actually gets shown on
    // the screen !
    window.renderFrame();
    frameCount++;
  }
}

void draw()
{
  window.clearPixels();
  for (int i = 0; i < WIDTH * HEIGHT; i++)
  {
    depthBuffer[i] = 0;
  }
}

int moveStage = 0;

void update(mat4 &cam)
{
  // Function for performing animation (shifting artifacts or moving the camera)
  float dt = deltaTime();

  if (moveStage != -1)
  {
    int prevStage = moveStage > 0 ? moveStage - 1 : 0;
    //vec3 camPos = -vec3(cam[3].x, cam[3].y, cam[3].z);
    //vec3 delta = dt * (cameraPositions[prevStage] - cameraPositions[moveStage]);
    mat4 delta = -dt * (cameraTransforms[prevStage] - cameraTransforms[moveStage]);
    //cout << camPos << endl;
    mat4 newTransform = mat4();
    for (int i = 0; i < 4; i++)
    {
      newTransform[i] = cam[i] + delta[i];
    }
    //would the move take us further from our goal
    float currentDist = mat4Dist(cam, cameraTransforms[moveStage]);
    float newDist = mat4Dist(newTransform, cameraTransforms[moveStage]);
    cout << "currentDist = " << currentDist << " newDist = " << newDist << endl;
    if (currentDist <= newDist)
    {
      //cam[3] = vec4(cameraPositions[moveStage].x, cameraPositions[moveStage].y, cameraPositions[moveStage].z, cam[3].w);
      if (moveStage + 1 < (int)cameraTransforms.size())
      {
        moveStage++;
      }
      else
      {
        moveStage = -1;
        savePPM("window.ppm", &window);
      }
      cout << "move stage is now " << moveStage << endl;
    }
    else
    {
      //cam[3] += vec4(delta.x, delta.y, delta.z, 0);
      for (int i = 0; i < 4; i++)
      {
        cam[i] += delta[i];
      }
    }
  }

  //cam[3] += dt * transpose(cam)[0] * 2.0f;
  //cam[3] += dt * vec4(1, 0, 0, 0);
  //cam = lookAt(, vec3(0.0f, 2.5f, -3.0f));
}

void handleEvent(SDL_Event event, mat4 &camToWorld)
{
  if (event.type == SDL_KEYDOWN)
  {
    if (event.key.keysym.sym == SDLK_LEFT)
    {
      cout << "LEFT" << endl;
      camToWorld[3] += (0.5f * transpose(camToWorld)[0]);
    }
    else if (event.key.keysym.sym == SDLK_RIGHT)
    {
      cout << "RIGHT" << endl;
      camToWorld[3] -= (0.5f * transpose(camToWorld)[0]);
    }
    else if (event.key.keysym.sym == SDLK_UP)
    {
      cout << "UP" << endl;
      camToWorld[3] += (0.5f * transpose(camToWorld)[2]);
    }
    else if (event.key.keysym.sym == SDLK_DOWN)
    {
      cout << "DOWN" << endl;
      camToWorld[3] -= (0.5f * transpose(camToWorld)[2]);
    }
    else if (event.key.keysym.sym == SDLK_LSHIFT)
    {
      cout << "LSHIFT" << endl;
      camToWorld[3] += (0.5f * transpose(camToWorld)[1]);
    }
    else if (event.key.keysym.sym == SDLK_SPACE)
    {
      cout << "SPACE" << endl;
      camToWorld[3] -= (0.5f * transpose(camToWorld)[1]);
    }
    cout << -camToWorld[3].x << ' ' << -camToWorld[3].y << ' ' << -camToWorld[3].z << ' ' << camToWorld[3].w << '\n';
  }
}

mat4 buildProjection(float fov, float near, float far)
{
  float aspect_ratio = WIDTH / (float)HEIGHT;
  float top = tan(radians(fov / 2)) * near;
  float bottom = -top;
  float right = top * aspect_ratio;
  float left = bottom;
  return mat4(vec4(2 * near / (right - left), 0, 0, 0),
              vec4(0, 2 * near / (top - bottom), 0, 0),
              vec4((right + left) / (right - left), (top + bottom) / (top - bottom), -(far + near) / (far - near), -1),
              vec4(0, 0, -(2 * far * near) / (far - near), 0));
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
          if (currentLine[x - xStart].depth > depthBuffer[x + y * WIDTH])
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
          if (currentLine[x - xStart].depth > depthBuffer[x + y * WIDTH])
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