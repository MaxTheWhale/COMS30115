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

#define WIDTH 640
#define HEIGHT 480
#define MOUSE_SENSITIVITY 0.0015f

void draw();
void line(CanvasPoint p, CanvasPoint q, int colour);
void triangle(CanvasTriangle t, int colour, bool filled = false);
int *loadPPM(string fileName, int &width, int &height);
void savePPM(string fileName, DrawingWindow *window);
void skipHashWS(ifstream &f);
void update(Camera &cam, vector<Model> models);
void handleEvent(SDL_Event event, Camera &cam);
float depthBuffer[WIDTH * HEIGHT];
bool wireframe;
vector<float> Interpolate(float a, float b, int n);
vector<vec3> Interpolate(vec3 a, vec3 b, int n);
vector<vec4> Interpolate(vec4 a, vec4 b, int n);

bool toRaytrace = false;
bool softShadows = false;

bool inClipSpace(vec3 point) {
  return (point.x > -1.0f && point.x < 1.0f && point.y > -1.0f && point.y < 1.0f && point.z > -1.0f && point.z < 1.0f);
}

// based on https://casual-effects.com/research/McGuire2011Clipping/McGuire-Clipping.pdf
int clipTriangle(vector<ModelTriangle>& tris, const vec4& normal) {
  vec4 temp;
  vector<int> toBeCulled;
  int n = tris.size();
  for (int i = 0; i < n; i++) {
    vector<float> distances;
    distances.push_back(dot(tris[i].vertices[0], normal) + 1.0f);
    distances.push_back(dot(tris[i].vertices[1], normal) + 1.0f);
    distances.push_back(dot(tris[i].vertices[2], normal) + 1.0f);
    if (distances[0] >= 0.0f && distances[1] >= 0.0f && distances[2] >= 0.0f) {
      continue;
    }
    if (distances[0] < 0.0f && distances[1] < 0.0f && distances[2] < 0.0f) {
      toBeCulled.push_back(i);
    }
    bool nextInside;
    if (distances[1] >= 0.0f && distances[0] < 0.0f) {
      nextInside = (distances[2] >= 0.0f);
      temp = tris[i].vertices[0];
      tris[i].vertices[0] = tris[i].vertices[1];
      tris[i].vertices[1] = tris[i].vertices[2];
      tris[i].vertices[2] = temp;
      rotate(distances.begin(),distances.begin()+1,distances.end());
    }
    else if (distances[2] >= 0.0f && distances[1] < 0.0f) {
      nextInside = (distances[0] >= 0.0f);
      temp = tris[i].vertices[2];
      tris[i].vertices[2] = tris[i].vertices[1];
      tris[i].vertices[1] = tris[i].vertices[0];
      tris[i].vertices[0] = temp;
      rotate(distances.begin(),distances.begin()+2,distances.end());
    }
    else {
      nextInside = (distances[1] >= 0.0f);
    }
    temp = (tris[i].vertices[0] + (tris[i].vertices[2] - tris[i].vertices[0]) * (distances[0] / (distances[0] - distances[2])));
    if (nextInside) {
      tris[i].vertices[2] = (tris[i].vertices[1] + (tris[i].vertices[2] - tris[i].vertices[1]) * (distances[1] / (distances[1] - distances[2])));
      tris.push_back(ModelTriangle(tris[i].vertices[0], tris[i].vertices[2], temp, tris[i].colour));
    }
    else {
      tris[i].vertices[1] = (tris[i].vertices[0] + (tris[i].vertices[1] - tris[i].vertices[0]) * (distances[0] / (distances[0] - distances[1])));
      tris[i].vertices[2] = temp;
    }
  }
  for (auto i : toBeCulled) {
    tris.erase(tris.begin() + i);
  }
  return tris.size();
}

int clipToView(vector<ModelTriangle>& tris) {
  const vec4 normals[6] = {vec4(1, 0, 0, 0), vec4(0, 1, 0, 0), vec4(0, 0, 1, 0), vec4(-1, 0, 0, 0), vec4(0, -1, 0, 0), vec4(0, 0, -1, 0)};
  for (auto n : normals) {
    int num = clipTriangle(tris, n);
    if (num == 0) {
      return 0;
    }
  }
  return tris.size();
}

void drawTriangles(Model &model, Camera &cam)
{
  mat4 MVP = cam.projection * cam.worldToCamera() * model.transform;
  for (auto tri : model.tris)
  {
    tri.vertices[0] = MVP * tri.vertices[0];
    tri.vertices[0] /= tri.vertices[0].w;
    tri.vertices[1] = MVP * tri.vertices[1];
    tri.vertices[1] /= tri.vertices[1].w;
    tri.vertices[2] = MVP * tri.vertices[2];
    tri.vertices[2] /= tri.vertices[2].w;
    
    vector<ModelTriangle> clippedTris;
    clippedTris.push_back(tri);
    clipToView(clippedTris);

    for (auto t : clippedTris) {
      CanvasPoint v1 = CanvasPoint(
        (t.vertices[0].x + 1.0f) * 0.5f * WIDTH,
        (1 - (t.vertices[0].y + 1.0f) * 0.5f) * HEIGHT,
        (99.9f / 2) * t.vertices[0].z + (100.1f / 2));
      CanvasPoint v2 = CanvasPoint(
        (t.vertices[1].x + 1.0f) * 0.5f * WIDTH,
        (1 - (t.vertices[1].y + 1.0f) * 0.5f) * HEIGHT,
        (99.9f / 2) * t.vertices[1].z + (100.1f / 2));
      CanvasPoint v3 = CanvasPoint(
        (t.vertices[2].x + 1.0f) * 0.5f * WIDTH,
        (1 - (t.vertices[2].y + 1.0f) * 0.5f) * HEIGHT,
        (99.9f / 2) * t.vertices[2].z + (100.1f / 2));
      triangle(CanvasTriangle(v1, v2, v3), tri.colour.toPackedInt(), wireframe);
    }
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


void handleMouse(Camera& cam) {
  int motion_x = 0;
  int motion_y = 0;
  SDL_GetRelativeMouseState(&motion_x, &motion_y);
  if (motion_x || motion_y) {
    cam.rotate(vec3(motion_y * MOUSE_SENSITIVITY, motion_x * MOUSE_SENSITIVITY, 0));
  }
}

int darkenColour(Colour colour, float brightness, float specular) {
  colour.red *= brightness;
  colour.green *= brightness;
  colour.blue *= brightness;

  colour.red += specular * 255;
  colour.blue += specular * 255;
  colour.green += specular * 255;

  colour.red  = colour.red > 255 ? 255 : colour.red;
  colour.green  = colour.green > 255 ? 255 : colour.green;
  colour.blue  = colour.blue > 255 ? 255 : colour.blue;

  return colour.toPackedInt();
}

float vectorLength(vec4 v) {
  return sqrt(pow(v.x, 2) + pow(v.y, 2) + pow(v.z, 2) + pow(v.w, 2));
}

vec3 toThree(vec4 v) {
  return vec3(v.x, v.y, v.z);
}

void texturedTriangle(CanvasTriangle screenTri, CanvasTriangle texTri,
                      Texture tex);

DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);

vector<ModelTriangle> getLights(Model model) {
  vector<ModelTriangle> lights = vector<ModelTriangle>();
  for(ModelTriangle triangle : model.tris) {
    if(triangle.name == "light") lights.push_back(triangle);
  }

  return lights;
}

void raytrace(Camera camera, Model model, int softness) {
  float fov = 90;
  float aspectRatio = WIDTH/(float)HEIGHT;
  float angle = tan(0.5 * fov * M_PI / 180.0); // just fov*0.5 converted to radians

  for(int j = 0; j < HEIGHT; j++) {
    for(int i = 0; i < WIDTH; i++) {
      vec2 NDC = vec2((i + 0.5) * (1 / (float) WIDTH), (j + 0.5) * (1 / (float) HEIGHT));
      float x = (2 * (NDC.x) - 1) * angle * aspectRatio;
      float y = (1 - 2 * (NDC.y)) * angle;

      vec4 rayDirection = camera.transform * vec4(x, y, -1, 0);
      float minDistance = std::numeric_limits<float>::infinity();
      ModelTriangle intersection = ModelTriangle();
      bool foundIntersection = false;
      float t = 0;

      for(ModelTriangle triangle : model.tris) {
        vec4 e0 = triangle.vertices[1] - triangle.vertices[0];
        vec4 e1 = triangle.vertices[2] - triangle.vertices[0];
        vec4 SPVector = camera.transform[3] - triangle.vertices[0];
        mat4 DEMatrix(-rayDirection, e0, e1, vec4(1,1,1,1));
        vec4 possibleSolution = glm::inverse(DEMatrix) * SPVector;

        //check if ray intersects triangle and not just triangle plane
        if(possibleSolution.y >= 0 && possibleSolution.y <= 1 && possibleSolution.z >= 0 && possibleSolution.z <= 1 && possibleSolution.y + possibleSolution.z <= 1) {
          if(possibleSolution.x < minDistance) {
            foundIntersection = true;
            t = possibleSolution.x;
            intersection = triangle;
            minDistance = possibleSolution.x;
          }
        }
      }

      if(foundIntersection) {
        float ambience = 0.25f;
        float intensity = 1.0f / softness;
        float shadow = 0.1f / softness;
        float shadowCount = 1.0f;
        float brightnessCount = 0.0f;
        float angleCount = 0.0f;
        float specularCount = 0.0f;
        bool inShadow = false;


        vector<ModelTriangle> lights = getLights(model);
        vec4 intersectionPoint = camera.transform[3] + (t * rayDirection);
        intersectionPoint.w = 1;

        for(ModelTriangle light : lights) {
          vector<vec4> points = vector<vec4>();
          vector<vec4> a = Interpolate(light.vertices[0], light.vertices[1], softness);
          vector<vec4> b = Interpolate(light.vertices[0], light.vertices[2], softness);
          vector<vec4> c = Interpolate(light.vertices[1], light.vertices[2], softness);
          points.insert(points.end(), a.begin(), a.end());
          points.insert(points.end(), b.begin(), b.end());
          points.insert(points.end(), c.begin(), c.end());

          for(vec4 vertex : points) {
            vec4 shadowRayDirection = vertex - intersectionPoint;

            //cout << shadowRayDirection << " with length " << vectorLength(shadowRayDirection) << endl;

            vec3 intersectionNormal = glm::cross(toThree(intersection.vertices[1] - intersection.vertices[0]), toThree(intersection.vertices[2] - intersection.vertices[0]));
            float angleOfIncidence = glm::dot(toThree(glm::normalize(shadowRayDirection)), glm::normalize(intersectionNormal));
            angleOfIncidence = angleOfIncidence < 0 ? ambience : angleOfIncidence;

            float tolerance = 0.01f;
            float shadowBias = 1e-4;

            for(ModelTriangle triangle : model.tris) {
              vec4 e0 = triangle.vertices[1] - triangle.vertices[0];
              vec4 e1 = triangle.vertices[2] - triangle.vertices[0];
              vec4 SPVector = (intersectionPoint + shadowBias) - triangle.vertices[0];
              mat4 DEMatrix(-shadowRayDirection, e0, e1, vec4(1,1,1,1));
              vec4 possibleSolution = glm::inverse(DEMatrix) * SPVector;


              //check if ray intersects triangle and not just triangle plane
              if(possibleSolution.y >= 0 && possibleSolution.y <= 1 && possibleSolution.z >= 0 && possibleSolution.z <= 1 && possibleSolution.y + possibleSolution.z <= 1) {
                //cout << "intersection length " << possibleSolution.x << endl;
                if(possibleSolution.x < 1 - tolerance && possibleSolution.x > tolerance) {
                  //cout << "in shadow for " << i << "," << j << endl;
                  inShadow = true;
                }
              }
            }

            float brightness = intensity/pow(vectorLength(shadowRayDirection),2);
            brightness = brightness < 0 ? 0 : brightness;
            brightness = brightness > 1 ? 1 : brightness;
            shadowCount -= inShadow ? shadow : 0;
            brightnessCount += brightness;
            angleCount += angleOfIncidence;

            vec3 reflection = toThree(-shadowRayDirection) - 2 * (glm::dot(toThree(-shadowRayDirection), glm::normalize(intersectionNormal))) * glm::normalize(intersectionNormal);
            float specular = pow(glm::dot(glm::normalize(toThree(-rayDirection)), glm::normalize(reflection)), 128);

            specularCount += specular;

          }
        }

        shadowCount = shadowCount < 0 ? 0 : shadowCount;
        shadowCount = shadowCount > 1 ? 1 : shadowCount;
        brightnessCount = brightnessCount < 0 ? 0 : brightnessCount;
        brightnessCount = brightnessCount > 1 ? 1 : brightnessCount;
        angleCount = angleCount < 0 ? 0 : angleCount;
        angleCount = angleCount > 1 ? 1 : angleCount;
        specularCount = specularCount < 0 ? 0 : specularCount;
        specularCount = specularCount > 1 ? 1 : specularCount;

        if(intersection.name == "light") window.setPixelColour(i, j, intersection.colour.toPackedInt());
        else window.setPixelColour(i, j, inShadow ? darkenColour(intersection.colour, angleCount * shadowCount * brightnessCount, specularCount) : darkenColour(intersection.colour, angleCount * brightnessCount, specularCount));
      } else {
        window.setPixelColour(i, j, 0);
      }
    }
  }
}

//vector<vec3> cameraPositions{ vec3(5.0f, 2.5f, 3.0f), vec3(5.0f, 0.0f, 3.0f), vec3(5.0f, 0.0f, 6.0f) };
vector<mat4> cameraTransforms = vector<mat4>();
int main(int argc, char *argv[])
{
  SDL_Event event;
  //SDL_SetRelativeMouseMode(SDL_TRUE);

  Model cornell = Model("cornell-box");

  Camera cam;
  cam.setProjection(90.0f, WIDTH / (float)HEIGHT, 0.1f, 100.0f);
  // cam.lookAt(vec3(5.0f, 2.5f, 3.0f), vec3(0.0f, 2.5f, 0.0f));
  // cam.moves.push(Movement(cam.transform));
  // cam.moves.top().lookAt(cam.getPosition(), vec3(0, -2.5f, 0));

  Times::init();

  while (true)
  {
    //cout << "camera transform = " << cam.transform << endl;
    Times::update();
    //cout << "deltaTime: " << Times::deltaTime() << endl;
    // We MUST poll for events - otherwise the window will freeze !
    if (window.pollForInputEvents(&event))
      handleEvent(event, cam);
    //handleMouse(cam);
    //cout << "deltaTime = " << Times::deltaTime() << endl;
    update(cam, vector<Model>{cornell});

    if(toRaytrace) {
      raytrace(cam, cornell, softShadows ? 2 : 1);
    } else {
      draw();
      drawTriangles(cornell, cam);
    }
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
      cam.move(-0.5f * vec3(cam.transform[0].x, cam.transform[0].y, cam.transform[0].z));
    }
    else if (event.key.keysym.sym == SDLK_RIGHT)
    {
      cout << "RIGHT" << endl;
      cam.move(0.5f * vec3(cam.transform[0].x, cam.transform[0].y, cam.transform[0].z));
    }
    else if (event.key.keysym.sym == SDLK_UP)
    {
      cout << "UP" << endl;
      cam.move(-0.5f * vec3(cam.transform[2].x, cam.transform[2].y, cam.transform[2].z));
    }
    else if (event.key.keysym.sym == SDLK_DOWN)
    {
      cout << "DOWN" << endl;
      cam.move(0.5f * vec3(cam.transform[2].x, cam.transform[2].y, cam.transform[2].z));
    }
    else if (event.key.keysym.sym == SDLK_LSHIFT)
    {
      cout << "LSHIFT" << endl;
      cam.move(-0.5f * vec3(cam.transform[1].x, cam.transform[1].y, cam.transform[1].z));
    }
    else if (event.key.keysym.sym == SDLK_SPACE)
    {
      cout << "SPACE" << endl;
      cam.move(0.5f * vec3(cam.transform[1].x, cam.transform[1].y, cam.transform[1].z));
    }
    else if (event.key.keysym.sym == SDLK_w)
    {
      cout << "W" << endl;
      wireframe = !wireframe;
    }
    else if (event.key.keysym.sym == SDLK_r){
      toRaytrace = !toRaytrace;
      cout << "R = " << toRaytrace << endl;
    }
    else if (event.key.keysym.sym == SDLK_o) {
      cout << "O" << endl;
      cam.setPosition(vec3(0.0f,0.0f,0.0f));
    }
    else if (event.key.keysym.sym == SDLK_s) {
      softShadows = !softShadows;
      cout << "S = " << softShadows << endl;
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

vector<vec4> Interpolate(vec4 a, vec4 b, int n)
{
  vector<vec4> result;
  if (n == 1)
  {
    result.push_back((a + b) / 2.0f);
    return result;
  }
  else if (n > 1)
  {
    vec4 step = (b - a) / (n - 1.0f);
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

inline float edgeFunction(const CanvasPoint& v0, const CanvasPoint& v1, const CanvasPoint& p) {
  return (p.x - v0.x) * (v1.y - v0.y) - (p.y - v0.y) * (v1.x - v0.x);
}

void triangle(CanvasTriangle t, int colour, bool filled)
{
  if (filled)
  {
    int x_min = glm::min(t.vertices[0].x, t.vertices[1].x);
    x_min = glm::min((float)x_min, t.vertices[2].x);
    int x_max = glm::max(t.vertices[0].x, t.vertices[1].x);
    x_max = glm::max((float)x_max, t.vertices[2].x);
    int y_min = glm::min(t.vertices[0].y, t.vertices[1].y);
    y_min = glm::min((float)y_min, t.vertices[2].y);
    int y_max = glm::max(t.vertices[0].y, t.vertices[1].y);
    y_max = glm::max((float)y_max, t.vertices[2].y);
    if (x_min < 0) x_min = 0;
    if (y_min < 0) y_min = 0;
    if (x_max >= WIDTH) x_max = WIDTH - 1;
    if (y_max >= HEIGHT) y_max = HEIGHT - 1; 
    for (int y = y_min; y <= y_max; y++) {
      for (int x = x_min; x <= x_max; x++) {
        CanvasPoint p = CanvasPoint(x + 0.5f, y + 0.5f);
        float w0 = edgeFunction(t.vertices[1], t.vertices[2], p);
        float w1 = edgeFunction(t.vertices[2], t.vertices[0], p);
        float w2 = edgeFunction(t.vertices[0], t.vertices[1], p);
        if (w0 >= 0 && w1 >= 0 && w2 >= 0) {
          float area = edgeFunction(t.vertices[0], t.vertices[1], t.vertices[2]);
          w0 /= area;
          w1 /= area;
          w2 /= area;
          float depth = w0 * t.vertices[0].depth + w1 * t.vertices[1].depth + w2 * t.vertices[2].depth;
          if (depth < depthBuffer[y * WIDTH + x]) {
            depthBuffer[y * WIDTH + x] = depth;
            window.setPixelColour(x, y, colour);
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