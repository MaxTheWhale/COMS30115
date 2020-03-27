#include <ModelTriangle.h>
#include <CanvasTriangle.h>
#include <RayTriangleIntersection.h>
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
#include "Rigidbody.hpp"
using namespace std;
using namespace glm;

#define WIDTH 640
#define HEIGHT 480
#define SSAA true
#define SSAA_SCALE 4
#define S_WIDTH (WIDTH*SSAA_SCALE)
#define S_HEIGHT (HEIGHT*SSAA_SCALE)
#define MOUSE_SENSITIVITY 0.0015f
#define AMBIENCE 0.1f
#define ASPECT_RATIO WIDTH/(float)HEIGHT

enum CLIP_CODE {TOP = 1, RIGHT = 2, BOTTOM = 4, LEFT = 8};

void draw();
void line(CanvasPoint p, CanvasPoint q, int colour, uint32_t *buffer, int width, int height);
void triangle(CanvasTriangle t, int colour, bool filled, int width, int height, uint32_t *buffer);
int *loadPPM(string fileName, int &width, int &height);
void savePPM(string fileName, DrawingWindow *window);
void skipHashWS(ifstream &f);
void update(Camera &cam, vector<Updatable*> updatables);
void handleEvent(SDL_Event event, Camera &cam);
float depthBuffer[S_WIDTH * S_HEIGHT];
uint32_t imageBuffer[S_WIDTH * S_HEIGHT];
bool wireframe;
vector<float> Interpolate(float a, float b, int n);
vector<vec3> Interpolate(vec3 a, vec3 b, int n);
vector<vec4> Interpolate(vec4 a, vec4 b, int n);

bool toRaytrace = false;
bool softShadows = false;
DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);

inline float vectorLength(vec4 v) {
  return sqrt(v.x*v.x + v.y*v.y + v.z*v.z + v.w*v.w);
}

inline vec3 toThree(vec4 v) {
  return vec3(v.x, v.y, v.z);
}

inline vec4 cross(const vec4& a, const vec4& b) {
  return vec4(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x, 0);
}

// based on https://casual-effects.com/research/McGuire2011Clipping/McGuire-Clipping.pdf
int clipTriangle(vector<ModelTriangle>& tris, const vec4& normal) {
  vec4 temp;
  float tempBright;
  vector<int> toBeCulled;
  int n = tris.size();
  for (int i = 0; i < n; i++) {
    vector<float> distances;
    distances.push_back(dot(tris[i].vertices[0], normal));
    distances.push_back(dot(tris[i].vertices[1], normal));
    distances.push_back(dot(tris[i].vertices[2], normal));
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
      tempBright = tris[i].brightness[0];
      tris[i].vertices[0] = tris[i].vertices[1];
      tris[i].brightness[0] = tris[i].brightness[1];
      tris[i].vertices[1] = tris[i].vertices[2];
      tris[i].brightness[1] = tris[i].brightness[2];
      tris[i].vertices[2] = temp;
      tris[i].brightness[2] = tempBright;
      rotate(distances.begin(),distances.begin()+1,distances.end());
    }
    else if (distances[2] >= 0.0f && distances[1] < 0.0f) {
      nextInside = (distances[0] >= 0.0f);
      temp = tris[i].vertices[2];
      tempBright = tris[i].brightness[2];
      tris[i].vertices[2] = tris[i].vertices[1];
      tris[i].brightness[2] = tris[i].brightness[1];
      tris[i].vertices[1] = tris[i].vertices[0];
      tris[i].brightness[1] = tris[i].brightness[0];
      tris[i].vertices[0] = temp;
      tris[i].brightness[0] = tempBright;
      rotate(distances.begin(),distances.begin()+2,distances.end());
    }
    else {
      nextInside = (distances[1] >= 0.0f);
    }
    temp = mix(tris[i].vertices[0], tris[i].vertices[2], (distances[0] / (distances[0] - distances[2])));
    tempBright = mix(tris[i].brightness[0], tris[i].brightness[2], (distances[0] / (distances[0] - distances[2])));
    if (nextInside) {
      tris[i].vertices[2] = mix(tris[i].vertices[1], tris[i].vertices[2], (distances[1] / (distances[1] - distances[2])));
      tris[i].brightness[2] = mix(tris[i].brightness[1], tris[i].brightness[2], (distances[1] / (distances[1] - distances[2])));
      tris.push_back(ModelTriangle(tris[i].vertices[0], tris[i].vertices[2], temp, tris[i].brightness[0], tris[i].brightness[2], tempBright, tris[i].colour, tris[i].normal));
    }
    else {
      tris[i].vertices[1] = mix(tris[i].vertices[0], tris[i].vertices[1], (distances[0] / (distances[0] - distances[1])));
      tris[i].brightness[1] = mix(tris[i].brightness[0], tris[i].brightness[1], (distances[0] / (distances[0] - distances[1])));
      tris[i].vertices[2] = temp;
      tris[i].brightness[2] = tempBright;
    }
  }
  for (auto i : toBeCulled) {
    tris.erase(tris.begin() + i);
  }
  return tris.size();
}

int clipToView(vector<ModelTriangle>& tris) {
  const vec4 normals[2] = {vec4(0, 0, 1, 1), vec4(0, 0, -1, 1)}; // swap these two lines for full frustrum clipping
  // const vec4 normals[6] = {vec4(1, 0, 0, 1), vec4(-1, 0, 0, 1), vec4(0, 1, 0, 1), vec4(0, -1, 0, 1), vec4(0, 0, 1, 1), vec4(0, 0, -1, 1)};
  for (auto n : normals) {
    int num = clipTriangle(tris, n);
    if (num == 0) {
      return 0;
    }
  }
  return tris.size();
}

void drawTriangles(Camera &cam, std::vector<Model *> models)
{
  int width = (SSAA) ? S_WIDTH : WIDTH;
  int height = (SSAA) ? S_HEIGHT : HEIGHT;
  uint32_t *buffer = (SSAA) ? imageBuffer : window.pixelBuffer;
  for (unsigned int i = 0; i < models.size(); i++)
  {
    Model &model = *models[i];
    mat4 MVP = cam.projection * cam.worldToCamera() * model.transform;
    vec4 eye = vec4(cam.getPosition(), 0);
    for (auto tri : model.tris)
    {
      if (dot((model.transform * tri.vertices[0]) - eye, tri.normal) >= 0.0f) continue;
      tri.brightness[0] = glm::max(dot(normalize(eye - (model.transform * tri.vertices[0])), tri.normal), 0.0f);
      tri.brightness[1] = glm::max(dot(normalize(eye - (model.transform * tri.vertices[1])), tri.normal), 0.0f);
      tri.brightness[2] = glm::max(dot(normalize(eye - (model.transform * tri.vertices[2])), tri.normal), 0.0f);
      tri.vertices[0] = MVP * tri.vertices[0];
      tri.vertices[1] = MVP * tri.vertices[1];
      tri.vertices[2] = MVP * tri.vertices[2];
      vector<ModelTriangle> clippedTris;
      clippedTris.push_back(tri);
      clipToView(clippedTris);

      for (auto t : clippedTris) {
        t.vertices[0] /= t.vertices[0].w;
        t.vertices[1] /= t.vertices[1].w;
        t.vertices[2] /= t.vertices[2].w;
        CanvasPoint v1 = CanvasPoint(
          (t.vertices[0].x + 1.0f) * 0.5f * width,
          (1 - (t.vertices[0].y + 1.0f) * 0.5f) * height,
          ((cam.far - cam.near) / 2.0f) * t.vertices[0].z + ((cam.far + cam.near) / 2.0f),
          t.brightness[0]);
        CanvasPoint v2 = CanvasPoint(
          (t.vertices[1].x + 1.0f) * 0.5f * width,
          (1 - (t.vertices[1].y + 1.0f) * 0.5f) * height,
          ((cam.far - cam.near) / 2.0f) * t.vertices[1].z + ((cam.far + cam.near) / 2.0f),
          t.brightness[1]);
        CanvasPoint v3 = CanvasPoint(
          (t.vertices[2].x + 1.0f) * 0.5f * width,
          (1 - (t.vertices[2].y + 1.0f) * 0.5f) * height,
          ((cam.far - cam.near) / 2.0f) * t.vertices[2].z + ((cam.far + cam.near) / 2.0f),
          t.brightness[2]);
        triangle(CanvasTriangle(v1, v2, v3), tri.colour.toPackedInt(), wireframe, width, height, buffer);
      }
    }
  }
}

void downsample() {
  uint32_t pixel, subpixel, red, green, blue;
  int num_samples = SSAA_SCALE * SSAA_SCALE;
  for (int y = 0; y < HEIGHT; y++) {
    for (int x = 0; x < WIDTH; x++) {
      red = 0;
      green = 0;
      blue = 0;
      for (int sy = 0; sy < SSAA_SCALE; sy++) {
        for (int sx = 0; sx < SSAA_SCALE; sx++) {
          subpixel = imageBuffer[(x * SSAA_SCALE + sx) + (y * SSAA_SCALE + sy) * S_WIDTH];
          red += (subpixel & 0x00ff0000);
          green += (subpixel & 0x0000ff00);
          blue += (subpixel & 0x000000ff);
        }
      }
      pixel = 0xff000000 | ((red / num_samples) & 0x00ff0000) | ((green / num_samples) & 0x0000ff00) | ((blue / num_samples) & 0x000000ff);
      window.pixelBuffer[x + y * WIDTH] = pixel;
    }
  }
}

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

vector<ModelTriangle> getLights(Model model) {
  vector<ModelTriangle> lights = vector<ModelTriangle>();
  for(ModelTriangle triangle : model.tris) {
    if(triangle.name == "light") lights.push_back(triangle);
  }

  return lights;
}

vector<vec4> getLightPoints(vector<ModelTriangle> lights) {
  vector<vec4> vertices = vector<vec4>();
  vec4 average = vec4(0, 0, 0, 0);
  for(ModelTriangle light : lights) {
    vertices.push_back(light.vertices[0]);
    vertices.push_back(light.vertices[1]);
    vertices.push_back(light.vertices[2]);

    average += light.vertices[0];
    average += light.vertices[1];
    average += light.vertices[2];
  }

  average.x /= vertices.size();
  average.y /= vertices.size();
  average.z /= vertices.size();
  average.w = 1;

  vector<vec4> result = vector<vec4>();
  result.push_back(average);

  return result;
}

void raytrace(Camera camera, std::vector<Model*> models, int softness) {
  //lighting options. Make ambience global and the others properties of the lights
  float intensity = 10.0f / softness;
  float shadow = 0.1f / softness;

  //idk what this does
  for (unsigned int i = 1; i < models.size(); i++) {
    for (auto tri : (*models[i]).tris) {
      (*models[0])
          .tris.push_back(
              ModelTriangle((*models[i]).transform * tri.vertices[0],
                            (*models[i]).transform * tri.vertices[1],
                            (*models[i]).transform * tri.vertices[2],
                            tri.colour, tri.normal));
    }
  }
  Model &model = (*models[0]);

  //gets all the triangles that are lights in the scene
  vector<ModelTriangle> lights = getLights(model);

  vector<vec4> lightPoints = getLightPoints(lights);

  //loop through each pixel in image plane
  for(int j = 0; j < HEIGHT; j++) {
    for(int i = 0; i < WIDTH; i++) {
      float angle = tanf(0.5f * glm::radians(camera.fov)); // just fov*0.5 converted to radians
      //convert image plane cordinates into world space
      vec2 NDC = vec2((i + 0.5) * (1 / (float) WIDTH), (j + 0.5) * (1 / (float) HEIGHT));
      float x = (2 * (NDC.x) - 1) * angle * ASPECT_RATIO;
      float y = (1 - 2 * (NDC.y)) * angle;

      // the main camera ray
      vec4 rayDirection = camera.transform * vec4(x, y, -1, 0);

      float minDistance = std::numeric_limits<float>::infinity();
      
      RayTriangleIntersection intersection = RayTriangleIntersection();
      bool foundIntersection = false;

      //calculate closest intersection by looping through each of the triangles
      for(ModelTriangle triangle : model.tris) {
        vec4 e0 = triangle.vertices[1] - triangle.vertices[0];
        vec4 e1 = triangle.vertices[2] - triangle.vertices[0];
        vec4 SPVector = camera.transform[3] - triangle.vertices[0];
        mat4 DEMatrix(-rayDirection, e0, e1, vec4(1, 1, 1, 1));
        vec4 possibleSolution = glm::inverse(DEMatrix) * SPVector;

        // check if ray intersects triangle and not just triangle plane
        if (possibleSolution.y >= 0 && possibleSolution.y <= 1 &&
            possibleSolution.z >= 0 && possibleSolution.z <= 1 &&
            possibleSolution.y + possibleSolution.z <= 1) {
          if (possibleSolution.x < minDistance) {
            foundIntersection = true;
            intersection = RayTriangleIntersection(camera.transform[3] + (possibleSolution.x * rayDirection) , possibleSolution.x, triangle);
            intersection.intersectionPoint.w = 1;
            minDistance = possibleSolution.x;
          }
        }
      }

      if(foundIntersection) {
        //TODO: make special float that automatically binds between 0 and 1
        //these values lighten the pixel, so they go from 0 (dark) to 1 (fully in light)
        float brightnessCount = 0.0f;
        float angleCount = 0.0f;
        float specularCount = 0.0f;

        float shadowCount = shadow;
        
        bool inShadow = false;

        //fires a shadow ray to each light point
        for(vec4 light : lightPoints) {
          vec4 shadowRayDirection = light - intersection.intersectionPoint;
          vec4 shadowRayNormalised = glm::normalize(shadowRayDirection);

          //cross and dot only work on vec4s
          vec4 intersectionNormal = glm::normalize(cross(intersection.intersectedTriangle.vertices[1] - intersection.intersectedTriangle.vertices[0],
                                    intersection.intersectedTriangle.vertices[2] - intersection.intersectedTriangle.vertices[0]));

          //calculate the angleOfIncidence between 0 and 1
          float angleOfIncidence = glm::dot(shadowRayNormalised, intersectionNormal);
          angleOfIncidence = angleOfIncidence < 0 ? AMBIENCE : angleOfIncidence;
          angleCount += angleOfIncidence;

          //adjust brightness for proximity lighting
          float brightness = intensity/pow(vectorLength(shadowRayDirection),2);
          brightnessCount += brightness;

          //128 will later have to be paramaterised to reflect each material
          vec4 reflection = glm::normalize((-shadowRayDirection) - 2.0f * (glm::dot((-shadowRayDirection), intersectionNormal) * intersectionNormal));
          float specular = pow(glm::dot(glm::normalize((-rayDirection)), reflection), 128);

          specularCount += specular;

          float shadowBias = 1e-4;

          //check if the ray is in shadow. 
          for(ModelTriangle triangle : model.tris) {
            vec4 e0 = triangle.vertices[1] - triangle.vertices[0];
            vec4 e1 = triangle.vertices[2] - triangle.vertices[0];
            vec4 SPVector = (intersection.intersectionPoint) - triangle.vertices[0];
            mat4 DEMatrix(-shadowRayDirection, e0, e1, vec4(1,1,1,1));
            vec4 possibleSolution = glm::inverse(DEMatrix) * SPVector;

            //check if ray intersects triangle and not just triangle plane
            if(possibleSolution.y >= 0 && possibleSolution.y <= 1 && possibleSolution.z >= 0 && possibleSolution.z <= 1 && possibleSolution.y + possibleSolution.z <= 1) {
              //I genuinely have no idea why doing the shadowBias this way works. without it the shadows are just everywhere???
              if(possibleSolution.x < 1 - shadowBias && possibleSolution.x > shadowBias) {
                inShadow = true;
                break;
              }
            }
          }

          shadowCount -= inShadow ? shadow/10 : 0;
        }

        //adjust the totalled lighting values
        shadowCount = clamp<float>(shadowCount, 0, 1);
        brightnessCount = clamp<float>(brightnessCount, 0, 1);
        angleCount = clamp<float>(angleCount, 0, 1);
        specularCount = clamp<float>(specularCount, 0, 1);

        //set the final pixels
        if(intersection.intersectedTriangle.name == "light") window.setPixelColour(i, j, intersection.intersectedTriangle.colour.toPackedInt());
        else window.setPixelColour(i, j, darkenColour(intersection.intersectedTriangle.colour, clamp<float>(angleCount * (inShadow ? shadowCount : 1.0f) * brightnessCount, AMBIENCE, 1), inShadow? 0 : specularCount));
      } else {
        window.setPixelColour(i, j, 0);
      }
    }
  }
}

//vector<mat4> cameraTransforms = vector<mat4>();

int main(int argc, char *argv[])
{
  SDL_Event event;
  //SDL_SetRelativeMouseMode(SDL_TRUE);

  vector<Model*> renderQueue = vector<Model*>();
  vector<Updatable*> updateQueue = vector<Updatable*>();

  Model cornell = Model("cornell-box");
  renderQueue.push_back(&cornell);
  // std::cout << "cornell address = " << &cornell << std::endl;
  Rigidbody cornellRB = Rigidbody(&cornell);
  cornellRB.hasGravity = false;
  updateQueue.push_back(&cornellRB);

  Model hs_logo = Model("HackspaceLogo/logo");
  hs_logo.scale(vec3(0.005f, 0.005f, 0.005f));
  hs_logo.move(vec3(-1.1f, 1.21f, -1.8f));
  renderQueue.push_back(&hs_logo);
  // Model cornell2 = Model("cornell-box");
  // // cornell2.move(glm::vec3(0,1,0));
  // cornell2.rotate(glm::vec3(0,1,0));
  // renderQueue.push_back(&cornell2);
  // Rigidbody cornellRB2 = Rigidbody(&cornell2);
  // cornellRB2.hasGravity = false;
  // updateQueue.push_back(&cornellRB2);

  Model sphere = Model("blob");
  sphere.setPosition(vec3(0,5.5f,-3));
  renderQueue.push_back(&sphere);
  Rigidbody sphereRB = Rigidbody(&sphere);
  updateQueue.push_back(&sphereRB);

  // Model tri1 = Model("triangle");
  // renderQueue.push_back(&tri1);
  // Rigidbody tri1RB = Rigidbody(&tri1);
  // tri1RB.hasGravity = false;
  // updateQueue.push_back(&tri1RB);

  // Model tri2 = Model("triangle2");
  // tri2.move(glm::vec3(0,0.5f,0));
  // tri2.rotate(glm::vec3(1, 1, 0));
  // renderQueue.push_back(&tri2);
  // Rigidbody tri2RB = Rigidbody(&tri2);
  // tri2RB.hasGravity = false;
  // updateQueue.push_back(&tri2RB);

  //std::cout << "address stored as " << cornellRB.model << std::endl;

  Camera cam;
  cam.setProjection(90.0f, WIDTH / (float)HEIGHT, 0.1f, 100.0f);
  cam.lookAt(vec3(0.0f, 2.5f, 3.0f), vec3(0.0f, 2.5f, 0.0f));
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
    //std::cout << "sphere transform = " << sphere.transform << std::endl;
    //update(cam, vector<Updatable*>{&cornell, &cornellRB, &sphereRB});
    update(cam, updateQueue);
    //std::vector<Model*> models{&cornell, &sphere};
    //std::cout << "about to render" << std::endl;
    if(toRaytrace) {
      raytrace(cam, renderQueue, softShadows ? 2 : 1);
    } else {
      draw();
      drawTriangles(cam, renderQueue);
      if (SSAA) downsample();
    }
    // Need to render the frame at the end, or nothing actually gets shown on
    // the screen !
    window.renderFrame();
  }
}

void draw()
{
  if (!SSAA) window.clearPixels();
  for (int i = 0; i < S_WIDTH * S_HEIGHT; i++)
  {
    depthBuffer[i] = std::numeric_limits<float>::infinity();
    imageBuffer[i] = 0;
  }
}

int moveStage = 0;

void update(Camera &cam, vector<Updatable*> updatables)
{
  // Function for performing animation (shifting artifacts or moving the camera)
  cam.update();
  for (unsigned int i = 0; i < updatables.size(); i++)
  {
    updatables[i]->update();
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
    cout << cam.getPosition() << '\n';
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

// Line clipping based on https://www.geeksforgeeks.org/line-clipping-set-1-cohen-sutherland-algorithm/
int clipCode(CanvasPoint& p, int width, int height) {
  int code = 0;
  if (p.x < 0.0f) code |= LEFT;
  if (p.y < 0.0f) code |= TOP;
  if (p.x > width - 1.0f) code |= RIGHT;
  if (p.y > height - 1.0f) code |= BOTTOM;
  return code;
}

bool clipLine(CanvasPoint& p, CanvasPoint& q, int width, int height) {
  int p_code = clipCode(p, width, height);
  int q_code = clipCode(q, width, height);
  while (true) {
    if (p_code == 0 && q_code == 0) return true;
    if (p_code & q_code) return false;
    int current_code = p_code ? p_code : q_code;
    float x = 0.0f, y = 0.0f;
    if (current_code & TOP) {
      x = mix(p.x, q.x, (-p.y) / (q.y - p.y));
      y = 0.0f;
    }
    else if (current_code & RIGHT) {
      y = mix(p.y, q.y, ((width - 1.0f) - p.x) / (q.x - p.x));
      x = width - 1.0f;
    }
    else if (current_code & BOTTOM) {
      x = mix(p.x, q.x, ((height - 1.0f) - p.y) / (q.y - p.y));
      y = height - 1.0f;
    }
    else if (current_code & LEFT) {
      y = mix(p.y, q.y, (-p.x) / (q.x - p.x));
      x = 0.0f;
    }
    if (current_code == p_code) {
      p.x = x;
      p.y = y;
      p_code = clipCode(p, width, height);
    }
    else {
      q.x = x; 
      q.y = y; 
      q_code = clipCode(q, width, height);
    }
  }
}

void line(CanvasPoint p, CanvasPoint q, int colour, uint32_t *buffer, int width, int height)
{
  if (!clipLine(p, q, width, height)) return;
  float x_diff = p.x - q.x;
  float y_diff = p.y - q.y;
  float max_diff = std::max(abs(x_diff), abs(y_diff));
  vector<float> interpolate_x = Interpolate(p.x, q.x, max_diff + 1);
  vector<float> interpolate_y = Interpolate(p.y, q.y, max_diff + 1);
  for (int i = 0; i < max_diff; i++)
  {
    buffer[(int)(round(interpolate_y[i]) * width + round(interpolate_x[i]))] = colour;
  }
}

inline float edgeFunction(const CanvasPoint& v0, const CanvasPoint& v1, const CanvasPoint& p) {
  return (p.x - v0.x) * (v1.y - v0.y) - (p.y - v0.y) * (v1.x - v0.x);
}

inline int scaleColour(int colour, float scale) {
  unsigned char red = (colour & 0x00ff0000) >> 16;
  red *= scale;
  unsigned char green = (colour & 0x0000ff00) >> 8;
  green *= scale;
  unsigned char blue = (colour & 0x000000ff);
  blue *= scale;
  return (colour & 0xff000000) | (red << 16) | (green << 8) | blue;
}

void triangle(CanvasTriangle t, int colour, bool filled, int width, int height, uint32_t *buffer)
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
    x_min = clamp<int>(x_min, 0, width - 1);
    x_max = clamp<int>(x_max, 0, width - 1);
    y_min = clamp<int>(y_min, 0, height - 1);
    y_max = clamp<int>(y_max, 0, height - 1);
    float area_inv = 1.0f / edgeFunction(t.vertices[0], t.vertices[1], t.vertices[2]);
    float w0_step_x = (t.vertices[2].y - t.vertices[1].y) * area_inv;
    float w1_step_x = (t.vertices[0].y - t.vertices[2].y) * area_inv;
    float w2_step_x = (t.vertices[1].y - t.vertices[0].y) * area_inv;
    float w0_step_y = (t.vertices[1].x - t.vertices[2].x) * area_inv;
    float w1_step_y = (t.vertices[2].x - t.vertices[0].x) * area_inv;
    float w2_step_y = (t.vertices[0].x - t.vertices[1].x) * area_inv;
    CanvasPoint p = CanvasPoint(x_min + 0.5f, y_min + 0.5f);
    float w0_line = edgeFunction(t.vertices[1], t.vertices[2], p) * area_inv;
    float w1_line = edgeFunction(t.vertices[2], t.vertices[0], p) * area_inv;
    float w2_line = edgeFunction(t.vertices[0], t.vertices[1], p) * area_inv;
    float w0, w1, w2;
    for (int y = y_min; y <= y_max; y++) {
      w0 = w0_line;
      w1 = w1_line;
      w2 = w2_line;
      for (int x = x_min; x <= x_max; x++) {
        if (w0 >= 0 && w1 >= 0 && w2 >= 0) {
          float depth = w0 * t.vertices[0].depth + w1 * t.vertices[1].depth + w2 * t.vertices[2].depth;
          float brightness = w0 * t.vertices[0].brightness + w1 * t.vertices[1].brightness + w2 * t.vertices[2].brightness;
          if (depth < depthBuffer[y * width + x]) {
            depthBuffer[y * width + x] = depth;
            buffer[y * width + x] = scaleColour(colour, brightness);
          }
        }
        w0 += w0_step_x;
        w1 += w1_step_x;
        w2 += w2_step_x;
      }
      w0_line += w0_step_y;
      w1_line += w1_step_y;
      w2_line += w2_step_y;
    }
  }
  else
  {
    line(t.vertices[0], t.vertices[1], colour, buffer, width, height);
    line(t.vertices[1], t.vertices[2], colour, buffer, width, height);
    line(t.vertices[2], t.vertices[0], colour, buffer, width, height);
  }
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