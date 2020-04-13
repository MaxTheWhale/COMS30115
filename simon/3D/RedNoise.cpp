#include <ModelTriangle.h>
#include <RayTriangleIntersection.h>
#include <DrawingWindow.h>
#include <Utils.h>
#include <fstream>
#include <glm/glm.hpp>
#include <unordered_map>
#include <vector>
#include <list>
#include <chrono>
#include <sys/time.h>
#include "Camera.hpp"
#include "Model.hpp"
#include "Times.hpp"
#include "VectorOutput.hpp"
#include "Rigidbody.hpp"
#include "Light.h"
using namespace std;
using namespace glm;

#define WIDTH 640
#define HEIGHT 480
#define IMG_SIZE (WIDTH*HEIGHT)
#define SSAA true
#define SSAA_SCALE 3
#define SSAA_SAMPLES (SSAA_SCALE*SSAA_SCALE)
#define MOUSE_SENSITIVITY 0.0015f
#define AMBIENCE 0.1f
#define ASPECT_RATIO WIDTH/(float)HEIGHT

enum CLIP_CODE {TOP = 1, RIGHT = 2, BOTTOM = 4, LEFT = 8};
enum COLOUR_MASK {ALPHA = 0xff000000, RED = 0x00ff0000, GREEN = 0x0000ff00, BLUE = 0x000000ff};

class Vertex {
  public:
    vec4 pos;
    float brightness;
    float u, v;
    Vertex operator+=(const Vertex& rhs)
    {
      pos += rhs.pos;
      brightness += rhs.brightness;
      u += rhs.u;
      v += rhs.v;
      return *this;
    }

    friend Vertex operator+(Vertex lhs, const Vertex& rhs)
    {
      lhs += rhs;
      return lhs;
    }

    Vertex operator-=(const Vertex& rhs)
    {
      pos -= rhs.pos;
      brightness -= rhs.brightness;
      u -= rhs.u;
      v -= rhs.v;
      return *this;
    }

    friend Vertex operator-(Vertex lhs, const Vertex& rhs)
    {
      lhs -= rhs;
      return lhs;
    }

    Vertex operator*=(float rhs)
    {
      pos *= rhs;
      brightness *= rhs;
      u *= rhs;
      v *= rhs;
      return *this;
    }

    friend Vertex operator*(Vertex lhs, float rhs)
    {
      lhs *= rhs;
      return lhs;
    }

    Vertex operator/=(float rhs)
    {
      pos /= rhs;
      brightness /= rhs;
      u /= rhs;
      v /= rhs;
      return *this;
    }

    friend Vertex operator/(Vertex lhs, float rhs)
    {
      lhs /= rhs;
      return lhs;
    }
};

class Triangle {
  public:
    Vertex vertices[3];
    int colour;
    vec4 normal;
    Triangle(ModelTriangle &tri) {
      for (int i = 0; i < 3; i++) {
        vertices[i].pos.x = tri.vertices[i].x;
        vertices[i].pos.y = tri.vertices[i].y;
        vertices[i].pos.z = tri.vertices[i].z;
        vertices[i].pos.w = tri.vertices[i].w;
        vertices[i].u = tri.uvs[i].x;
        vertices[i].v = tri.uvs[i].y;
        vertices[i].brightness = tri.brightness[i];
      }
      normal = tri.normal;
      colour = tri.material.diffuse.toPackedInt();
    }
    Triangle(const Vertex &v0, const Vertex &v1, const Vertex &v2, const int &tColour, const vec4 &tNormal) {
      vertices[0] = v0;
      vertices[1] = v1;
      vertices[2] = v2;
      colour = tColour;
      normal = tNormal;
    }
};

void draw();
void line(vec4 p, vec4 q, int colour, uint32_t *buffer, vec2 &offset);
void triangle(Triangle &t, Texture &tex, int colour, bool filled, uint32_t *buffer, float *depthBuff, vec2 offset);
int *loadPPM(string fileName, int &width, int &height);
void savePPM(string fileName, DrawingWindow *window);
void skipHashWS(ifstream &f);
void update(Camera &cam, vector<Updatable*> updatables);
void handleEvent(SDL_Event event, Camera &cam);
#if SSAA
float depthBuffer[IMG_SIZE * SSAA_SAMPLES];
uint32_t imageBuffer[IMG_SIZE * SSAA_SAMPLES];
#else
float depthBuffer[IMG_SIZE];
uint32_t imageBuffer[IMG_SIZE];
#endif
bool wireframe;
bool bilinear;
bool perspective;
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

inline Vertex mixVertex(const Vertex &p, const Vertex &q, float a) {
  return p * (1.0f - a) + q * a;
}

// based on https://casual-effects.com/research/McGuire2011Clipping/McGuire-Clipping.pdf
int clipTriangle(list<Triangle>& tris, const vec4& normal) {
  Vertex temp;
  list<Triangle>::iterator tri = tris.begin();
  int n = tris.size();
  for (int i = 0; (i < n) && tri != tris.end(); i++) {
    vector<float> distances;
    distances.push_back(dot((*tri).vertices[0].pos, normal));
    distances.push_back(dot((*tri).vertices[1].pos, normal));
    distances.push_back(dot((*tri).vertices[2].pos, normal));
    if (distances[0] >= 0.0f && distances[1] >= 0.0f && distances[2] >= 0.0f) {
      continue;
    }
    if (distances[0] < 0.0f && distances[1] < 0.0f && distances[2] < 0.0f) {
      tris.erase(tri);
    }
    bool nextInside;
    if (distances[1] >= 0.0f && distances[0] < 0.0f) {
      nextInside = (distances[2] >= 0.0f);
      temp = (*tri).vertices[0];
      (*tri).vertices[0] = (*tri).vertices[1];
      (*tri).vertices[1] = (*tri).vertices[2];
      (*tri).vertices[2] = temp;
      rotate(distances.begin(),distances.begin()+1,distances.end());
    }
    else if (distances[2] >= 0.0f && distances[1] < 0.0f) {
      nextInside = (distances[0] >= 0.0f);
      temp = (*tri).vertices[2];
      (*tri).vertices[2] = (*tri).vertices[1];
      (*tri).vertices[1] = (*tri).vertices[0];
      (*tri).vertices[0] = temp;
      rotate(distances.begin(),distances.begin()+2,distances.end());
    }
    else {
      nextInside = (distances[1] >= 0.0f);
    }
    temp = mixVertex((*tri).vertices[0], (*tri).vertices[2], (distances[0] / (distances[0] - distances[2])));
    if (nextInside) {
      (*tri).vertices[2] = mixVertex((*tri).vertices[1], (*tri).vertices[2], (distances[1] / (distances[1] - distances[2])));
      Triangle newTri = Triangle((*tri).vertices[0], (*tri).vertices[2], temp, (*tri).colour, (*tri).normal);
      tris.push_back(newTri);
    }
    else {
      (*tri).vertices[1] = mixVertex((*tri).vertices[0], (*tri).vertices[1], (distances[0] / (distances[0] - distances[1])));
      (*tri).vertices[2] = temp;
    }
    tri++;
  }
  return tris.size();
}

int clipToView(list<Triangle>& tris) {
  const vec4 normals[6] = {vec4(1, 0, 0, 1), vec4(-1, 0, 0, 1), vec4(0, 1, 0, 1), vec4(0, -1, 0, 1), vec4(0, 0, 1, 1), vec4(0, 0, -1, 1)};
  for (auto n : normals) {
    int num = clipTriangle(tris, n);
    if (num == 0) {
      return 0;
    }
  }
  return tris.size();
}

vector<vec2> generateRotatedGrid(int gridSize) {
  vector<vec2> result;
  const int numSamples = gridSize * gridSize;
  const float step = 1.0f / numSamples;
  float x = (step * 0.5f) + step * (numSamples - gridSize);
  for (float y = step * 0.5f; y < 1.0f; y += step) {
    result.push_back(vec2(x, y));
    x -= step * gridSize;
    if (x < 0)
      x += step * (numSamples + 1);
  }
  return result;
}

void drawTriangles(Camera &cam, std::vector<Model *> models)
{
  uint32_t *buffer = (SSAA) ? imageBuffer : window.pixelBuffer;
  vector<vec2> offsets = generateRotatedGrid(SSAA_SCALE);
  for (unsigned int i = 0; i < models.size(); i++)
  {
    Model &model = *models[i];
    mat4 MVP = cam.projection * cam.worldToCamera() * model.transform;
    vec4 eye = vec4(cam.getPosition(), 0);
    for (auto& modelTri : model.tris)
    {
      Triangle tri = Triangle(modelTri);
      if (dot((model.transform * tri.vertices[0].pos) - eye, model.transform * tri.normal) >= 0.0f) continue;
      tri.vertices[0].brightness = glm::max(dot(normalize(eye - (model.transform * tri.vertices[0].pos)), normalize(model.transform * tri.normal)), 0.0f);
      tri.vertices[1].brightness = glm::max(dot(normalize(eye - (model.transform * tri.vertices[1].pos)), normalize(model.transform * tri.normal)), 0.0f);
      tri.vertices[2].brightness = glm::max(dot(normalize(eye - (model.transform * tri.vertices[2].pos)), normalize(model.transform * tri.normal)), 0.0f);
      tri.vertices[0].pos = MVP * tri.vertices[0].pos;
      tri.vertices[1].pos = MVP * tri.vertices[1].pos;
      tri.vertices[2].pos = MVP * tri.vertices[2].pos;
      list<Triangle> clippedTris;
      clippedTris.push_back(tri);
      clipToView(clippedTris);
      for (auto t : clippedTris) {
        for (int v = 0; v < 3; v++) {
          t.vertices[v].pos.w = 1.0f / t.vertices[v].pos.w;
          t.vertices[v].pos.x *= t.vertices[v].pos.w;
          t.vertices[v].pos.y *= t.vertices[v].pos.w;
          t.vertices[v].pos.z *= t.vertices[v].pos.w;
          t.vertices[v].pos.x = (t.vertices[v].pos.x + 1.0f) * 0.5f * WIDTH;
          t.vertices[v].pos.y = (1 - (t.vertices[v].pos.y + 1.0f) * 0.5f) * HEIGHT;
          t.vertices[v].pos.z = ((cam.far - cam.near) / 2.0f) * t.vertices[v].pos.z + ((cam.far + cam.near) / 2.0f);
        }
        #if SSAA
        for (int s = 0; s < SSAA_SAMPLES; s++) {
          triangle(t, model.texture, tri.colour, wireframe, buffer + (IMG_SIZE * s), depthBuffer + (IMG_SIZE * s), offsets[s]);
        }
        #else
        triangle(t, model.texture, tri.colour, wireframe, buffer, depthBuffer, vec2(0.5f, 0.5f));
        #endif
      }
    }
  }
}

void downsample(uint32_t *source, uint32_t *dest, int width, int height, int num_samples) {
  uint32_t pixel, sample, red, green, blue;
  for (int y = 0; y < height; y++) {
    for (int x = 0; x < width; x++) {
      red = 0;
      green = 0;
      blue = 0;
      for (int s = 0; s < num_samples; s++) {
        sample = source[(s * width * height) + x + y * width];
        red += (sample & RED);
        green += (sample & GREEN);
        blue += (sample & BLUE);
      }
      pixel = ALPHA | ((red / num_samples) & RED) | ((green / num_samples) & GREEN) | ((blue / num_samples) & BLUE);
      dest[x + y * width] = pixel;
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

Colour addColours(Colour base, Colour add) {
  base.red += add.red;
  base.green += add.green;
  base.blue += add.blue;
  return base;
}

RayTriangleIntersection findClosestIntersection(Camera camera, Model model, vec4 rayDirection) {
  RayTriangleIntersection intersection = RayTriangleIntersection();
  float minDistance = std::numeric_limits<float>::infinity();

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
        intersection = RayTriangleIntersection(camera.transform[3] + (possibleSolution.x * rayDirection) , possibleSolution.x, triangle);
        intersection.intersectionPoint.w = 1;
        minDistance = possibleSolution.x;
      }
    }
  }

  return intersection;
}

bool inShadow(Model model, vec4 shadowRayDirection, RayTriangleIntersection intersection) {
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
        return true;
      }
    }
  }

  return false;
}

void raytrace(Camera camera, std::vector<Model*> models) {
  //idk what this does
  Model model = Model(*models[0]);
  for (unsigned int i = 1; i < models.size(); i++) {
    for (auto tri : (*models[i]).tris) {
      //if (dot(((*models[i]).transform * tri.vertices[0]) - vec4(camera.getPosition(), 0), tri.normal) >= 0.0f) continue;
      model
          .tris.push_back(
              ModelTriangle((*models[i]).transform * tri.vertices[0],
                            (*models[i]).transform * tri.vertices[1],
                            (*models[i]).transform * tri.vertices[2],
                            tri.material, tri.normal));
    }
  }

  //setup the lights
  vector<ModelTriangle> lights = vector<ModelTriangle>();
  for(ModelTriangle triangle : model.tris) {
    if(triangle.name == "light") lights.push_back(triangle);
  }

  Light mainLight = Light("light", lights);
  mainLight.calculateCentre();
  mainLight.centre.y += 0.1f;

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

      RayTriangleIntersection intersection = findClosestIntersection(camera, model, rayDirection);

      Colour colour = Colour(0, 0, 0);

      if(intersection.intersectedTriangle.name == mainLight.name) {
        window.setPixelColour(i, j, intersection.intersectedTriangle.material.diffuse.toPackedInt());
      } else if(intersection.intersectedTriangle.name != "") {
        colour = intersection.intersectedTriangle.material.diffuse;

        vec4 shadowRayDirection = mainLight.centre - intersection.intersectionPoint;
        bool isInShadow = inShadow(model, shadowRayDirection, intersection);

        if(isInShadow) colour = colour * mainLight.shadow;
        else {
          //calculate the angleOfIncidence between 0 and 1
          float angleOfIncidence = glm::dot(glm::normalize(shadowRayDirection), intersection.intersectedTriangle.normal);
          colour = colour * clamp<float>(angleOfIncidence, mainLight.shadow, 1);

          //adjust brightness for proximity lighting
          float brightness = mainLight.intensity/pow(vectorLength(shadowRayDirection),2);
          colour = colour * clamp<float>(brightness, mainLight.shadow, 1);
        }

        window.setPixelColour(i, j, colour.toPackedInt());
        
      } else {
        window.setPixelColour(i, j, colour.toPackedInt());
      }
    }
  }
  cout << "Finished one frame!" << endl;
}

// void raytraceOld(Camera camera, std::vector<Model*> models, int softness) {
//   //lighting options. Make ambience global and the others properties of the lights
//   float intensity = 10.0f / softness;
//   float shadow = 0.1f / softness;

//   //idk what this does
//   Model model = Model(*models[0]);
//   for (unsigned int i = 1; i < models.size(); i++) {
//     for (auto tri : (*models[i]).tris) {
//       //if (dot(((*models[i]).transform * tri.vertices[0]) - vec4(camera.getPosition(), 0), tri.normal) >= 0.0f) continue;
//       model
//           .tris.push_back(
//               ModelTriangle((*models[i]).transform * tri.vertices[0],
//                             (*models[i]).transform * tri.vertices[1],
//                             (*models[i]).transform * tri.vertices[2],
//                             tri.material, tri.normal));
//     }
//   }

//   //gets all the triangles that are lights in the scene
//   vector<ModelTriangle> lights = getLights(model);

//   vector<vec4> lightPoints = getLightPoints(lights);

//   //loop through each pixel in image plane
//   for(int j = 0; j < HEIGHT; j++) {
//     for(int i = 0; i < WIDTH; i++) {
//       float angle = tanf(0.5f * glm::radians(camera.fov)); // just fov*0.5 converted to radians
//       //convert image plane cordinates into world space
//       vec2 NDC = vec2((i + 0.5) * (1 / (float) WIDTH), (j + 0.5) * (1 / (float) HEIGHT));
//       float x = (2 * (NDC.x) - 1) * angle * ASPECT_RATIO;
//       float y = (1 - 2 * (NDC.y)) * angle;

//       // the main camera ray
//       vec4 rayDirection = camera.transform * vec4(x, y, -1, 0);

//       float minDistance = std::numeric_limits<float>::infinity();
      
//       RayTriangleIntersection intersection = RayTriangleIntersection();
//       bool foundIntersection = false;
//       Colour transparentColour = Colour(0, 0, 0);

//       //calculate closest intersection by looping through each of the triangles
//       for(ModelTriangle triangle : model.tris) {
//         vec4 e0 = triangle.vertices[1] - triangle.vertices[0];
//         vec4 e1 = triangle.vertices[2] - triangle.vertices[0];
//         vec4 SPVector = camera.transform[3] - triangle.vertices[0];
//         mat4 DEMatrix(-rayDirection, e0, e1, vec4(1, 1, 1, 1));
//         vec4 possibleSolution = glm::inverse(DEMatrix) * SPVector;

//         // check if ray intersects triangle and not just triangle plane
//         if (possibleSolution.y >= 0 && possibleSolution.y <= 1 &&
//           possibleSolution.z >= 0 && possibleSolution.z <= 1 &&
//           possibleSolution.y + possibleSolution.z <= 1) {
//           if(triangle.name == "short_box") {
//             transparentColour.red += 200;
//           } else if (possibleSolution.x < minDistance) {
//             foundIntersection = true;
//             intersection = RayTriangleIntersection(camera.transform[3] + (possibleSolution.x * rayDirection) , possibleSolution.x, triangle);
//             intersection.intersectionPoint.w = 1;
//             minDistance = possibleSolution.x;
//           }
//         }
//       }

//       if(foundIntersection) {
//         bool mirror = false;

//         if(intersection.intersectedTriangle.name == "tall_box") {
//           // the main camera ray
//           vec4 mirrorRayDirection = glm::normalize((intersection.intersectedTriangle.normal) - 2.0f * (glm::dot((intersection.intersectedTriangle.normal), rayDirection) * rayDirection));
//           mirrorRayDirection.w = 0;

//           float minDistance = std::numeric_limits<float>::infinity();
          
//           RayTriangleIntersection mirrorIntersection = RayTriangleIntersection();
//           bool foundIntersection = false;

//           //calculate closest intersection by looping through each of the triangles
//           for(ModelTriangle triangle : model.tris) {
//             vec4 e0 = triangle.vertices[1] - triangle.vertices[0];
//             vec4 e1 = triangle.vertices[2] - triangle.vertices[0];
//             vec4 SPVector = camera.transform[3] - triangle.vertices[0];
//             mat4 DEMatrix(-mirrorRayDirection, e0, e1, vec4(1, 1, 1, 1));
//             vec4 possibleSolution = glm::inverse(DEMatrix) * SPVector;

//             // check if ray intersects triangle and not just triangle plane
//             if (possibleSolution.y >= 0 && possibleSolution.y <= 1 &&
//                 possibleSolution.z >= 0 && possibleSolution.z <= 1 &&
//                 possibleSolution.y + possibleSolution.z <= 1) {
//               if (possibleSolution.x < minDistance) {
//                 foundIntersection = true;
//                 mirrorIntersection = RayTriangleIntersection(camera.transform[3] + (possibleSolution.x * rayDirection) , possibleSolution.x, triangle);
//                 mirrorIntersection.intersectionPoint.w = 1;
//                 minDistance = possibleSolution.x;
//               }
//             }
//           }

//           if(foundIntersection) {
//             mirror = true;
//             intersection = mirrorIntersection;
//           } else {
//             break;
//           }
//         }


//         //TODO: make special float that automatically binds between 0 and 1
//         //these values lighten the pixel, so they go from 0 (dark) to 1 (fully in light)
//         float brightnessCount = 0.0f;
//         float angleCount = 0.0f;
//         float specularCount = 0.0f;

//         float shadowCount = shadow;
        
//         bool inShadow = false;

//         //fires a shadow ray to each light point
//         for(vec4 light : lightPoints) {
//           vec4 shadowRayDirection = light - intersection.intersectionPoint;
//           vec4 shadowRayNormalised = glm::normalize(shadowRayDirection);

//           //calculate the angleOfIncidence between 0 and 1
//           float angleOfIncidence = glm::dot(shadowRayNormalised, intersection.intersectedTriangle.normal);
//           angleOfIncidence = angleOfIncidence < 0 ? AMBIENCE : angleOfIncidence;
//           angleCount += angleOfIncidence;

//           //adjust brightness for proximity lighting
//           float brightness = intensity/pow(vectorLength(shadowRayDirection),2);
//           brightnessCount += brightness;

//           //128 will later have to be paramaterised to reflect each material
//           vec4 reflection = glm::normalize((-shadowRayDirection) - 2.0f * (glm::dot((-shadowRayDirection), intersection.intersectedTriangle.normal) * intersection.intersectedTriangle.normal));
//           float specular = pow(glm::dot(glm::normalize((-rayDirection)), reflection), 128);

//           specularCount += specular;

//           float shadowBias = 1e-4;

//           //check if the ray is in shadow. 
//           for(ModelTriangle triangle : model.tris) {
//             vec4 e0 = triangle.vertices[1] - triangle.vertices[0];
//             vec4 e1 = triangle.vertices[2] - triangle.vertices[0];
//             vec4 SPVector = (intersection.intersectionPoint) - triangle.vertices[0];
//             mat4 DEMatrix(-shadowRayDirection, e0, e1, vec4(1,1,1,1));
//             vec4 possibleSolution = glm::inverse(DEMatrix) * SPVector;

//             //check if ray intersects triangle and not just triangle plane
//             if(possibleSolution.y >= 0 && possibleSolution.y <= 1 && possibleSolution.z >= 0 && possibleSolution.z <= 1 && possibleSolution.y + possibleSolution.z <= 1) {
//               //I genuinely have no idea why doing the shadowBias this way works. without it the shadows are just everywhere???
//               if(possibleSolution.x < 1 - shadowBias && possibleSolution.x > shadowBias) {
//                 inShadow = true;
//                 break;
//               }
//             }
//           }

//           shadowCount -= inShadow ? shadow/10 : 0;
//         }

//         //adjust the totalled lighting values
//         shadowCount = glm::clamp<float>(shadowCount, 0, 1);
//         brightnessCount = glm::clamp<float>(brightnessCount, 0, 1);
//         angleCount = glm::clamp<float>(angleCount, 0, 1);
//         specularCount = glm::clamp<float>(specularCount, 0, 1);

//         Colour c = addColours(intersection.intersectedTriangle.material.diffuse, transparentColour);

//         //set the final pixels
//         if(intersection.intersectedTriangle.name == "light") window.setPixelColour(i, j, c.toPackedInt());
//         else window.setPixelColour(i, j, darkenColour(c, glm::clamp<float>(angleCount * (inShadow ? shadowCount : 1.0f) * brightnessCount, AMBIENCE, 1), inShadow || mirror ? 0 : specularCount));
//       } else {
//         window.setPixelColour(i, j, 0);
//       }
//     }
//   }
// }

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
  sphere.setPosition(vec3(0,10.0f,-3));
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
  auto start = std::chrono::high_resolution_clock::now();
  int frameCount = 0;
  while (true)
  {
    auto elapsed = std::chrono::high_resolution_clock::now() - start;
    long long millis = std::chrono::duration_cast<std::chrono::milliseconds>(elapsed).count();
    if (millis >= 1000) {
      cout << "FPS: " << frameCount << '\n';
      start = std::chrono::high_resolution_clock::now();
      frameCount = 0;
    }
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
      draw();
      raytrace(cam, renderQueue);
    } else {
      draw();
      drawTriangles(cam, renderQueue);
      if (SSAA) downsample(imageBuffer, window.pixelBuffer, WIDTH, HEIGHT, SSAA_SAMPLES);
    }
    // Need to render the frame at the end, or nothing actually gets shown on
    // the screen !
    window.renderFrame();
    frameCount++;
  }
}

void draw()
{
  if (!SSAA) window.clearPixels();
  size_t img_size = SSAA ? (IMG_SIZE * SSAA_SAMPLES) : IMG_SIZE;
  for (size_t i = 0; i < img_size; i++)
  {
    depthBuffer[i] = std::numeric_limits<float>::infinity();
  }
  memset(imageBuffer, 0, img_size * sizeof(uint32_t));
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
    else if (event.key.keysym.sym == SDLK_b) {
      cout << "B" << endl;
      bilinear = !bilinear;
    }
    else if (event.key.keysym.sym == SDLK_p) {
      cout << "P" << endl;
      perspective = !perspective;
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

// Line clipping based on https://www.geeksforgeeks.org/line-clipping-set-1-cohen-sutherland-algorithm/
int clipCode(vec4& p, int width, int height) {
  int code = 0;
  if (p.x < 0.0f) code |= LEFT;
  if (p.y < 0.0f) code |= TOP;
  if (p.x > width - 1.0f) code |= RIGHT;
  if (p.y > height - 1.0f) code |= BOTTOM;
  return code;
}

bool clipLine(vec4& p, vec4& q, int width, int height) {
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

void line(vec4 p, vec4 q, int colour, uint32_t *buffer, vec2 &offset)
{
  if (!clipLine(p, q, WIDTH, HEIGHT)) return;
  float x_diff = p.x - q.x;
  float y_diff = p.y - q.y;
  float max_diff = std::max(abs(x_diff), abs(y_diff));
  vector<float> interpolate_x = Interpolate(p.x + offset.x, q.x + offset.x, max_diff + 1);
  vector<float> interpolate_y = Interpolate(p.y + offset.y, q.y + offset.y, max_diff + 1);
  for (int i = 0; i < max_diff; i++)
  {
    buffer[(int)(interpolate_y[i]) * WIDTH + (int)interpolate_x[i]] = colour;
  }
}

inline float edgeFunction(const float v0_x, const float v0_y, const float v1_x, const float v1_y, const float p_x, const float p_y) {
  return (p_x - v0_x) * (v1_y - v0_y) - (p_y - v0_y) * (v1_x - v0_x);
}

inline int scaleColour(int colour, float scale) {
  unsigned char red = (colour & RED) >> 16;
  red *= scale;
  unsigned char green = (colour & GREEN) >> 8;
  green *= scale;
  unsigned char blue = (colour & BLUE);
  blue *= scale;
  return (colour & ALPHA) | (red << 16) | (green << 8) | blue;
}

inline int bilinearColour(int tl, int tr, int bl, int br, vec2 pos) {
  float xy = pos.x * pos.y;
  float a0 = xy - pos.x - pos.y + 1.0f;
  float a1 = pos.y - xy;
  float a2 = pos.x - xy;
  float a3 = xy;
  tl = scaleColour(tl, a0);
  bl = scaleColour(bl, a1);
  tr = scaleColour(tr, a2);
  br = scaleColour(br, a3);
  return ALPHA | (tl + bl + tr + br);
}

void triangle(Triangle &t, Texture &tex, int colour, bool filled, uint32_t *buffer, float *depthBuff, vec2 offset)
{
  if (filled)
  {
    bool textured = (t.vertices[0].u >= 0.0f);
    int x_min = glm::min(t.vertices[0].pos.x, t.vertices[1].pos.x);
    x_min = glm::min((float)x_min, t.vertices[2].pos.x);
    int x_max = glm::max(t.vertices[0].pos.x, t.vertices[1].pos.x);
    x_max = glm::max((float)x_max, t.vertices[2].pos.x);
    int y_min = glm::min(t.vertices[0].pos.y, t.vertices[1].pos.y);
    y_min = glm::min((float)y_min, t.vertices[2].pos.y);
    int y_max = glm::max(t.vertices[0].pos.y, t.vertices[1].pos.y);
    y_max = glm::max((float)y_max, t.vertices[2].pos.y);
    x_min = glm::clamp<int>(x_min, 0, WIDTH - 1);
    x_max = glm::clamp<int>(x_max, 0, WIDTH - 1);
    y_min = glm::clamp<int>(y_min, 0, HEIGHT - 1);
    y_max = glm::clamp<int>(y_max, 0, HEIGHT - 1);
    float area_inv = 1.0f / edgeFunction(t.vertices[0].pos.x, t.vertices[0].pos.y, t.vertices[1].pos.x, t.vertices[1].pos.y, t.vertices[2].pos.x, t.vertices[2].pos.y);
    float w0_step_x = (t.vertices[2].pos.y - t.vertices[1].pos.y) * area_inv;
    float w1_step_x = (t.vertices[0].pos.y - t.vertices[2].pos.y) * area_inv;
    float w2_step_x = (t.vertices[1].pos.y - t.vertices[0].pos.y) * area_inv;
    float w0_step_y = (t.vertices[1].pos.x - t.vertices[2].pos.x) * area_inv;
    float w1_step_y = (t.vertices[2].pos.x - t.vertices[0].pos.x) * area_inv;
    float w2_step_y = (t.vertices[0].pos.x - t.vertices[1].pos.x) * area_inv;
    vec2 p = vec2(x_min + offset.x, y_min + offset.y);
    float w0_line = edgeFunction(t.vertices[1].pos.x, t.vertices[1].pos.y, t.vertices[2].pos.x, t.vertices[2].pos.y, p.x, p.y) * area_inv;
    float w1_line = edgeFunction(t.vertices[2].pos.x, t.vertices[2].pos.y, t.vertices[0].pos.x, t.vertices[0].pos.y, p.x, p.y) * area_inv;
    float w2_line = edgeFunction(t.vertices[0].pos.x, t.vertices[0].pos.y, t.vertices[1].pos.x, t.vertices[1].pos.y, p.x, p.y) * area_inv;
    float w0, w1, w2;
    for (int y = y_min; y <= y_max; y++) {
      w0 = w0_line;
      w1 = w1_line;
      w2 = w2_line;
      for (int x = x_min; x <= x_max; x++) {
        if (w0 >= 0 && w1 >= 0 && w2 >= 0) {
          float depth = w0 * t.vertices[0].pos.z + w1 * t.vertices[1].pos.z + w2 * t.vertices[2].pos.z;
          if (depth < depthBuff[y * WIDTH + x]) {
            depthBuff[y * WIDTH + x] = depth;
            float q0, q1, q2;
            if (perspective) {
              q0 = w0 * t.vertices[0].pos.w;
              q1 = w1 * t.vertices[1].pos.w;
              q2 = w2 * t.vertices[2].pos.w;
            }
            else {
              q0 = w0; q1 = w1; q2 = w2;
            }
            float brightness = (q0 * t.vertices[0].brightness + q1 * t.vertices[1].brightness + q2 * t.vertices[2].brightness) / (q0 + q1 + q2);
            if (textured) {
              float u = (q0 * t.vertices[0].u + q1 * t.vertices[1].u + q2 * t.vertices[2].u) / (q0 + q1 + q2);
              float v = (q0 * t.vertices[0].v + q1 * t.vertices[1].v + q2 * t.vertices[2].v) / (q0 + q1 + q2);
              u = mod(u, 1.0f);
              v = mod(v, 1.0f);
              if (bilinear) {
                u *= tex.width - 1;
                v *= tex.height - 1;
                int tl = (int)u + (int)v * tex.width;
                int tr = tl + 1;
                int bl = tl + tex.width;
                int br = bl + 1;
                int biColour = bilinearColour(tex.data[tl], tex.data[tr], tex.data[bl], tex.data[br], vec2(mod(u, 1.0f), mod(v, 1.0f)));
                buffer[y * WIDTH + x] = scaleColour(biColour, brightness);
              }
              else {
                u *= tex.width;
                v *= tex.height;
                buffer[y * WIDTH + x] = scaleColour(tex.data[(int)u + (int)v * tex.width], brightness);
              }
            }
            else
              buffer[y * WIDTH + x] = scaleColour(colour, brightness);
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
    line(t.vertices[0].pos, t.vertices[1].pos, colour, buffer, offset);
    line(t.vertices[1].pos, t.vertices[2].pos, colour, buffer, offset);
    line(t.vertices[2].pos, t.vertices[0].pos, colour, buffer, offset);
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