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
#include "SSAA.hpp"
#include "Vertex.hpp"
#include "Triangle.hpp"
#include "Clipping.hpp"
#include "Lines.hpp"
#include "Camera.hpp"
#include "Model.hpp"
#include "Times.hpp"
#include "VectorUtil.hpp"
#include "Rigidbody.hpp"
#include "Light.h"
#include "Magnet.hpp"

using namespace std;
using namespace glm;

#define WIDTH 640
#define HEIGHT 480
#define IMG_SIZE (WIDTH*HEIGHT)
#define SSAA true
#define SSAA_SCALE 2
#define SSAA_SAMPLES (SSAA_SCALE*SSAA_SCALE)
#define MOUSE_SENSITIVITY 0.0015f
#define AMBIENCE 0.1f
#define ASPECT_RATIO WIDTH/(float)HEIGHT

#ifndef M_PIf
#define M_PIf 3.14159265358979323846f
#endif

enum COLOUR_MASK {ALPHA = 0xff000000, RED = 0x00ff0000, GREEN = 0x0000ff00, BLUE = 0x000000ff};

void draw();
void triangle(Triangle &t, bool filled, uint32_t *buffer, float *depthBuff, vec2 offset, vec4 &eye_pos);
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
bool bilinear = true;
bool perspective = true;
bool toRaytrace = false;
bool softShadows = false;
DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);

vec3 bilinearInterpolate(vec3 tl, vec3 tr, vec3 bl, vec3 br, vec2 pos) {
  float xy = pos.x * pos.y;
  float a0 = xy - pos.x - pos.y + 1.0f;
  float a1 = pos.y - xy;
  float a2 = pos.x - xy;
  float a3 = xy;
  tl *= a0;
  bl *= a1;
  tr *= a2;
  br *= a3;
  return (tl + bl + tr + br);
}

vec3 getTexPoint(float u, float v, Texture& tex, bool bilinear) {
  if (bilinear) {
    float u_bi = u * (tex.width - 1);
    float v_bi = v * (tex.height - 1);
    int tl = (int)u_bi + (int)v_bi * tex.width;
    int tr = tl + 1;
    int bl = tl + tex.width;
    int br = bl + 1;
    return bilinearInterpolate(tex.dataVec[tl], tex.dataVec[tr], tex.dataVec[bl], tex.dataVec[br], vec2(mod(u_bi, 1.0f), mod(v_bi, 1.0f)));
  }
  else {
    return tex.dataVec[(int)(u * tex.width) + (int)(v * tex.height) * tex.width];
  }
}

void drawTriangles(Camera &cam, std::vector<Model *> models)
{
  uint32_t *buffer = (SSAA) ? imageBuffer : window.pixelBuffer;
  vector<vec2> offsets = generateRotatedGrid(SSAA_SCALE);
  mat4 viewProjection = cam.projection * cam.worldToCamera();
  vec4 eye = vec4(cam.getPosition(), 1.0f);
  for (unsigned int i = 0; i < models.size(); i++)
  {
    Model &model = *models[i];
    for (auto& modelTri : model.tris)
    {
      Triangle tri = Triangle(modelTri);
      tri.normal = normalize(model.transform * tri.normal);
      tri.vertices[0].normal = normalize(model.transform * tri.vertices[0].normal);
      tri.vertices[1].normal = normalize(model.transform * tri.vertices[1].normal);
      tri.vertices[2].normal = normalize(model.transform * tri.vertices[2].normal);
      tri.vertices[0].pos_3d = model.transform * tri.vertices[0].pos;
      if (dot(tri.vertices[0].pos_3d - eye, tri.normal) >= 0.0f) continue;
      tri.vertices[1].pos_3d = model.transform * tri.vertices[1].pos;
      tri.vertices[2].pos_3d = model.transform * tri.vertices[2].pos;
      if (tri.mat.normal_map.dataVec != nullptr) {
        tri.tangent = model.transform * tri.tangent;
        tri.TBN = mat3(toThree(tri.tangent), toThree(cross(tri.normal, tri.tangent)), toThree(tri.normal));
      }
      tri.vertices[0].pos = viewProjection * tri.vertices[0].pos_3d;
      tri.vertices[1].pos = viewProjection * tri.vertices[1].pos_3d;
      tri.vertices[2].pos = viewProjection * tri.vertices[2].pos_3d;
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
          t.vertices[v].pos.z = ((cam.far - cam.near) * 0.5f) * t.vertices[v].pos.z + ((cam.far + cam.near) * 0.5f);
        }
        #if SSAA
        #pragma omp parallel for
        for (int s = 0; s < SSAA_SAMPLES; s++) {
          triangle(t, wireframe, buffer + (IMG_SIZE * s), depthBuffer + (IMG_SIZE * s), offsets[s], eye);
        }
        #else
        triangle(t, wireframe, buffer, depthBuffer, vec2(0.5f, 0.5f), eye);
        #endif
      }
    }
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

RayTriangleIntersection findClosestIntersection(vec4 start, Model model, vec4 rayDirection) {
  RayTriangleIntersection intersection = RayTriangleIntersection();
  float minDistance = std::numeric_limits<float>::infinity();

  //calculate closest intersection by looping through each of the triangles
  for(ModelTriangle triangle : model.tris) {
    vec4 e0 = triangle.vertices[1] - triangle.vertices[0];
    vec4 e1 = triangle.vertices[2] - triangle.vertices[0];
    vec4 SPVector = start - triangle.vertices[0];
    mat4 DEMatrix(-rayDirection, e0, e1, vec4(1, 1, 1, 1));
    vec4 possibleSolution = glm::inverse(DEMatrix) * SPVector;

    // check if ray intersects triangle and not just triangle plane
    if (possibleSolution.y >= 0 && possibleSolution.y <= 1 &&
      possibleSolution.z >= 0 && possibleSolution.z <= 1 &&
      possibleSolution.y + possibleSolution.z <= 1) {
      if (possibleSolution.x < minDistance && possibleSolution.x > 0) {
        intersection = RayTriangleIntersection(start + (possibleSolution.x * rayDirection) , possibleSolution.x, triangle);
        intersection.intersectionPoint.w = 1;
        intersection.u = possibleSolution.y;
        intersection.v = possibleSolution.z;
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

vec4 refract(vec4 I, vec4 N, float ior) {
  vec4 refractedRay;

  float cosi = glm::clamp<float>(glm::dot(I, N), -1.0f, 1.0f);
  float etai = 1.0f, etat = ior;
  vec4 n = N;
  if (cosi < 0) {
    cosi = -cosi;
  } else {
    std::swap(etai, etat);
    n = -N;
  }
  float eta = etai / etat;
  float k = 1 - eta * eta * (1 - cosi * cosi);
  return k < 0 ? vec4(0, 0, 0, 0) : eta * I + (eta * cosi - sqrtf(k)) * n;
}

Colour getPixelColour(RayTriangleIntersection intersection, Light mainLight, vec4 rayDirection, Model model, int depth, int i, int j) {
  Texture tex = intersection.intersectedTriangle.material.texture;
  Texture bumps = intersection.intersectedTriangle.material.normal_map;
  vec4 normal;

  if(bumps.dataVec != nullptr) {
          ModelTriangle t = intersection.intersectedTriangle;
          float q0 = intersection.u;
          float q1 = intersection.v;
          float q2 = 1 - q0 - q1;
          float u = q2 * t.uvs[0].x + q0 * t.uvs[1].x + q1 * t.uvs[2].x;
          float v = q2 * t.uvs[0].y + q0 * t.uvs[1].y + q1 * t.uvs[2].y;
          u = mod(u, 1.0f);
          v = mod(v, 1.0f);
          vec3 bumpVec = getTexPoint(u, v, bumps, bilinear);
          normal = vec4(bumpVec.x, bumpVec.y, bumpVec.z, 0.0f);
        } else {
          normal = intersection.intersectedTriangle.normal;
        }

  Colour colour = Colour(0, 0, 0);

      if(intersection.intersectedTriangle.material.dissolve < 1 && depth < 2) {
        // cout << "trying with normal " << normal.x << ", " << normal.y << ", " << normal.z << endl;
        vec4 refractedRay;
        float kr;

        float cosi = glm::clamp<float>(glm::dot(rayDirection, normal), -1.0f, 1.0f);
        float etai = 1.0f, etat = 1.5f;

        if(cosi > 0) {
          std::swap(etai, etat);
        }

        float sint = etai / etat * sqrtf(std::max(0.0f, 1 - cosi * cosi));

        if (sint >= 1) {
          kr = 1.0f;
        } else {
          float cost = sqrtf(std::max(0.0f, 1 - sint * sint));
          cosi = fabsf(cosi);
          float Rs = ((etat * cosi) - (etai * cost)) / ((etat * cosi) + (etai * cost));
          float Rp = ((etai * cosi) - (etat * cost)) / ((etai * cosi) + (etat * cost));
          kr = (Rs * Rs + Rp * Rp) / 2;
        }

        bool outside = glm::dot(rayDirection, normal) < 0;
        float bias = 1e-4;

        Colour refractionColour = intersection.intersectedTriangle.material.diffuse;
        if (kr < 1) {
          vec4 refractionDirection = refract(rayDirection, normal, 1.5f);
          RayTriangleIntersection refractIntersection = findClosestIntersection(outside ? intersection.intersectionPoint - (normal * bias) : intersection.intersectionPoint + (normal * bias), model, refractionDirection);
          refractionColour = getPixelColour(refractIntersection, mainLight, rayDirection, model, depth + 1, i, j);
        }

        Colour reflectedColour = intersection.intersectedTriangle.material.diffuse;
        if(depth < 2) {
          vec4 mirrorRayDirection = glm::normalize(rayDirection - 2.0f * (glm::dot(rayDirection, normal) * normal));
          mirrorRayDirection.w = 0;

          RayTriangleIntersection mirrorIntersection = findClosestIntersection(outside ? intersection.intersectionPoint - (normal * bias) : intersection.intersectionPoint + (normal * bias), model, mirrorRayDirection);
          
          reflectedColour = getPixelColour(mirrorIntersection, mainLight, rayDirection, model, depth + 1, i, j) + (intersection.intersectedTriangle.material.specular * 0.8f);
        } else {
          reflectedColour = intersection.intersectedTriangle.material.specular;
        }

        Colour glassColour;
        if(tex.dataVec != nullptr) {
          ModelTriangle t = intersection.intersectedTriangle;
          float q0 = intersection.u;
          float q1 = intersection.v;
          float q2 = 1 - q0 - q1;
          float u = q2 * t.uvs[0].x + q0 * t.uvs[1].x + q1 * t.uvs[2].x;
          float v = q2 * t.uvs[0].y + q0 * t.uvs[1].y + q1 * t.uvs[2].y;
          u = mod(u, 1.0f);
          v = mod(v, 1.0f);
          vec3 texVec = getTexPoint(u, v, tex, bilinear);
          glassColour = Colour(texVec.r * 255, texVec.g * 255, texVec.b * 255);
        } else {
          glassColour = intersection.intersectedTriangle.material.diffuse;
        }
        
        return glassColour * intersection.intersectedTriangle.material.dissolve + reflectedColour * kr + refractionColour * (1 - kr);
      }

      if(intersection.intersectedTriangle.material.specular.red >= 0) {
        if(depth < 2) {
          vec4 mirrorRayDirection = glm::normalize(rayDirection - 2.0f * (glm::dot(rayDirection, normal) * normal));
          mirrorRayDirection.w = 0;

          RayTriangleIntersection mirrorIntersection = findClosestIntersection(intersection.intersectionPoint + (normal * 0.1f), model, mirrorRayDirection);
          
          return getPixelColour(mirrorIntersection, mainLight, rayDirection, model, depth + 1, i, j) + (intersection.intersectedTriangle.material.specular * 0.8f);
        } else {
          return intersection.intersectedTriangle.material.specular;
        }
      }

      if(intersection.intersectedTriangle.name == mainLight.name) {
        return intersection.intersectedTriangle.material.diffuse;
      } else {
        if(tex.dataVec != nullptr) {
          ModelTriangle t = intersection.intersectedTriangle;
          float q0 = intersection.u;
          float q1 = intersection.v;
          float q2 = 1 - q0 - q1;
          float u = q2 * t.uvs[0].x + q0 * t.uvs[1].x + q1 * t.uvs[2].x;
          float v = q2 * t.uvs[0].y + q0 * t.uvs[1].y + q1 * t.uvs[2].y;
          u = mod(u, 1.0f);
          v = mod(v, 1.0f);
          vec3 texVec = getTexPoint(u, v, tex, bilinear);
          colour = Colour(texVec.r * 255, texVec.g * 255, texVec.b * 255);
        } else {
          colour = intersection.intersectedTriangle.material.diffuse;
        }

        vec4 shadowRayDirection = mainLight.centre - intersection.intersectionPoint;
        bool isInShadow = inShadow(model, shadowRayDirection, intersection);

        if(isInShadow) colour = colour * mainLight.shadow;
        else {
          //calculate the angleOfIncidence between 0 and 1
          float angleOfIncidence = glm::dot(glm::normalize(shadowRayDirection), normal);
          colour = colour * glm::clamp<float>(angleOfIncidence, mainLight.shadow, 1);

          //adjust brightness for proximity lighting
          float brightness = mainLight.intensity/pow(length(shadowRayDirection),2);
          colour = colour * glm::clamp<float>(brightness, mainLight.shadow, 1);

          if(intersection.intersectedTriangle.material.highlights > 0) {
            //TODO: make the colour of the highlights match the specular material colour
            vec4 reflection = glm::normalize((-shadowRayDirection) - 2.0f * (glm::dot((-shadowRayDirection), normal) * normal));
            float specular = pow(glm::dot(glm::normalize((-rayDirection)), reflection), intersection.intersectedTriangle.material.highlights);
            colour = colour + (glm::clamp<float>(specular, 0, 1) * 255);
          }
        }

        return colour; 
      }
      return colour;
}

void raytrace(Camera camera, std::vector<Model*> models) {
  //idk what this does
  Model model = Model(*models[0]);
  for (unsigned int i = 1; i < models.size(); i++) {
    for (auto tri : (*models[i]).tris) {
      ModelTriangle newTri = ModelTriangle((*models[i]).transform * tri.vertices[0],
                            (*models[i]).transform * tri.vertices[1],
                            (*models[i]).transform * tri.vertices[2],
                            tri.material, (*models[i]).transform * tri.normal);
      newTri.uvs[0] = tri.uvs[0];
      newTri.uvs[1] = tri.uvs[1];
      newTri.uvs[2] = tri.uvs[2];
      if (newTri.material.normal_map.dataVec != nullptr) {
        newTri.tangent = (*models[i]).transform * tri.tangent;
        newTri.TBN = mat3(toThree(newTri.tangent), toThree(cross(newTri.normal, newTri.tangent)), toThree(newTri.normal));
      }
      model.tris.push_back(newTri);
    }
  }

  //setup the lights
  vector<ModelTriangle> lights = vector<ModelTriangle>();
  for(ModelTriangle triangle : model.tris) {
    if(triangle.name == "light") lights.push_back(triangle);
  }

  Light mainLight = Light("light", lights);
  mainLight.calculateCentre();

  //loop through each pixel in image plane
  #pragma omp parallel for
  for(int j = 0; j < HEIGHT; j++) {
    for(int i = 0; i < WIDTH; i++) {
      float angle = tanf(0.5f * glm::radians(camera.fov)); // just fov*0.5 converted to radians
      //convert image plane cordinates into world space
      vec2 NDC = vec2((i + 0.5) * (1 / (float) WIDTH), (j + 0.5) * (1 / (float) HEIGHT));
      float x = (2 * (NDC.x) - 1) * angle * ASPECT_RATIO;
      float y = (1 - 2 * (NDC.y)) * angle;

      // the main camera ray
      vec4 rayDirection = camera.transform * vec4(x, y, -1, 0);

      RayTriangleIntersection intersection = findClosestIntersection(camera.transform[3], model, rayDirection);

      window.setPixelColour(i, j, getPixelColour(intersection, mainLight, rayDirection, model, 0, i, j).toPackedInt());  
    }
  }
  savePPM("window.ppm", &window);
  cout << "Finished one frame!" << endl;
}

//vector<mat4> cameraTransforms = vector<mat4>();

int main(int argc, char *argv[])
{
  SDL_Event event;
  //SDL_SetRelativeMouseMode(SDL_TRUE);

  vector<Model*> renderQueue = vector<Model*>();
  vector<Updatable*> updateQueue = vector<Updatable*>();

  Model cornell = Model("cornell-box");
  //cornell.rotate(glm::vec3(45,0,0));
  renderQueue.push_back(&cornell);
  // std::cout << "cornell address = " << &cornell << std::endl;
  // Rigidbody cornellRB = Rigidbody(&cornell);
  // cornellRB.hasGravity = false;
  // updateQueue.push_back(&cornellRB);

  // Model hs_logo = Model("HackspaceLogo/logo");
  // hs_logo.scale(vec3(0.005f, 0.005f, 0.005f));
  // hs_logo.move(vec3(-1.1f, 1.21f, -1.8f));
  // renderQueue.push_back(&hs_logo);

  // Model glass = Model("GlassSphere/glassSphere");
  // glass.scale(vec3(0.75f, 0.75f, 0.75f));
  // glass.move(vec3(0.75f, 2.3f, -2.0f));
  // renderQueue.push_back(&glass);

  // for (int i = -5; i < 5; i++) {
    // Model logo = Model("HackspaceLogo/Logo");
    // logo.scale(vec3(0.005f, 0.005f, 0.005f));
    // logo.move(vec3(-1, 1, 0));
    // logo.center = vec3(300, 300, 0);
    // renderQueue.push_back(&logo);
    // Rigidbody rb = Rigidbody(&logo);
    // rb.hasGravity = false;
    // rb.collisionEnabled = false;
    // rb.applyForce(vec3(-0.2f, 0, 0), vec3(0,-10,0));
    // updateQueue.push_back(&rb);
  // }

  // Model miku = Model("sc");
  // miku.rotate(vec3(0.0f, 0.3f, 0.0f));
  // miku.move(vec3(1.0f, 0.0f, -4.5f));
  // renderQueue.push_back(&miku);

  Model bumpy = Model("bumpy");
  bumpy.scale(vec3(0.75f, 0.75f, 0.75f));
  bumpy.move(vec3(0.75f, 2.3f, -2.0f));
  renderQueue.push_back(&bumpy);

  // Model cornell2 = Model("cornell-box");
  // // cornell2.move(glm::vec3(0,1,0));
  // cornell2.rotate(glm::vec3(0,1,0));
  // renderQueue.push_back(&cornell2);
  // Rigidbody cornellRB2 = Rigidbody(&cornell2);
  // cornellRB2.hasGravity = false;
  // updateQueue.push_back(&cornellRB2);

  Model sphere = Model("blob");
  sphere.setPosition(vec3(3,0.0f,0));
  //sphere.rotate(vec3(1,0,0));
  renderQueue.push_back(&sphere);
  Rigidbody sphereRB = Rigidbody(&sphere);
  updateQueue.push_back(&sphereRB);
  sphereRB.positionFixed = false;
  sphereRB.hasGravity = false;
  sphereRB.collisionEnabled = true;
  sphereRB.elasticity = 0.9f;
  // sphereRB.applyForce(vec3(0,0.05f,0), vec3(0.5f,0,0));
  sphereRB.applyForce(vec3(0.5f,0,0), vec3(0,0,0));
  cout << "sphereRB address = " << &sphereRB << endl;

  Model angle = Model("tilted");
  angle.scale(vec3(3,3,3));
  angle.setPosition(vec3(0,0,0));
  angle.furthestExtent = angle.calcExtent();
  renderQueue.push_back(&angle);
  Magnet mag = Magnet(&angle);
  updateQueue.push_back(&mag);
  // Rigidbody angleRB = Rigidbody(&angle);
  // updateQueue.push_back(&angleRB);
  // angleRB.positionFixed = true;
  // angleRB.hasGravity = false;
  // angleRB.collisionEnabled = true;
  // cout << "angleRB address = " << &angleRB << endl;
  
  Model angle2 = Model("tilted");
  angle2.scale(vec3(3,3,3));
  angle2.setPosition(vec3(0,0,10));
  angle2.furthestExtent = angle2.calcExtent();
  renderQueue.push_back(&angle2);
  Magnet mag2 = Magnet(&angle2);
  updateQueue.push_back(&mag2);

  // cout << "normals: " << endl;
  // for (int i = 0; i < angle.tris.size(); i++) {
  //   cout << angle.tris[i].normal << endl;;
  // }

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
  // for (unsigned int i = 0; i < renderQueue.size(); i++) {
  //   for (auto tri : (*renderQueue[i]).tris) {
  //     cout << "UVs for " << tri.name << ": " << tri.uvs[0].x << "," << tri.uvs[0].y << "  " << tri.uvs[1].x << "," << tri.uvs[1].y << "  " << tri.uvs[2].x << "," << tri.uvs[2].y << endl;
  //   }
  // }
  Camera cam;
  cam.setProjection(90.0f, WIDTH / (float)HEIGHT, 0.1f, 100.0f);
  cam.lookAt(vec3(5.0f, 20.0f, 0.0f), vec3(0.0f, 0, 0));
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
      // raytraceOld(cam, renderQueue, 1);
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

inline float edgeFunction(const float v0_x, const float v0_y, const float v1_x, const float v1_y, const float p_x, const float p_y) {
  return (p_x - v0_x) * (v1_y - v0_y) - (p_y - v0_y) * (v1_x - v0_x);
}

inline int vec3ToPackedInt(vec3 colour) {
  return ALPHA | (int(colour.r * 255.0f) << 16) | (int(colour.g * 255.0f) << 8) | int(colour.b * 255.0f);
}

// As described here: https://en.wikipedia.org/wiki/Phong_reflection_model
inline vec3 phongReflection(vec3 &Ks, vec3 &Kd, vec3 &Ka, int &alpha, vec3 &Is, vec3 &Id, vec3 &Ia, vec3 &Lm, vec3 &N, vec3 &Rm, vec3 &V) {
  return (Kd * glm::max(dot(Lm, N), 0.0f) * Id) + (Ks * powf(dot(Rm, V), alpha) * Is) + Ka * Ia;
}

void triangle(Triangle &t, bool filled, uint32_t *buffer, float *depthBuff, vec2 offset, vec4 &eye_pos)
{
  vec4 light_pos = vec4(-0.234f, 5.2f, -3.043f, 1.0f);
  vec3 Ia = vec3(0.2f, 0.2f, 0.2f);
  vec3 Is = vec3(1.0f, 1.0f, 1.0f);
  vec3 Kd = t.mat.diffuseVec;
  vec3 Ks = t.mat.specularVec;
  vec3 Ka = t.mat.ambientVec;
  int alpha = t.mat.highlights;
  if (filled)
  {
    bool textured = (t.mat.texture.dataVec != nullptr);
    bool bump_map = (t.mat.normal_map.dataVec != nullptr);
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
              float q_inv = 1.0f / (q0 + q1 + q2);
              q0 *= q_inv;
              q1 *= q_inv;
              q2 *= q_inv;
            }
            else {
              q0 = w0; q1 = w1; q2 = w2;
            }
            float u, v;
            vec3 N;
            if (textured || bump_map) {
              u = q0 * t.vertices[0].u + q1 * t.vertices[1].u + q2 * t.vertices[2].u;
              v = q0 * t.vertices[0].v + q1 * t.vertices[1].v + q2 * t.vertices[2].v;
              u = mod(u, 1.0f);
              v = mod(v, 1.0f);
            }
            if (textured) {
              Kd = getTexPoint(u, v, t.mat.texture, bilinear);
              Ka = Kd;
            }
            if (bump_map) {
              N = getTexPoint(u, v, t.mat.normal_map, bilinear);
              N = t.TBN * N;
            }
            else {
              N = toThree(q0 * t.vertices[0].normal + q1 * t.vertices[1].normal + q2 * t.vertices[2].normal);
            }
            vec4 pos_3d = q0 * t.vertices[0].pos_3d + q1 * t.vertices[1].pos_3d + q2 * t.vertices[2].pos_3d;
            float radius = distance(light_pos, pos_3d);
            vec3 Id = vec3(200.0f, 200.0f, 200.0f) / (4.0f * M_PIf * radius * radius);
            vec3 V = toThree(normalize(eye_pos - pos_3d));
            vec3 Lm = toThree(normalize(light_pos - pos_3d));
            vec3 Rm = normalize(2.0f * N * dot(Lm, N) - Lm);
            if (t.mat.illum < 2) {
              Ks = vec3(0.0f);
            }
            vec3 reflectedLight = glm::min(phongReflection(Ks, Kd, Ka, alpha, Is, Id, Ia, Lm, N, Rm, V), 1.0f);
            buffer[y * WIDTH + x] = vec3ToPackedInt(reflectedLight);
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
    line(t.vertices[0].pos, t.vertices[1].pos, t.mat.diffuse.toPackedInt(), buffer, WIDTH, HEIGHT, offset);
    line(t.vertices[1].pos, t.vertices[2].pos, t.mat.diffuse.toPackedInt(), buffer, WIDTH, HEIGHT, offset);
    line(t.vertices[2].pos, t.vertices[0].pos, t.mat.diffuse.toPackedInt(), buffer, WIDTH, HEIGHT, offset);
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
      buffer[x * 3 + 2] = value & 0xff;
      buffer[x * 3 + 1] = (value & 0xff00) >> 8;
      buffer[x * 3] = (value & 0xff0000) >> 16;
    }
    f.write((char *)buffer, window->width * 3);
  }
  f.close();
}