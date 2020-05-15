#pragma once

#include "Transformable.hpp"
#include "Animatable.hpp"
#include <ModelTriangle.h>
#include <vector>
#include <unordered_map>
#include <string>
#include "Material.h"

//class for loading and storing OBJ files as distinct models in the scene
class Model : public Animatable {
  public:
    std::vector<ModelTriangle> tris;
    std::unordered_map<std::string, Material> palette;
    Model(std::string filename);
    Model(const Model& original);
    vec3 center;
    float furthestExtent = -1.0f;
    bool castShadow = true; //does this model interrupt lighting
    bool fullBright = false; //if true the model is rendered at maximum brightness at all times
    float calcExtent();
    vec4 calcTangent(ModelTriangle& tri);
  protected:
    vec3 centerOfMass();
    std::unordered_map<std::string, Material> loadMTL(std::string filename);
    std::vector<ModelTriangle> loadOBJ(std::string filename,
                              std::unordered_map<std::string, Material> palette);
};

glm::vec3 *loadPPM(string fileName, int &width, int &height);