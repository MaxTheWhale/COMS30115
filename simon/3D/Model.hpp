#pragma once

#include "Transformable.hpp"
#include "Animatable.hpp"
#include <ModelTriangle.h>
#include <vector>
#include <unordered_map>
#include <string>
#include "Material.h"

class Model : public Animatable {
  public:
    std::vector<ModelTriangle> tris;
    std::unordered_map<std::string, Material> palette;
    Model(std::string filename);
    vec3 center;
  protected:
    vec3 centerOfMass();
    float furthestExtent = -1.0f;
    float calcExtent();
    std::unordered_map<std::string, Material> loadMTL(std::string filename);
    std::vector<ModelTriangle> loadOBJ(std::string filename,
                              std::unordered_map<std::string, Material> palette);
};