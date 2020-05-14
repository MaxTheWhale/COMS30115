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
    Model(const Model& original);
    vec3 center;
    float furthestExtent = -1.0f;
    float calcExtent();
    vec4 calcTangent(ModelTriangle& tri);
  protected:
    vec3 centerOfMass();
    std::unordered_map<std::string, Material> loadMTL(std::string filename);
    std::vector<ModelTriangle> loadOBJ(std::string filename,
                              std::unordered_map<std::string, Material> palette);
};