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
    Texture texture;
    Model(std::string filename);
  protected:
    std::unordered_map<std::string, Material> loadMTL(std::string filename, int*& data, int& width, int& height);
    std::vector<ModelTriangle> loadOBJ(std::string filename,
                              std::unordered_map<std::string, Material> palette);
};