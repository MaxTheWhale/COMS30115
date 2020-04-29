#pragma once

#include "Vertex.hpp"
#include "Material.h"
#include "ModelTriangle.h"

class Triangle {
  public:
    Vertex vertices[3];
    Material mat;
    glm::vec4 normal;
    Triangle(ModelTriangle &tri);
    Triangle(const Vertex &v0, const Vertex &v1, const Vertex &v2, const Material &tMat, const glm::vec4 &tNormal);
};