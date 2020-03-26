#pragma once

#include <glm/glm.hpp>
#include "Colour.h"
#include <string>

class ModelTriangle
{
  public:
    glm::vec4 vertices[3];
    glm::vec4 normal;
    float brightness[3];
    Colour colour;
    std::string name;

    ModelTriangle()
    {
    }

    ModelTriangle(glm::vec4 v0, glm::vec4 v1, glm::vec4 v2, Colour trigColour)
    {
      vertices[0] = v0;
      vertices[1] = v1;
      vertices[2] = v2;
      colour = trigColour;
    }

    ModelTriangle(glm::vec4 v0, glm::vec4 v1, glm::vec4 v2, Colour trigColour, glm::vec4 trigNormal)
    {
      vertices[0] = v0;
      vertices[1] = v1;
      vertices[2] = v2;
      colour = trigColour;
      normal = trigNormal;
    }

    ModelTriangle(glm::vec4 v0, glm::vec4 v1, glm::vec4 v2, float b0, float b1, float b2, Colour trigColour, glm::vec4 trigNormal)
    {
      vertices[0] = v0;
      vertices[1] = v1;
      vertices[2] = v2;
      brightness[0] = b0;
      brightness[1] = b1;
      brightness[2] = b2;
      colour = trigColour;
      normal = trigNormal;
    }

    ModelTriangle(glm::vec3 v0, glm::vec3 v1, glm::vec3 v2, Colour trigColour)
    {
      vertices[0] = glm::vec4(v0.x, v0.y, v0.z, 1);
      vertices[1] = glm::vec4(v1.x, v1.y, v1.z, 1);
      vertices[2] = glm::vec4(v2.x, v2.y, v2.z, 1);
      normal = glm::vec4(glm::normalize(glm::cross(v0 - v1, v0 - v2)), 0);
      colour = trigColour;
    }

    ModelTriangle(glm::vec3 v0, glm::vec3 v1, glm::vec3 v2, Colour trigColour, std::string n)
    {
      vertices[0] = glm::vec4(v0.x, v0.y, v0.z, 1);
      vertices[1] = glm::vec4(v1.x, v1.y, v1.z, 1);
      vertices[2] = glm::vec4(v2.x, v2.y, v2.z, 1);
      normal = glm::vec4(glm::normalize(glm::cross(v0 - v1, v0 - v2)), 0);
      colour = trigColour;
      name = n;
    }
};

std::ostream& operator<<(std::ostream& os, const ModelTriangle& triangle);
