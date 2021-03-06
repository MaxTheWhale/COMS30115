#pragma once

#include <glm/glm.hpp>
#include "../../Material.h"
#include <string>

class ModelTriangle
{
  public:
    glm::vec4 vertices[3];
    glm::vec2 uvs[3];
    glm::vec4 normal;
    glm::vec4 tangent;
    glm::vec4 normals[3];
    glm::mat3 TBN;
    float brightness[3];
    bool fullBright = false;
    Material material;
    std::string name = "";

    ModelTriangle()
    {
    }

    ModelTriangle(glm::vec4 v0, glm::vec4 v1, glm::vec4 v2, Material trigMaterial)
    {
      vertices[0] = v0;
      vertices[1] = v1;
      vertices[2] = v2;
      material = trigMaterial;
      uvs[0] = glm::vec2(-1.0f, -1.0f);
      uvs[1] = glm::vec2(-1.0f, -1.0f);
      uvs[2] = glm::vec2(-1.0f, -1.0f);
    }

    ModelTriangle(glm::vec4 v0, glm::vec4 v1, glm::vec4 v2, Material trigMaterial, glm::vec4 trigNormal)
    {
      vertices[0] = v0;
      vertices[1] = v1;
      vertices[2] = v2;
      material = trigMaterial;
      normal = trigNormal;
      uvs[0] = glm::vec2(-1.0f, -1.0f);
      uvs[1] = glm::vec2(-1.0f, -1.0f);
      uvs[2] = glm::vec2(-1.0f, -1.0f);
      normals[0] = normal;
      normals[1] = normal;
      normals[2] = normal;
    }

    ModelTriangle(glm::vec4 v0, glm::vec4 v1, glm::vec4 v2, float b0, float b1, float b2, Material trigMaterial, glm::vec4 trigNormal)
    {
      vertices[0] = v0;
      vertices[1] = v1;
      vertices[2] = v2;
      brightness[0] = b0;
      brightness[1] = b1;
      brightness[2] = b2;
      material = trigMaterial;
      normal = trigNormal;
      uvs[0] = glm::vec2(-1.0f, -1.0f);
      uvs[1] = glm::vec2(-1.0f, -1.0f);
      uvs[2] = glm::vec2(-1.0f, -1.0f);
      normals[0] = normal;
      normals[1] = normal;
      normals[2] = normal;
    }

    ModelTriangle(glm::vec3 v0, glm::vec3 v1, glm::vec3 v2, Material trigMaterial)
    {
      vertices[0] = glm::vec4(v0.x, v0.y, v0.z, 1);
      vertices[1] = glm::vec4(v1.x, v1.y, v1.z, 1);
      vertices[2] = glm::vec4(v2.x, v2.y, v2.z, 1);
      normal = glm::vec4(glm::normalize(glm::cross(v0 - v1, v0 - v2)), 0);
      material = trigMaterial;
      uvs[0] = glm::vec2(-1.0f, -1.0f);
      uvs[1] = glm::vec2(-1.0f, -1.0f);
      uvs[2] = glm::vec2(-1.0f, -1.0f);
      normals[0] = normal;
      normals[1] = normal;
      normals[2] = normal;
    }

    ModelTriangle(glm::vec3 v0, glm::vec3 v1, glm::vec3 v2, Material trigMaterial, std::string n)
    {
      vertices[0] = glm::vec4(v0.x, v0.y, v0.z, 1);
      vertices[1] = glm::vec4(v1.x, v1.y, v1.z, 1);
      vertices[2] = glm::vec4(v2.x, v2.y, v2.z, 1);
      normal = glm::vec4(glm::normalize(glm::cross(v0 - v1, v0 - v2)), 0);
      normals[0] = normal;
      normals[1] = normal;
      normals[2] = normal;
      material = trigMaterial;
      name = n;
      uvs[0] = glm::vec2(-1.0f, -1.0f);
      uvs[1] = glm::vec2(-1.0f, -1.0f);
      uvs[2] = glm::vec2(-1.0f, -1.0f);
    }
};

std::ostream& operator<<(std::ostream& os, const ModelTriangle& triangle);
