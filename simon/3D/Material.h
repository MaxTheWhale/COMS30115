#pragma once

#include <iostream>

struct Texture {
  glm::vec3* dataVec;
  int width, height;
};

class Material
{
  public:
    std::string name;
    glm::vec3 ambientVec = glm::vec3(-1.0f, -1.0f, -1.0f);
    glm::vec3 diffuseVec;
    glm::vec3 specularVec = glm::vec3(1.0f, 1.0f, 1.0f);
    float highlights = 0.0f;
    int illum = 1;
    float dissolve = 1.0f;
    Texture texture = {nullptr, 0, 0};
    Texture normal_map = {nullptr, 0, 0};

    Material()
    {
    }

    Material(std::string n)
    {
        name = n;
    }

    Material(std::string n, glm::vec3 a, glm::vec3 d, glm::vec3 s)
    {
      name = n;
      ambientVec = a;
      diffuseVec = d;
      specularVec = s;
    }

};