#pragma once

#include <iostream>
#include <Colour.h>

struct Texture {
  int* data;
  glm::vec3* dataVec;
  int width, height;
};

class Material
{
  public:
    std::string name;
    Colour ambient = Colour(-1, -1, -1);
    Colour diffuse = Colour(-1, -1, -1);
    Colour specular = Colour(-1, -1, -1);
    glm::vec3 ambientVec;
    glm::vec3 diffuseVec;
    glm::vec3 specularVec = glm::vec3(1.0f, 1.0f, 1.0f);
    float highlights = 0.0f;
    int illum = 1;
    int dissolve = 1;
    Texture texture = {nullptr, 0, 0};

    Material()
    {
    }

    Material(std::string n)
    {
        name = n;
    }

    Material(std::string n, Colour a, Colour d, Colour s)
    {
      name = n;
      ambient = a;
      diffuse = d;
      specular = d;
    }

};

std::ostream& operator<<(std::ostream& os, const Colour& diffuse);