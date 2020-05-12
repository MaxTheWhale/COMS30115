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
    float dissolve = 1.0f;
    Texture texture = {nullptr, nullptr, 0, 0};
    Texture normal_map = {nullptr, nullptr, 0, 0};

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
      specular = s;
      ambientVec = glm::vec3(a.red / 255.0f, a.green / 255.0f, a.blue / 255.0f);
      diffuseVec = glm::vec3(d.red / 255.0f, d.green / 255.0f, d.blue / 255.0f);
      specularVec = glm::vec3(s.red / 255.0f, s.green / 255.0f, s.blue / 255.0f);
    }

};

std::ostream& operator<<(std::ostream& os, const Colour& diffuse);