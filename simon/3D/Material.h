#pragma once

#include <iostream>
#include <Colour.h>

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