#pragma once

#include <iostream>
#include <Colour.h>

class Material
{
  public:
    std::string name;
    Colour ambient = Colour();
    Colour diffuse = Colour();
    Colour specular = Colour();
    int highlights = 0;
    int illum = 1;

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