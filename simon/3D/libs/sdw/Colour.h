#pragma once

#include <iostream>
#include <glm/glm.hpp>

class Colour
{
  public:
    std::string name;
    int red;
    int green;
    int blue;

    Colour()
    {
    }

    Colour(int r, int g, int b)
    {
      name = "";
      red = r;
      green = g;
      blue = b;
    }

    Colour(std::string n, int r, int g, int b)
    {
      name = n;
      red = r;
      green = g;
      blue = b;
    }

    Colour operator * (float x)
    {
      Colour temp = Colour();
      temp.red = glm::clamp<int>(red * x, 0, 255);
      temp.green = glm::clamp<int>(green * x, 0, 255);
      temp.blue = glm::clamp<int>(blue * x, 0, 255);

      return temp;
    }

    Colour operator * (Colour x)
    {
      Colour temp = Colour();
      temp.red = glm::clamp<int>(red * (x.red/255), 0, 255);
      temp.green = glm::clamp<int>(green * (x.green/255), 0, 255);
      temp.blue = glm::clamp<int>(blue * (x.blue/255), 0, 255);

      return temp;
    }

    Colour operator + (float x)
    {
      Colour temp = Colour();
      temp.red = glm::clamp<int>(red + x, 0, 255);
      temp.green = glm::clamp<int>(green + x, 0, 255);
      temp.blue = glm::clamp<int>(blue + x, 0, 255);

      return temp;
    }

    Colour operator + (Colour x)
    {
      Colour temp = Colour();
      temp.red = glm::clamp<int>(red + x.red, 0, 255);
      temp.green = glm::clamp<int>(green + x.green, 0, 255);
      temp.blue = glm::clamp<int>(blue + x.blue, 0, 255);

      return temp;
    }

    int toPackedInt() {
      return (0xff000000) | ((red & 0xff) << 16) | ((green & 0xff) << 8) | (blue & 0xff);
    }
};

std::ostream& operator<<(std::ostream& os, const Colour& colour);