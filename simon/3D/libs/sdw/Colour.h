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
      temp.red = red * x;
      glm::clamp<int>(temp.red, 0, 255);
      temp.green = green * x;
      glm::clamp<int>(temp.green, 0, 255);
      temp.blue = blue * x;
      glm::clamp<int>(temp.blue, 0, 255);

      return temp;
    }

    int toPackedInt() {
      return (0xff000000) | ((red & 0xff) << 16) | ((green & 0xff) << 8) | (blue & 0xff);
    }
};

std::ostream& operator<<(std::ostream& os, const Colour& colour);