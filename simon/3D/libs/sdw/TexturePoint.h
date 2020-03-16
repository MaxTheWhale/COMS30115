#pragma once

#include <iostream>

class TexturePoint
{
  public:
    float x;
    float y;

    TexturePoint()
    {
    }

    TexturePoint(float xPos, float yPos)
    {
      x = xPos;
      y = yPos;
    }

    TexturePoint operator+=(const TexturePoint& rhs)
    {
      x += rhs.x;
      y += rhs.y;
      return *this;
    }

    friend TexturePoint operator+(TexturePoint lhs, const TexturePoint& rhs)
    {
      lhs += rhs;
      return lhs;
    }

    TexturePoint operator-=(const TexturePoint& rhs)
    {
      x -= rhs.x;
      y -= rhs.y;
      return *this;
    }

    friend TexturePoint operator-(TexturePoint lhs, const TexturePoint& rhs)
    {
      lhs -= rhs;
      return lhs;
    }

    TexturePoint operator*=(float rhs)
    {
      x *= rhs;
      y *= rhs;
      return *this;
    }

    friend TexturePoint operator*(TexturePoint lhs, float rhs)
    {
      lhs *= rhs;
      return lhs;
    }

    TexturePoint operator/=(float rhs)
    {
      x /= rhs;
      y /= rhs;
      return *this;
    }

    friend TexturePoint operator/(TexturePoint lhs, float rhs)
    {
      lhs /= rhs;
      return lhs;
    }
};

std::ostream& operator<<(std::ostream& os, const TexturePoint& point)
{
    os << "(" << point.x << ", " << point.y << ")" << std::endl;
    return os;
}
