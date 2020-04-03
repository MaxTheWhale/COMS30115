#include "TexturePoint.h"
#include <iostream>

class CanvasPoint
{
  public:
    float x;
    float y;
    float depth;
    float brightness;

    CanvasPoint()
    {
    }

    CanvasPoint(float xPos, float yPos)
    {
      x = xPos;
      y = yPos;
      depth = 0.0;
      brightness = 1.0;
    }

    CanvasPoint(float xPos, float yPos, float pointDepth)
    {
      x = xPos;
      y = yPos;
      depth = pointDepth;
      brightness = 1.0;
    }

    CanvasPoint(float xPos, float yPos, float pointDepth, float pointBrightness)
    {
      x = xPos;
      y = yPos;
      depth = pointDepth;
      brightness = pointBrightness;
    }

    CanvasPoint operator+=(const CanvasPoint& rhs)
    {
      x += rhs.x;
      y += rhs.y;
      depth += rhs.depth;
      brightness += rhs.brightness;
      return *this;
    }

    friend CanvasPoint operator+(CanvasPoint lhs, const CanvasPoint& rhs)
    {
      lhs += rhs;
      return lhs;
    }

    CanvasPoint operator-=(const CanvasPoint& rhs)
    {
      x -= rhs.x;
      y -= rhs.y;
      depth -= rhs.depth;
      brightness -= rhs.brightness;
      return *this;
    }

    friend CanvasPoint operator-(CanvasPoint lhs, const CanvasPoint& rhs)
    {
      lhs -= rhs;
      return lhs;
    }

    CanvasPoint operator*=(float rhs)
    {
      x *= rhs;
      y *= rhs;
      depth *= rhs;
      brightness *= rhs;
      return *this;
    }

    friend CanvasPoint operator*(CanvasPoint lhs, float rhs)
    {
      lhs *= rhs;
      return lhs;
    }

    CanvasPoint operator/=(float rhs)
    {
      x /= rhs;
      y /= rhs;
      depth /= rhs;
      brightness /= rhs;
      return *this;
    }

    friend CanvasPoint operator/(CanvasPoint lhs, float rhs)
    {
      lhs /= rhs;
      return lhs;
    }
};

std::ostream& operator<<(std::ostream& os, const CanvasPoint& point)
{
    os << "(" << point.x << ", " << point.y << ", " << point.depth << ") " << point.brightness << std::endl;
    return os;
}
