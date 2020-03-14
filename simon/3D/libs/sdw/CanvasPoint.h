#include "TexturePoint.h"
#include <iostream>

class CanvasPoint
{
  public:
    float x;
    float y;
    float depth;
    float brightness;
    TexturePoint texturePoint;

    CanvasPoint()
    {
        texturePoint = TexturePoint(-1,-1);
    }

    CanvasPoint(float xPos, float yPos)
    {
      x = xPos;
      y = yPos;
      depth = 0.0;
      brightness = 1.0;
      texturePoint = TexturePoint(-1,-1);
    }

    CanvasPoint(float xPos, float yPos, float pointDepth)
    {
      x = xPos;
      y = yPos;
      depth = pointDepth;
      brightness = 1.0;
      texturePoint = TexturePoint(-1,-1);
    }

    CanvasPoint(float xPos, float yPos, float pointDepth, float pointBrightness)
    {
      x = xPos;
      y = yPos;
      depth = pointDepth;
      brightness = pointBrightness;
      texturePoint = TexturePoint(-1,-1);
    }

    CanvasPoint(float xPos, float yPos, float pointDepth, float pointBrightness, TexturePoint& texPoint)
    {
      x = xPos;
      y = yPos;
      depth = pointDepth;
      brightness = pointBrightness;
      texturePoint = texPoint;
    }

    CanvasPoint operator+=(const CanvasPoint& rhs)
    {
      x += rhs.x;
      y += rhs.y;
      depth += rhs.depth;
      brightness += rhs.brightness;
      texturePoint += rhs.texturePoint;
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
      texturePoint -= rhs.texturePoint;
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
      texturePoint *= rhs;
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
      texturePoint /= rhs;
      return *this;
    }

    friend CanvasPoint operator/(CanvasPoint lhs, float rhs)
    {
      lhs /= rhs;
      return lhs;
    }

    friend bool operator<(const CanvasPoint& lhs, const CanvasPoint& rhs)
    {
      return lhs.y < rhs.y;
    }
    friend bool operator==(const CanvasPoint& lhs, const CanvasPoint& rhs)
    {
      return lhs.y == rhs.y;
    }
    friend bool operator!=(const CanvasPoint& lhs, const CanvasPoint& rhs)
    {
      return !(lhs == rhs);
    }
};

std::ostream& operator<<(std::ostream& os, const CanvasPoint& point)
{
    os << "(" << point.x << ", " << point.y << ", " << point.depth << ") " << point.brightness << std::endl;
    return os;
}
