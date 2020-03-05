#include <glm/glm.hpp>
#include "Colour.h"
#include <string>

class ModelTriangle
{
  public:
    glm::vec4 vertices[3];
    Colour colour;

    ModelTriangle()
    {
    }

    ModelTriangle(glm::vec4 v0, glm::vec4 v1, glm::vec4 v2, Colour trigColour)
    {
      v0.w = 1;
      v1.w = 1;
      v2.w = 1;
      vertices[0] = v0;
      vertices[1] = v1;
      vertices[2] = v2;
      colour = trigColour;
    }

    ModelTriangle(glm::vec3 v0, glm::vec3 v1, glm::vec3 v2, Colour trigColour)
    {
      vertices[0] = glm::vec4(v0.x, v0.y, v0.z, 1);
      vertices[1] = glm::vec4(v1.x, v1.y, v1.z, 1);
      vertices[2] = glm::vec4(v2.x, v2.y, v2.z, 1);
      colour = trigColour;
    }
};

std::ostream& operator<<(std::ostream& os, const ModelTriangle& triangle)
{
    os << "(" << triangle.vertices[0].x << ", " << triangle.vertices[0].y << ", " << triangle.vertices[0].z << ")" << std::endl;
    os << "(" << triangle.vertices[1].x << ", " << triangle.vertices[1].y << ", " << triangle.vertices[1].z << ")" << std::endl;
    os << "(" << triangle.vertices[2].x << ", " << triangle.vertices[2].y << ", " << triangle.vertices[2].z << ")" << std::endl;
    os << std::endl;
    return os;
}
