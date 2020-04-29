#pragma once

#include <glm/glm.hpp>

class Vertex {
  public:
    glm::vec4 pos;
    glm::vec4 pos_3d;
    glm::vec4 normal;
    float brightness;
    float u, v;

    Vertex operator+=(const Vertex& rhs);
    friend Vertex operator+(Vertex lhs, const Vertex& rhs);
    Vertex operator-=(const Vertex& rhs);
    friend Vertex operator-(Vertex lhs, const Vertex& rhs);
    Vertex operator*=(float rhs);
    friend Vertex operator*(Vertex lhs, float rhs);
    Vertex operator/=(float rhs);
    friend Vertex operator/(Vertex lhs, float rhs);
};