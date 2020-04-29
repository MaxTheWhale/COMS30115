#pragma once

#include <iostream>
#include <vector>
#include <glm/glm.hpp>

using namespace std;
using namespace glm;

std::ostream &operator<<(std::ostream &os, const std::vector<float> vector);

std::ostream &operator<<(std::ostream &os, const vec3 vec3);

std::ostream &operator<<(std::ostream &os, const vec4 vec4);

std::ostream &operator<<(std::ostream &os, const mat3 mat3);

std::ostream &operator<<(std::ostream &os, const mat4 mat4);

inline vec3 toThree(vec4 v) {
  return vec3(v.x, v.y, v.z);
}

inline vec4 cross(const vec4& a, const vec4& b) {
  return vec4(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x, 0);
}
