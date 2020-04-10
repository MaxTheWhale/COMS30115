#include "VectorOutput.hpp"

using namespace std;
using namespace glm;

std::ostream &operator<<(std::ostream &os, const std::vector<float> vector)
{
  os << '[';
  for (int i = 0; i < (int)vector.size(); i++)
  {
    os << vector[i];
    if (i != (int)vector.size() - 1)
    {
      os << ',';
    }
  }
  os << ']';
  return os;
}
std::ostream &operator<<(std::ostream &os, const vec3 vec3)
{
  std::vector<float> vec{vec3.x, vec3.y, vec3.z};
  os << vec;
  return os;
}
std::ostream &operator<<(std::ostream &os, const vec4 vec4)
{
  std::vector<float> vec{vec4.x, vec4.y, vec4.z, vec4.w};
  os << vec;
  return os;
}
std::ostream &operator<<(std::ostream &os, const mat3 mat3)
{
  for (int i = 0; i < 3; i++)
  {
    os << mat3[i];
  }
  return os;
}
std::ostream &operator<<(std::ostream &os, const mat4 mat4)
{
  for (int i = 0; i < 4; i++)
  {
    os << mat4[i];
  }
  return os;
}