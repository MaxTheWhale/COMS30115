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
