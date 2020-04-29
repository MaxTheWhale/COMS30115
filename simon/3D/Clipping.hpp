#pragma once

#include "Triangle.hpp"
#include <list>

int clipTriangle(std::list<Triangle>& tris, const glm::vec4& normal);
int clipToView(std::list<Triangle>& tris);