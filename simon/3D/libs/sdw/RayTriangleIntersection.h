#pragma once

#include <glm/glm.hpp>
#include <iostream>
#include <ModelTriangle.h>
#include "../../VectorUtil.hpp"

class RayTriangleIntersection
{
  public:
    glm::vec4 intersectionPoint;
    float distanceFromCamera;
    ModelTriangle intersectedTriangle;
    float u;
    float v;

    RayTriangleIntersection()
    {
    }

    RayTriangleIntersection(glm::vec4 point, float distance, ModelTriangle triangle)
    {
        intersectionPoint = point;
        distanceFromCamera = distance;
        intersectedTriangle = triangle;
    }
};

std::ostream& operator<<(std::ostream& os, const RayTriangleIntersection& intersection)
{
    os << "Intersection is at " << intersection.intersectionPoint << " on triangle " << intersection.intersectedTriangle << " at a distance of " << intersection.distanceFromCamera << std::endl;
    return os;
}
