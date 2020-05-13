#pragma once

#include <iostream>
#include <Colour.h>
#include <vector>
#include <ModelTriangle.h>
#include <glm/glm.hpp>

class LightOld
{
  public:
    std::string name;
    Colour colour = Colour();
    std::vector<ModelTriangle> triangles = std::vector<ModelTriangle>();
    vec4 centre = vec4();
    float intensity = 10.0f;
    float shadow = 0.1f;

    LightOld()
    {
    }

    LightOld(std::string n, vector<ModelTriangle> tris)
    {
        name = n;
        triangles = tris;
    }

    void calculateCentre() {
        vector<vec4> vertices = vector<vec4>();
        vec4 average = vec4(0, 0, 0, 0);
        for(ModelTriangle light : triangles) {
            vertices.push_back(light.vertices[0]);
            vertices.push_back(light.vertices[1]);
            vertices.push_back(light.vertices[2]);

            average += light.vertices[0];
            average += light.vertices[1];
            average += light.vertices[2];
        }

        average.x /= vertices.size();
        average.y /= vertices.size();
        average.z /= vertices.size();
        average.w = 1;

        centre = average;
    }

};

std::ostream& operator<<(std::ostream& os, const Colour& colour);