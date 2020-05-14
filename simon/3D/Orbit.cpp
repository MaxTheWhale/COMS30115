#include "Orbit.hpp"
#include "VectorUtil.hpp"
#include "Animatable.hpp"

using namespace glm;

#ifndef M_PIf
#define M_PIf 3.14159265358979323846f
#endif

Orbit::Orbit(mat4 transform) {
    this->transform = transform;
}

Orbit::Orbit(mat4 transform, float startPosition) {
    this->transform = transform;
    this->progress = startPosition;
    this->startPosition = startPosition;
}

bool Orbit::execute(Animatable* parent, Movement& previous) {
    if (radius < 0) {
        radius = distance(getPosition(), parent->getPosition());
    }
    float x = radius * cos(progress);
    float y = radius * sin(progress);
    float scale = parent->timeStep() / time;
    parent->transform[3] = transform[3] + normalize(transform[0]) * x + normalize(transform[1]) * y;
    //parameter increased by timeStep fraction of 2pi scaled by period
    progress += scale * 2 * M_PIf;
    if (progress > 2 * M_PIf + startPosition) {
        progress = 0;
        return true;
    }
    return false;
}