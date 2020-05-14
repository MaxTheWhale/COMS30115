#include "Orbit.hpp"
#include "VectorUtil.hpp"
#include "Animatable.hpp"

using namespace glm;

#ifndef M_PIf
#define M_PIf 3.14159265358979323846f
#endif

Orbit::Orbit(Transformable* center) {
    this->center = center;
}

Orbit::Orbit(Transformable* center, float startPosition) {
    this->center = center;
    this->progress = startPosition;
    this->startPosition = startPosition;
}

bool Orbit::execute(Animatable* parent, Movement& previous) {
    mat4 t = center->transform * Transformable::rotationFromEuler(rotation);
    if (radius < 0) {
        radius = distance(center->getPosition(), parent->getPosition());
    }
    float x = radius * cos(progress);
    float y = radius * sin(progress);
    float scale = parent->timeStep() / time;
    parent->transform[3] = t[3] + normalize(t[0]) * x + normalize(t[1]) * y;
    //parameter increased by timeStep fraction of 2pi scaled by period
    progress += scale * 2 * M_PIf;
    if (progress > 2 * M_PIf + startPosition) {
        progress = startPosition;
        return true;
    }
    return false;
}