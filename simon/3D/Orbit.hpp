#pragma once

#include "Movement.hpp"
#include <glm/glm.hpp>

using namespace glm;

class Orbit : public Movement {
    public:
        Orbit(Transformable* center);
        Orbit(Transformable* center, float startPosition); //start position is between 0 and 2π
        bool execute(Animatable* parent, Movement& previous) override;
    private:
        float startPosition = 0;
        float progress = 0;
        float radius = -1.0f;
        Transformable* center;
};