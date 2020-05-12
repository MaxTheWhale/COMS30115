#pragma once
#include "Movement.hpp"
#include <glm/glm.hpp>

using namespace glm;

class Orbit : public Movement {
    public:
        Orbit(mat4 transform);
        Orbit(mat4 transform, float startPosition);
        bool execute(Animatable* parent, Movement& previous) override;
    private:
        float startPosition = 0;
        float progress = 0;
        float radius = -1.0f;
};