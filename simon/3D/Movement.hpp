#pragma once

#include "Transformable.hpp"
#include <glm/glm.hpp>

class Animatable; //I hate this damn language

class Movement : public Transformable {
    public:
        Movement();
        Movement(glm::mat4 transform);
        Movement(float time);
        Movement(glm::mat4 transform, float time);
        float time;
        int repeats = 1;
        virtual bool execute(Animatable* parent, Movement& previous);
        bool stareAt = false;
        glm::vec3 stareTarget;
};