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
        Movement(glm::vec3 rotation, float time);
        float time;
        glm::vec3 rotation;
        glm::vec3 prevRotation;
        int repeats = 1;
        virtual bool execute(Animatable* parent, Movement& previous);
        bool stareAt = false;
        bool isRotation = false;
        Transformable* stareTarget;
    private:
        void finish(Animatable* parent);
        float elapsed = 0;
};