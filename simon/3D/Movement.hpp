#pragma once

#include "Transformable.hpp"
#include <glm/glm.hpp>

class Animatable; // I hate this damn language

//class to represent a single animated movement
class Movement : public Transformable {
    public:
        Movement();
        Movement(glm::mat4 transform);
        Movement(float time);
        Movement(glm::mat4 transform, float time);
        Movement(glm::vec3 rotation, float time);
        float time;
        bool isRotation = false;
        glm::vec3 rotation; //rotation is handled as a seperate vec3 to improve performance and simplify the maths
        glm::vec3 prevRotation;
        int repeats = 1;
        virtual bool execute(Animatable* parent, Movement& previous);
        bool stareAt = false; //sets the model's rotation to point at a specified transformable every frame
        Transformable* stareTarget;
    private:
        void finish(Animatable* parent);
        float elapsed = 0;
};