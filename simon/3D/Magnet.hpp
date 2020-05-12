#pragma once

#include "Transformable.hpp"
#include "Updatable.hpp"
#include "Temporal.hpp"
#include <glm/glm.hpp>

class Magnet : public Updatable, Temporal {
    public:
        Magnet(Transformable* center);
        float attractionStrength = 20.0f;
        void update() override;
        static glm::vec3 totalForce;
    protected:
        Transformable* center;
};