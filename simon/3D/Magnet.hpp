#pragma once

#include "Transformable.hpp"
#include "Updatable.hpp"
#include "Temporal.hpp"
#include "Rigidbody.hpp"
#include <glm/glm.hpp>

class Magnet : public Updatable, Temporal {
    public:
        Magnet(Transformable* center, vector<Rigidbody*>& rbList);
        float attractionStrength = 10.0f;
        vector<Rigidbody*> *allRBs;
        void update() override;
        static glm::vec3 totalForce;
    protected:
        Transformable* center;
};