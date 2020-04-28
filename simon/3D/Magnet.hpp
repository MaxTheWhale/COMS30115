#pragma once

#include "Transformable.hpp"
#include "Updatable.hpp"

class Magnet : public Updatable {
    public:
        Magnet(Transformable* center);
        float attractionStrength = 0.25f;
        void update() override;
    protected:
        Transformable* center;
};