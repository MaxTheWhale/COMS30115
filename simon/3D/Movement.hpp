#pragma once

#include "Transformable.hpp"

class Movement : public Transformable {
    public:
        Movement();
        Movement(float time);
        float time;
};