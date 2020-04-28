#pragma once
#include "Times.hpp"

class Temporal {
    public:
        bool realTime = false;
        float frameRate = 60.0f;
    protected:
        float timeStep();
};