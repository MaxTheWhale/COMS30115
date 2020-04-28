#pragma once
#include "Times.hpp"

class Temporal {
    public:
        bool realTime;
        float frameRate = 30.0f;
    protected:
        float timeStep();
};