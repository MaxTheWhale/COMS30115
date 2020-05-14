#pragma once

#include "Times.hpp"

class Temporal {
    public:
        bool realTime = false;
        float frameRate = 30.0f;
        float timeStep();
};