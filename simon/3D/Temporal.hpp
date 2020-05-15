#pragma once

#include "Times.hpp"

//class for synching all moving parts to the same framerate
//also allows for running in adjusted real-time like a game engine (ie Unity's Time.deltaTime)
class Temporal {
    public:
        bool realTime = false;
        float frameRate = 60.0f;
        float timeStep();
};