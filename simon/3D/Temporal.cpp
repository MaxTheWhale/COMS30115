#include "Temporal.hpp"

float Temporal::timeStep() {
    if (realTime) {
        return Times::deltaTime();
    }
    else {
        return 1.0f/frameRate;
    }
}