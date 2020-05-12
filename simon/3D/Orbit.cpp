#include "Orbit.hpp"
#include "VectorUtil.hpp"
#include "Animatable.hpp"

using namespace glm;

# define pi 3.14159265358979323846

Orbit::Orbit(mat4 transform) {
    this->transform = transform;
}

Orbit::Orbit(mat4 transform, float startPosition) {
    this->transform = transform;
    this->progress = startPosition;
    this->startPosition = startPosition;
}

bool Orbit::execute(Animatable* parent, Movement& previous) {
    if (radius < 0) {
        radius = distance(getPosition(), parent->getPosition());
        cout << "radius set to " << radius << endl;
    }
    float x = radius * cos(progress);
    float y = radius * sin(progress);
    cout << "x,y = " << x << "," << y << endl << transform << endl;
    // parent->setPosition(getPosition() + vec3(x,y,0));
    float scale = parent->timeStep() / time;
    parent->transform[3] = transform[3] + normalize(transform[0]) * x + normalize(transform[1]) * y;
    //parameter increased by timeStep fraction of 2pi scaled by period
    progress += scale * 2 * pi;
    if (progress > 2 * pi + startPosition) {
        progress = 0;
        // radius = -1;
        // this->setPosition(parent->getPosition()); //I know what I'm about son
        return true;
    }
    return false;
}