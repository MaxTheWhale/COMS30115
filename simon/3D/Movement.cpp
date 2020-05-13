#include "Movement.hpp"
#include "Animatable.hpp"

Movement::Movement() {
    time = 1.0f;
    transform = glm::mat4(1,0,0,0,
                                0,1,0,0,
                                0,0,1,0,
                                0,0,0,1);
}
Movement::Movement(glm::mat4 transform) {
    time = 1.0f;
    this->transform = transform;
}
Movement::Movement(float time) {
    this->time = time;
}
Movement::Movement(glm::mat4 transform, float time) {
    this->time = time;
    this->transform = transform;
}
Movement::Movement(glm::vec3 rotation, float time) {
    this->time = time;
    this->rotation = rotation;
    this->prevRotation = glm::vec3(0, 0, 0);
    this->isRotation = true;
}

float mat4Dist(mat4 &a, mat4 &b)
{
    float result = 0;
    for (int i = 0; i < 4; i++)
    {
        result += distance(a[i], b[i]);
    }
    return result;
}

void Movement::finish(Animatable* parent) {
    parent->transform = transform;
    parent->transform = rotationFromEuler(rotation) * parent->transform;
}

bool Movement::execute(Animatable* parent, Movement& previous) {
    // cout << "delta = deltaTime (" << Times::deltaTime() << ") * (previous transform (" << this->previous.transform << ") - next transform (" << moves.top().transform << "))" << endl;
    if (time <= 0) {
        finish(parent);
        return true;
    }
    elapsed += parent->timeStep();
    if (elapsed >= time) {
        // finish(parent);
        return true;
    }
    float scale = -parent->timeStep() / time;
    // cout << "scale = " << scale << '\n';

    if (isRotation) {
        vec3 rotDelta = -scale * rotation;
        if (((length(prevRotation + rotDelta) >= length(rotation)) && (length(prevRotation) < length(rotation))) ||
            ((length(prevRotation + rotDelta) <= length(rotation)) && (length(prevRotation) > length(rotation)))) {
                // finish(parent);
                // return true;
            }
        else {
            prevRotation += rotDelta;
            cout << "prevRot: " << prevRotation << '\n';
            cout << "parent transform: " << parent->transform << endl;
            parent->transform = parent->transform * rotationFromEuler(rotDelta);
            //return false;
        }
    }
    vec4 delta = scale * (previous.transform[3] - transform[3]);
    // cout << "previous = " << previous.transform << endl << "target = " << moves.top().transform << endl << "delta = " << delta << endl;
    mat4 newTransform = parent->transform;
    newTransform[3] += delta;
    if (stareAt) {
        // cout << "stare target: " << stareTarget << endl;
        parent->lookAt(parent->getPosition(), stareTarget->getPosition());
    }
    //would the move take us further from our goal
    float currentDist = distance(parent->transform[3], transform[3]);
    float newDist = distance(newTransform[3], transform[3]);
    // cout << "currentDist = " << currentDist << " newDist = " << newDist << endl;
    if (currentDist < newDist) {
        // if (stareAt) {
        //     parent->transform[3] = transform[3];
        // }
        // else {
        //     parent->transform = transform;
        // }
        // finish(parent);
        // return true;
    }
    else {
        parent->transform[3] += delta;
    }
    return false;

}