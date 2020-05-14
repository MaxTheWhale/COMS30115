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
    if (time == 0) {
        finish(parent);
        return true;
    }
    elapsed += parent->timeStep();
    if (elapsed >= time && time > 0) {
        cout << "movement completed" << endl;
        return true;
    }
    float scale = -parent->timeStep() / time;

    if (isRotation) {
        vec3 rotDelta = -scale * rotation;
        if (((length(prevRotation + rotDelta) >= length(rotation)) && (length(prevRotation) < length(rotation))) ||
            ((length(prevRotation + rotDelta) <= length(rotation)) && (length(prevRotation) > length(rotation)))) {
            }
        else {
            prevRotation += rotDelta;
            parent->transform = parent->transform * rotationFromEuler(rotDelta);
        }
    }
    vec4 delta = scale * (previous.transform[3] - transform[3]);
    mat4 newTransform = parent->transform;
    newTransform[3] += delta;
    if (stareAt) {
        parent->lookAt(parent->getPosition(), stareTarget->getPosition());
    }
    //would the move take us further from our goal
    float currentDist = distance(parent->transform[3], transform[3]);
    float newDist = distance(newTransform[3], transform[3]);
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