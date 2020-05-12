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

bool Movement::execute(Animatable* parent, Movement& previous) {
    // cout << "delta = deltaTime (" << Times::deltaTime() << ") * (previous transform (" << this->previous.transform << ") - next transform (" << moves.top().transform << "))" << endl;
    float scale = -parent->timeStep() / time;
    cout << "scale = " << scale << '\n';

    if (isRotation) {
        vec3 rotDelta = -scale * rotation;
        if (((length(prevRotation + rotDelta) >= length(rotation)) && (length(prevRotation) < length(rotation))) ||
            ((length(prevRotation + rotDelta) <= length(rotation)) && (length(prevRotation) > length(rotation)))) {
                return true;
            }
        else {
            prevRotation += rotDelta;
            cout << "prevRot: " << prevRotation << '\n';
            parent->transform = rotationFromEuler(rotDelta) * parent->transform;
            return false;
        }
    }
    else {
        mat4 delta = scale * (previous.transform - transform);
        if (stareAt) {
            for (int i = 0; i < 3; i++) {
                delta[i] = vec4(0,0,0,0);
            }
        }
        // cout << "previous = " << previous.transform << endl << "target = " << moves.top().transform << endl << "delta = " << delta << endl;
        mat4 newTransform = mat4();
        for (int i = 0; i < 4; i++)
        {
            newTransform[i] = parent->transform[i] + delta[i];
        }
        if (stareAt) {
            cout << "stare target: " << stareTarget << endl;
            parent->lookAt(parent->getPosition(), stareTarget);
        }
        //would the move take us further from our goal
        float currentDist;
        float newDist;
        if (stareAt) {
            currentDist = distance(parent->transform[3], transform[3]);
            newDist = distance(newTransform[3], transform[3]);
        }
        else {
            currentDist = mat4Dist(parent->transform, transform);
            newDist = mat4Dist(newTransform, transform);
        }
        // cout << "currentDist = " << currentDist << " newDist = " << newDist << endl;
        if (currentDist < newDist) {
            // if (stareAt) {
            //     parent->transform[3] = transform[3];
            // }
            // else {
            //     parent->transform = transform;
            // }
            return true;
        }
        else {
            for (int i = 0; i < 4; i++) {
                parent->transform[i] += delta[i];
            }
        }
        return false;
    }

}