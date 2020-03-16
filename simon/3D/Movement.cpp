#include "Movement.hpp"

Movement::Movement() {
    this->time = 0.0f;
}
Movement::Movement(float time) {
    this->time = time;
    this->transform = glm::mat4(1,0,0,0,
                                0,1,0,0,
                                0,0,1,0,
                                0,0,0,1);
}