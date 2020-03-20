#include "Movement.hpp"

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