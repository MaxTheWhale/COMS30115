#include "Transformable.hpp"
#include <iostream>

Transformable::Transformable() {
    std::cout << "Transformable constructor called" << std::endl;
    this->transform= glm::mat4(1, 0, 0, 0,
                            0, 1, 0, 0,
                            0, 0, 1, 0,
                            0, 0, 0, 1);
}

void Transformable::move(const glm::vec3& delta) {
    transform[3].x += delta.x;
    transform[3].y += delta.y;
    transform[3].z += delta.z;
}
void Transformable::rotate(const glm::vec3& delta) {
    transform *= rotationFromEuler(delta);
}
void Transformable::scale(const glm::vec3& delta) {
    transform[0] *= delta.x;
    transform[1] *= delta.y;
    transform[2] *= delta.z;
}
void Transformable::setPosition(const glm::vec3& new_position) {
    transform[3] = glm::vec4(new_position.x, new_position.y, new_position.z, 1);
}
void Transformable::setRotation(const glm::vec3& new_rotation) {
    glm::mat4 rotationMat = rotationFromEuler(new_rotation);
    for (int i = 0; i < 3; i++) {
        float scale = glm::length(transform[i]);
        transform[i] = glm::normalize(rotationMat[i]) * scale;
    }
}
void Transformable::setScale(const glm::vec3& new_scale) {
    for (int i = 0; i < 3; i++) {
        transform[i] = glm::normalize(transform[i]) * new_scale[i];
    }
}
// void Transformable::updateTransform() {
//     transform = positionMat * rotationMat * scaleMat;
// }

glm::mat4 Transformable::rotationFromEuler(const glm::vec3& rotation) {
    float sa, sb, sc, ca, cb, cc;
    sa = sinf(rotation.x);
    sb = sinf(rotation.y);
    sc = sinf(rotation.z);
    ca = cosf(rotation.x);
    cb = cosf(rotation.y);
    cc = cosf(rotation.z);
    return glm::mat4(cb * cc, sc, cc * -sb, 0,
                     sa * sb - ca * cb * sc, cc * ca, sb * sc * ca + cb * sa, 0,
                     cb * sc * sa + sb * ca, sa * -cc, cb * ca - sb * sc * sa, 0,
                     0, 0, 0, 1);
}