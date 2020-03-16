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
    positionMat[3].x += delta.x;
    positionMat[3].y += delta.y;
    positionMat[3].z += delta.z;
    updateTransform();
}
void Transformable::rotate(const glm::vec3& delta) {
    rotationMat *= rotationFromEuler(delta);
    updateTransform();
}
void Transformable::scale(const glm::vec3& delta) {
    scaleMat[0].x *= delta.x;
    scaleMat[1].y *= delta.y;
    scaleMat[2].z *= delta.z;
    updateTransform();
}
void Transformable::setPosition(const glm::vec3& new_position) {
    positionMat = glm::mat4(1, 0, 0, 0,
                            0, 1, 0, 0,
                            0, 0, 1, 0,
                            new_position.x, new_position.y, new_position.z, 1);
    updateTransform();
}
void Transformable::setRotation(const glm::vec3& new_rotation) {
    rotationMat = rotationFromEuler(new_rotation);
    updateTransform();
}
void Transformable::setScale(const glm::vec3& new_scale) {
    scaleMat = glm::mat4(new_scale.x, 0, 0, 0,
                         0, new_scale.y, 0, 0,
                         0, 0, new_scale.z, 0,
                         0, 0, 0, 1);
    updateTransform();
}
void Transformable::updateTransform() {
    transform = positionMat * rotationMat * scaleMat;
}

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