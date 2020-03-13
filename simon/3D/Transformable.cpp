#include "Transformable.hpp"

void Transformable::move(const glm::vec3& delta) {
    positionMat[3].x += delta.x;
    positionMat[3].y += delta.y;
    positionMat[3].z += delta.z;
    updateTransform();
}
void Transformable::rotate(const glm::vec3& delta) {
}
void Transformable::scale(const glm::vec3& delta) {
    
}
void Transformable::setPosition(const glm::vec3& new_position) {
    
}
void Transformable::setRotation(const glm::vec3& new_rotation) {
    float sa, sb, sc, ca, cb, cc;
    sa = sinf(new_rotation.x);
    sb = sinf(new_rotation.y);
    sc = sinf(new_rotation.z);
    ca = cosf(new_rotation.x);
    cb = cosf(new_rotation.y);
    cc = cosf(new_rotation.z);
    rotationMat = glm::mat4(cb * cc, sc, cc * -sb, 0,
                            sa * sb - ca * cb * sc, cc * ca, sb * sc * ca + cb * sa, 0,
                            cb * sc * sa + sb * ca, sa * -cc, cb * ca - sb * sc * sa, 0,
                            0, 0, 0, 1);
    updateTransform();
}
void Transformable::setScale(const glm::vec3& new_scale) {
    
}
void Transformable::updateTransform() {
    transform = positionMat * rotationMat * scaleMat;
}