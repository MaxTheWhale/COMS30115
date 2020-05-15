#pragma once

#include <glm/glm.hpp>

//class to wrap a homogeneous 4x4 matrix transform
class Transformable {
  public:
    glm::mat4 transform;
    void move(const glm::vec3& delta);
    void rotate(const glm::vec3& delta);
    void scale(const glm::vec3& delta);
    void setPosition(const glm::vec3& new_position);
    void setRotation(const glm::vec3& new_rotation);
    void setScale(const glm::vec3& new_scale);
    glm::vec3 getScale();
    glm::vec3 getPosition();
    void lookAt(const glm::vec3& from, const glm::vec3& to);
    Transformable();
    static glm::mat4 rotationFromEuler(const glm::vec3& rotation);
};