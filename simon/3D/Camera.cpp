#include "Camera.hpp"

void Camera::lookAt(const glm::vec3& from, const glm::vec3& to) {
  forward = glm::normalize(from - to);
  right = glm::normalize(glm::cross(glm::vec3(0, 1, 0), forward));
  up = glm::normalize(glm::cross(forward, right));
  setPosition(from);
  rotationMat = glm::mat4(right.x, right.y, right.z, 0,
                     up.x, up.y, up.z, 0,
                     forward.x, forward.y, forward.z, 0,
                     0, 0, 0, 1);
  updateTransform();
}

void Camera::setProjection(float fov, float aspect_ratio, float near, float far) {
  float top = tan(glm::radians(fov / 2)) * near;
  float bottom = -top;
  float right = top * aspect_ratio;
  float left = bottom;
  projection = glm::mat4(2 * near / (right - left), 0, 0, 0,
                         0, 2 * near / (top - bottom), 0, 0,
                         (right + left) / (right - left), (top + bottom) / (top - bottom), -(far + near) / (far - near), -1,
                         0, 0, -(2 * far * near) / (far - near), 0);
}

void Camera::updateTransform() {
  transform = rotationMat * positionMat;
  glm::vec4 pos = positionMat[3];
  worldToCamera = glm::transpose(rotationMat) *
                  glm::mat4(1, 0, 0, 0,
                            0, 1, 0, 0,
                            0, 0, 1, 0,
                            -pos.x, -pos.y, -pos.z, 1);
}