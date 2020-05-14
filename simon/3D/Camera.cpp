#include "Camera.hpp"

void Camera::setProjection(float fov, float aspect_ratio, float near, float far) {
  float top = tan(glm::radians(fov / 2)) * near;
  float right = top * aspect_ratio;
  projection = glm::mat4(near / right, 0, 0, 0,
                         0, near / top, 0, 0,
                         0, 0, -(far + near) / (far - near), -1,
                         0, 0, -(2 * far * near) / (far - near), 0);
  this->near = near;
  this->far = far;
  this->fov = fov;
}

glm::mat4 Camera::worldToCamera() {
  glm::vec4 pos = transform[3];
  glm::mat4 rotationMat = transform;
  rotationMat[3] = glm::vec4(0, 0, 0, 1);
  return glm::transpose(rotationMat) *
                  glm::mat4(1, 0, 0, 0,
                            0, 1, 0, 0,
                            0, 0, 1, 0,
                            -pos.x, -pos.y, -pos.z, 1);
}