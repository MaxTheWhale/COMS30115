#include "Camera.hpp"

// void Camera::lookAt(const glm::vec3& from, const glm::vec3& to) {
//   vec3 forward = glm::normalize(from - to);
//   right = glm::normalize(glm::cross(glm::vec3(0, 1, 0), forward));
//   up = glm::normalize(glm::cross(forward, right));
//   //setPosition(from);
//   transform = glm::mat4(right.x, right.y, right.z, 0,
//                      up.x, up.y, up.z, 0,
//                      forward.x, forward.y, forward.z, 0,
//                      from.x, from.y, from.z, 1);
// }

void Camera::setProjection(float fov, float aspect_ratio, float near, float far) {
  float top = tan(glm::radians(fov / 2)) * near;
  float right = top * aspect_ratio;
  projection = glm::mat4(near / right, 0, 0, 0,
                         0, near / top, 0, 0,
                         0, 0, -(far + near) / (far - near), -1,
                         0, 0, -(2 * far * near) / (far - near), 0);
  (*this).near = near;
  (*this).far = far;
  (*this).fov = fov;
}

// void Camera::updateTransform() {
//   //transform = rotationMat * positionMat;
//   glm::vec4 pos = transform[3];
//   glm::mat4 rotationMat = transform;
//   rotationMat[3] = glm::vec4(0,0,0,0);
//   worldToCamera = glm::transpose(rotationMat) *
//                   glm::mat4(1, 0, 0, 0,
//                             0, 1, 0, 0,
//                             0, 0, 1, 0,
//                             -pos.x, -pos.y, -pos.z, 1);
// }

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