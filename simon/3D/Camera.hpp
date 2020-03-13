#include <glm/glm.hpp>
#include "Transformable.hpp"

class Camera : public Transformable {
  public:
    glm::vec3 right, up, forward;
    glm::vec3 position;
    glm::mat4 transform, projection, worldToCamera;
    void lookAt(const glm::vec3& from, const glm::vec3& to);
    void setProjection(float fov, float aspect_ratio, float near, float far);
  private:
    glm::mat4 posMat, rotMat;
    void updateTransform();
};