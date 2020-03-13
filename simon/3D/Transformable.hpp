#include <glm/glm.hpp>

class Transformable {
  public:
    glm::mat4 transform;
    void move(const glm::vec3& delta);
    void rotate(const glm::vec3& delta);
    void scale(const glm::vec3& delta);
    void setPosition(const glm::vec3& new_position);
    void setRotation(const glm::vec3& new_rotation);
    void setScale(const glm::vec3& new_scale);
  protected:
    glm::mat4 positionMat, rotationMat, scaleMat;
    void updateTransform();
};