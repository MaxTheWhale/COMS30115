#include <glm/glm.hpp>
#include "Transformable.hpp"
#include "Animatable.hpp"

class Camera : public Animatable {
  public:
    glm::vec3 right, up, forward;
    glm::vec3 position;
    glm::mat4 projection;//, worldToCamera;
    glm::mat4 worldToCamera();
    float fov = 90.0f;
    float near, far;

    //void lookAt(const glm::vec3& from, const glm::vec3& to); //moved down to Transformable
    void setProjection(float fov, float aspect_ratio, float near, float far);
  private:
    //void updateTransform() override;
};