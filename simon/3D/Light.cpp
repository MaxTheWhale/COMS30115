#include "Light.hpp"

using namespace std;

Light::Light(const glm::vec3& diffuseIntensity, const glm::vec3& specularIntensity) {
    this->diffuseIntensity = diffuseIntensity;
    this->specularIntensity = specularIntensity;
}