#pragma once

#include "Transformable.hpp"
#include "Animatable.hpp"

class Light : public Animatable {
  public:
    glm::vec3 diffuseIntensity, specularIntensity;
    Light(const glm::vec3& diffuseIntensity, const glm::vec3& specularIntensity);
};