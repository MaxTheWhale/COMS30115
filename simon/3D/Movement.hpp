#pragma once

#include "Transformable.hpp"
#include <glm/glm.hpp>

class Movement : public Transformable {
    public:
        Movement();
        Movement(glm::mat4 transform);
        float time;
};