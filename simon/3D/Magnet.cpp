#include "Magnet.hpp"
#include "Rigidbody.hpp"
#include <glm/glm.hpp>

using namespace glm;

Magnet::Magnet(Transformable* center) {
    this->center = center;
}

void Magnet::update() {
    for (unsigned int i = 0; i < Rigidbody::allRBs.size(); i++) {
        Rigidbody* rb = Rigidbody::allRBs[i];
        vec3 dist = center->getPosition() - rb->model->getPosition();
        vec3 force = attractionStrength * normalize(dist) / length(dist);
        rb->applyForce(force, vec3(0,0,0));
    }
}