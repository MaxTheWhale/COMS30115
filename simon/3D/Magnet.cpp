#include "Magnet.hpp"
#include "Rigidbody.hpp"
#include <glm/glm.hpp>

using namespace glm;

vec3 Magnet::totalForce = vec3(0,0,0);

Magnet::Magnet(Transformable* center, vector<Rigidbody*>& rbList) {
    this->allRBs = &rbList;
    this->center = center;
}

void Magnet::update() {
    for (unsigned int i = 0; i < (*allRBs).size(); i++) {
        Rigidbody* rb = (*allRBs)[i];
        if (!rb->suckable || rb->model->transform == this->center->transform) { 
            continue;
        }
        vec3 dist = center->getPosition() - rb->model->getPosition();
        vec3 force = timeStep() * attractionStrength * dist / (length(dist) * length(dist));
        totalForce += force;
        rb->applyForce(force, vec3(0,0,0));
    }
}