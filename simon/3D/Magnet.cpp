#include "Magnet.hpp"
#include "Rigidbody.hpp"
#include <glm/glm.hpp>

using namespace glm;

vec3 Magnet::totalForce = vec3(0,0,0);

Magnet::Magnet(Transformable* center) {
    this->center = center;
}

void Magnet::update() {
    for (unsigned int i = 0; i < Rigidbody::getAllRBs().size(); i++) {
        Rigidbody* rb = Rigidbody::getAllRBs()[i];
        if (!rb->suckable || rb->model->transform == this->center->transform) { 
            continue;
        }
        vec3 dist = center->getPosition() - rb->model->getPosition();
        vec3 force = timeStep() * attractionStrength * dist / (length(dist) * length(dist));
        // cout << "force = " << force << endl;
        totalForce += force;
        rb->applyForce(force, vec3(0,0,0));
    }
}