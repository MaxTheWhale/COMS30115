#include "Magnet.hpp"
#include "Rigidbody.hpp"
#include <glm/glm.hpp>

using namespace glm;

vec3 Magnet::totalForce = vec3(0,0,0);

Magnet::Magnet(Transformable* center) {
    this->center = center;
}

void Magnet::update() {
    for (unsigned int i = 0; i < Rigidbody::allRBs.size(); i++) {
        // Rigidbody* rb = Rigidbody::allRBs[i];
        if (!Rigidbody::allRBs[i]->suckable || Rigidbody::allRBs[i]->model->transform == this->center->transform) { 
            continue;
        }
        vec3 dist = center->getPosition() - Rigidbody::allRBs[i]->model->getPosition();
        vec3 force = timeStep() * attractionStrength * dist / (length(dist) * length(dist));
        // cout << "force = " << force << endl;
        totalForce += force;
        Rigidbody::allRBs[i]->applyForce(force, vec3(0,0,0));
    }
}