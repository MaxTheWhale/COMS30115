#include "Rigidbody.hpp"
#include <iostream>
#include "Times.hpp"

using namespace glm;

std::vector<Rigidbody*> Rigidbody::allRBs = std::vector<Rigidbody*>();

Rigidbody::Rigidbody(Model* model) {
    allRBs.push_back(this);
    this->model = model;
    this->velocity = mat4(1,0,0,0,
                               0,1,0,0,
                               0,0,1,0,
                               0,0,0,1);
}

glm::vec3 Rigidbody::gravity = glm::vec3(0, -0.1f, 0);

void Rigidbody::update() {
    //std::cout << "rigidbody update called" << std::endl;
    //std::cout << "allRBs = " << allRBs.size() << std::endl;
    for(unsigned int i = 0; i < allRBs.size(); i++) {
        //std::cout << "this->model = " << this->model << std::endl;
        //std::cout << allRBs[i]->model << std::endl;
        if (allRBs[i]->model != this->model) {
            collide(*allRBs[i]);
        }
    }
    if (hasGravity) {
        vec3 grav = gravity * Times::deltaTime();
        mat4 gravTransform = mat4(1,0,0,0,
                                  0,1,0,0,
                                  0,0,1,0,
                                  grav.x, grav.y, grav.z, 1);
        //std::cout << "gravTransform = " << gravTransform << std::endl;
        velocity *= gravTransform;
    }
    model->transform *= velocity;
}

vec3 toVec3(vec4 in) {
    vec3 vec = vec3(in.x, in.y, in.z);
    return vec;
}

//calculates n
vec3 calcN(const vec3 verts[3]) {
    return cross((verts[1] - verts[0]),
        (verts[2] - verts[0]));
}

//returns false if triangles do not intersect, true if they might
bool checkSign(const float dists[3]) {
    if (dists[0] > 0 && dists[1] > 0 && dists[2] > 0) {
        //positive
        return false;
    }
    else if (dists[0] < 0 && dists[1] < 0 && dists[2] < 0) {
        //negative
        return false;
    }
    else if (dists[0] == 0 && dists[1] == 0 && dists[2] == 0) {
        //coplanar
        std::cout << "coplanar" << std::endl;
    }
    return true;
}

float calcInterval(vec3 D, const vec3 Va, const vec3 Vb, float da, float db) {
    float p0 = dot(D, Va);
    float p1 = dot(D, Vb);
    float t = p0 + (p1 - p0) * (da/(da - db));
    return t;
}

bool intersection(ModelTriangle localTri, ModelTriangle otherTri, mat4 localTransform, mat4 otherTransform) {
    vec3 verts1[3];
    vec3 verts2[3];
    for (int i = 0; i < 3; i++) {
        verts1[i] = toVec3(localTransform * localTri.vertices[i]);
        verts2[i] = toVec3(otherTransform * otherTri.vertices[i]);
    }
    vec3 n1 = calcN(verts1);
    vec3 n2 = calcN(verts2);
    // vec3 n2 = cross((toVec3(otherTransform * otherTri.vertices[1]) - toVec3(otherTransform * otherTri.vertices[0])),
    //     (toVec3(otherTransform * otherTri.vertices[2]) - toVec3(otherTransform * otherTri.vertices[0])));
    float d1 = dot(-n1, verts1[0]);
    float d2 = dot(-n2, verts2[0]);
    float dist1[3];
    float dist2[3];
    for(int i = 0; i < 3; i++) {
        float value1 = dot(n2, verts1[i] + d2);
        float value2 = dot(n1, verts2[i] + d1);
        if (value1 == 0 || value2 == 0) {
            return true;
        }
        dist1[i] = value1;
        dist2[i] = value2;
    }
    if (!checkSign(dist1) || !checkSign(dist2)) {
        return false;
    }
    vec3 D = cross(n1, n2);
    float t1 = calcInterval(D, verts1[0], verts1[1], dist1[0], dist1[1]);
    float t2 = calcInterval(D, verts1[1], verts1[2], dist1[1], dist1[2]);
    float tPrime1 = calcInterval(D, verts2[0], verts2[1], dist2[0], dist2[1]);
    float tPrime2 = calcInterval(D, verts2[1], verts2[2], dist2[1], dist2[2]);
    if (std::max(t1,t2) <= std::min(tPrime1, tPrime2)) { //may not be efficient since standard lib functions may not be inlined by compiler
        return true;
    }
    return false;
}

void Rigidbody::collide(Rigidbody other) {
    if (!collisionEnabled || !other.collisionEnabled) {
        return;
    }
    for (unsigned int i = 0; i < model->tris.size(); i++) {
        for (unsigned int j = 0; j < other.model->tris.size(); j++) {
            if (intersection(model->tris[i], other.model->tris[j], model->transform, other.model->transform)) {
                //std::cout << "collision detected" << std::endl;
                hasGravity = false;
                velocity = mat4();
            }
            // else {
            //     std::cout << "no collision detected" << std::endl;
            // }
        }
    }
}