#include "Rigidbody.hpp"
#include <iostream>
#include "Times.hpp"
#include <algorithm>

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

glm::vec3 Rigidbody::gravity = glm::vec3(0, -0.5f, 0);


void Rigidbody::update() {
    //if we're not moving and we're in contact with something then it can be assumed that we are resting on it
    if (hasGravity && !(velocity[3] == vec4(0,0,0,1) && !collidedWith.empty())) {
        float timescale = realTimeScale ? Times::deltaTime() : 1.0f/30.0f;
        vec3 grav = gravity * timescale;
        mat4 gravTransform = mat4(1,0,0,0,
                                  0,1,0,0,
                                  0,0,1,0,
                                  grav.x, grav.y, grav.z, 1);
        //std::cout << "gravTransform = " << gravTransform << std::endl;
        velocity *= gravTransform;
    }
    //clear the collision list before running collision checks again
    collidedWith.clear();

    mat4 oldTransform = model->transform;

    model->transform *= velocity;
    //disgusting hack to make linear velocity independent of angular
    model->transform[3] = oldTransform[3] + velocity[3];
    model->transform[3][3] = 1;

    for(unsigned int i = 0; i < allRBs.size(); i++) {
        if (allRBs[i]->model != this->model) {
            if(collide(*allRBs[i])) {
                // hasGravity = false;
                // velocity = mat4();
                model->transform = oldTransform;
            }
        }
    }
    cout << "velocity = " << velocity << endl;
    cout << "transform = " << model->transform << endl;
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
    if (dists[0] != 0 && dists[1] != 0 && dists[2] != 0) { //no point on plane
        if (dists[0] > 0 && dists[1] > 0 && dists[2] > 0) { //all on positive side
            //positive
            return false;
        }
        else if (dists[0] < 0 && dists[1] < 0 && dists[2] < 0) { // all on negative side
            //negative
            return false;
        }
    }
    else if (dists[0] == 0 && dists[1] == 0 && dists[2] == 0) { //all points on plane
        //coplanar
        std::cout << "coplanar" << std::endl;
    }
    return true;
}

float calcInterval(vec3 D, const vec3 Va, const vec3 Vb, float da, float db) {
    float p0 = dot(D, Va);
    float p1 = dot(D, Vb);
    //cout << "ratio = " << (da/(da - db)) << endl;
    float t = p0 + (p1 - p0) * (da/(da - db));
    return t;
}

bool intersection(ModelTriangle localTri, ModelTriangle otherTri, mat4 localTransform, mat4 otherTransform) {
    vec3 verts1[3];
    vec3 verts2[3];
    //cout << "verts1 = ";
    for (int i = 0; i < 3; i++) {
        verts1[i] = toVec3(localTransform * localTri.vertices[i]);
        //cout << verts1[i] << ", ";
        verts2[i] = toVec3(otherTransform * otherTri.vertices[i]);
    }
    //cout << endl << "verts2 = ";
    // for (int i = 0; i < 3; i++) {
    //     cout << verts2[i] << ", ";
    // }
    // cout << endl;
    vec3 n1 = calcN(verts1);
    vec3 n2 = calcN(verts2);

    float d1 = dot(-n1, verts1[0]);
    float d2 = dot(-n2, verts2[0]);

    //distances from each vertex to the plane of the other triangle
    float dist1[3];
    float dist2[3];
    for(int i = 0; i < 3; i++) {
        float value1 = dot(n2, verts1[i]) + d2;
        float value2 = dot(n1, verts2[i]) + d1;
        dist1[i] = value1;
        dist2[i] = value2;
    }
    //early check to weed out cases where one triangle is distinctly on one side of the other triangle's plane
    if (!checkSign(dist1) || !checkSign(dist2)) {
        //std::cout << "early return" << endl;
        return false;
    }
    vec3 D = cross(n1, n2);
    //cout << "D = "<< D << endl;
    float t1 = calcInterval(D, verts1[0], verts1[1], dist1[0], dist1[1]);
    float t2 = calcInterval(D, verts1[1], verts1[2], dist1[1], dist1[2]);
    float tPrime1 = calcInterval(D, verts2[0], verts2[1], dist2[0], dist2[1]);
    float tPrime2 = calcInterval(D, verts2[1], verts2[2], dist2[1], dist2[2]);
    //std::cout << "interval1: " << t1 << "," << t2 << " interval2: " << tPrime1 << "," << tPrime2 << endl;
    if (std::max(t1,t2) <= std::min(tPrime1, tPrime2) || std::min(t1,t2) <= std::max(tPrime1, tPrime2)) { //may not be efficient since standard lib functions may not be inlined by compiler
        return true;
    }
    //std::cout << "default" << endl;
    return false;
}

bool Rigidbody::collide(Rigidbody other) {
    // return false;
    if (positionFixed || !collisionEnabled || !other.collisionEnabled) {
        return false;
    }
    for (unsigned int i = 0; i < model->tris.size(); i++) {
        for (unsigned int j = 0; j < other.model->tris.size(); j++) {
            //this line is disgusting because C++'s iterator interface hurts my soul
            if (std::find_if(collidedWith.begin(), collidedWith.end(), [other](Rigidbody* rb) {return &other == rb;}) == collidedWith.end() &&
                    intersection(model->tris[i], other.model->tris[j], model->transform, other.model->transform)) {
                if (velocity[3] != vec4(0,0,0,1)) { //object is resting in contact with other object, so no reaction force is applied
                    std::cout << "collision detected with " << &other << std::endl;
                    vec4 normal = other.model->tris[j].normal * other.model->transform;
                    float combinedElasticity = this->elasticity + other.elasticity; //twice the average ratio of energy conserved
                    vec3 force = toVec3(-combinedElasticity * dot(normal, velocity[3]) * normal);
                    cout << "force = " << force << endl;
                    applyForce(force, vec3(0, 0, 0));
                    for (int i = 0; i < 3; i++) {
                        if (velocity[3][i] < 0.01f) {
                            velocity[3][i] = 0;
                        }
                    }
                }
                collidedWith.push_back(&other);
                other.collidedWith.push_back(this);
                return true;
            }
        }
    }
    return false;
}

glm::mat4 rotationFromEuler(const glm::vec3& rotation) {
    float sa, sb, sc, ca, cb, cc;
    sa = sinf(rotation.x);
    sb = sinf(rotation.y);
    sc = sinf(rotation.z);
    ca = cosf(rotation.x);
    cb = cosf(rotation.y);
    cc = cosf(rotation.z);
    return glm::mat4(cb * cc, sc, cc * -sb, 0,
                     sa * sb - ca * cb * sc, cc * ca, sb * sc * ca + cb * sa, 0,
                     cb * sc * sa + sb * ca, sa * -cc, cb * ca - sb * sc * sa, 0,
                     0, 0, 0, 1);
}

// vec3 normalizeSafe(vec3 in) {
//     if (in == vec3(0,0,0)) {
//         return in;
//     }
//     else {
//         return normalize(in);
//     }
// }

void Rigidbody::applyForce(vec3 force, vec3 position) {
    vec3 scaledCenter = model->center * model->getScale();
    // cout << "scaledCenter = " << scaledCenter << endl;
    position += scaledCenter;
    // cout << "adjusted position = " << position << endl;
    if (position == vec3(0,0,0)) {
        velocity[3] += vec4(force, 0);
        return;
    }
    vec3 linearForce = dot(force, normalize(position)) * normalize(position);
    // cout << "linear force = " << linearForce << endl;
    //vec3 arbitraryOrthoganolVector = cross(force, position);
    velocity[3] += vec4(linearForce/mass, 0);

    cout << "linearAcc = " << vec4(linearForce/mass, 1) << endl;

    vec3 angularForce = force - linearForce;
    vec3 torque = cross(position, angularForce);
    // cout << "torque = " << torque << endl;

    //where 1 is the side length ie extent
    mat3 tensor = ((mass * 1 * 1)/6) * mat3();

    mat4 conversionMat = model->transform;
    conversionMat[3] = vec4(0,0,0,1);

    //cout << "inverse = " << inverse(tensor) << endl;
    vec3 angAcc = inverse(tensor) * toVec3(conversionMat * vec4(torque, 1));
    //cout << "angAcc = " << angAcc << endl;
    mat4 rot = rotationFromEuler(angAcc);
    //cout << "rotation matrix = " << rot << endl;
    velocity *= transpose(rot);
    velocity[3][3] = 1;
    // std::cout << "Velocity = " << velocity << endl;
}