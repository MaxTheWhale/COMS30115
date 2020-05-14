#include "Rigidbody.hpp"
#include <iostream>
#include <algorithm>

using namespace glm;

// std::vector<Rigidbody*> Rigidbody::allRBs = std::vector<Rigidbody*>();
std::vector<Rigidbody*>& Rigidbody::getAllRBs() {
    static std::vector<Rigidbody*> allRBs = std::vector<Rigidbody*>();
    return allRBs;
}

const mat3 Rigidbody::collisionLayers = glm::mat3(1,1,1,
                                                  1,0,0,
                                                  1,0,0);

Rigidbody::Rigidbody(Model* model) {
    getAllRBs().push_back(this);
    this->model = model;
    this->velocity = mat4(1,0,0,0,
                               0,1,0,0,
                               0,0,1,0,
                               0,0,0,1);
}

glm::vec3 Rigidbody::gravity = glm::vec3(0, -0.9f, 0);


void Rigidbody::update() {
    if (positionFixed) {
        return;
    }
    //if we're not moving and we're in contact with something then it can be assumed that we are resting on it
    if (hasGravity){// && !(velocity[3] == vec4(0,0,0,1) && !collidedWith.empty())) {
        vec3 grav = gravity * timeStep();
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

    model->transform *= velocity;// * timeStep();
    //disgusting hack to make linear velocity independent of angular
    model->transform[3] = oldTransform[3] + velocity[3];
    model->transform[3][3] = 1;

    for(unsigned int i = 0; i < getAllRBs().size(); i++) {
        if (getAllRBs()[i]->model != this->model) {
            if(collide(*getAllRBs()[i])) {
                // hasGravity = false;
                // velocity = mat4();
                model->transform = oldTransform;
            }
        }
    }
    // cout << "velocity = " << velocity << endl;
    // cout << "transform = " << model->transform << endl;
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

//so this c++ compiler implements a mod operator that doesn't flip negative numbers positive
//this one does
int positiveMod(int value, int mod) {
    if (value >= 0) {
        return value % mod;
    }
    else {
        return (value % mod) + mod;
    }
}

//returns the index of the element in the array that has a different sign to the other two
int differentSign(float values[3]) {
    for (int i = 0; i < 3; i++) {
        int prev = (i + 1) % 3;
        int next = positiveMod(i - 1, 3);
        if ((values[i] >= 0 && values[next] < 0 && values[next] < 0) || 
        (values[i] < 0 && values[prev] >=0 && values[prev] >= 0)) {
            return i;
        }
    }
    return -1;
}

//tests if the intervals a to b and x to y overlap on a line
bool intervalOverlap(float a, float b, float x, float y) {
    //fully contained
    if ((std::min(a,b) > std::min(x,y) && std::max(a,b) < std::max(x,y)) ||
        (std::min(x,y) > std::min(a,b) && std::max(x,y) < std::max(a,b))) {
        return false;
    }
    return std::max(a,b) >= std::min(x, y) && std::min(a,b) <= std::max(x, y);
}

float intervalMid(float a, float b, float x, float y) {
    if (std::max(a,b) < std::max(x,y)) {
        return (std::max(a,b) + std::min(x,y)) / 2.0f;
    }
    else {
        return (std::max(x,y) + std::min(a,b)) / 2.0f;
    }
}

bool Rigidbody::intersection(ModelTriangle localTri, ModelTriangle otherTri, mat4 localTransform, mat4 otherTransform) {
    vec3 verts1[3];
    vec3 verts2[3];
    // cout << "verts1 = ";
    for (int i = 0; i < 3; i++) {
        verts1[i] = toVec3(localTransform * localTri.vertices[i]);
        // cout << verts1[i] << ", ";
        verts2[i] = toVec3(otherTransform * otherTri.vertices[i]);
    }
    // cout << endl << "verts2 = ";
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

    //now we have to check which point is on the opposite side of the plane to the other two
    int oppositeIndexOne = differentSign(dist1);
    int one0 = (oppositeIndexOne + 1) % 3;
    int one2 = positiveMod(oppositeIndexOne - 1, 3);

    int oppositeIndexTwo = differentSign(dist2);
    int two0 = (oppositeIndexTwo + 1) % 3;
    int two2 = positiveMod(oppositeIndexTwo - 1, 3);

    // cout << "vertex indices: " << oppositeIndex << " " << other0 << " " << other2 << endl;

    vec3 D = cross(n1, n2);
    // cout << "D = "<< D << endl;
    float t1 = calcInterval(D, verts1[one0], verts1[oppositeIndexOne], dist1[one0], dist1[oppositeIndexOne]);
    float t2 = calcInterval(D, verts1[oppositeIndexOne], verts1[one2], dist1[oppositeIndexOne], dist1[one2]);
    float tPrime1 = calcInterval(D, verts2[two0], verts2[oppositeIndexTwo], dist2[two0], dist2[oppositeIndexTwo]);
    float tPrime2 = calcInterval(D, verts2[oppositeIndexTwo], verts2[two2], dist2[oppositeIndexTwo], dist2[two2]);
    //std::cout << "interval1: " << t1 << "," << t2 << " interval2: " << tPrime1 << "," << tPrime2 << endl;
    if (intervalOverlap(t1,t2,tPrime1,tPrime2)) {
        // cout << "t1,t2 = (" << t1 << "," << t2 << ") tPrime1, tPrime2 = (" << tPrime1 << "," << tPrime2 << ")" << endl;
        this->lastCollision = intervalMid(t1,t2,tPrime1,tPrime2) * normalize(D);
        // this->lastCollision = verts1[oppositeIndexOne];
        return true;
    }
    //std::cout << "default" << endl;
    return false;
}

bool isBasis(vec4 vec) {
    if (glm::length(vec) == 1) {
        // cout << "length 1" << endl;
        for (int i = 0; i < 3; i++) {
            // cout << abs(vec[i]) << endl;
            if (abs(vec[i]) == 1) {
                return true;
            }
        }
    }
    return false;
}

bool Rigidbody::collide(Rigidbody other) {
    if (positionFixed || !collisionEnabled || !other.collisionEnabled ||
        collisionLayers[collisionLayer][other.collisionLayer] == 0) {
        return false;
    }
    float dist = glm::distance(toVec3(this->model->transform[3]), toVec3(other.model->transform[3]));
    float maxDist = this->model->furthestExtent + other.model->furthestExtent;
    // cout << "dist = " << dist << endl;
    // cout << "maxDist = " << maxDist << endl;
    if (dist > maxDist) { //no possible collisions, too far
        return false;
    }
    for (unsigned int i = 0; i < model->tris.size(); i++) {
        for (unsigned int j = 0; j < other.model->tris.size(); j++) {
            //this line is disgusting because C++'s iterator interface hurts my soul
            if (std::find_if(collidedWith.begin(), collidedWith.end(), [other](Rigidbody* rb) {return &other == rb;}) == collidedWith.end() &&
                    intersection(model->tris[i], other.model->tris[j], model->transform, other.model->transform)) {
                //if (velocity[3] != vec4(0,0,0,1)) { //object is resting in contact with other object, so no reaction force is applied
                    vec4 normal = other.model->tris[j].normal * other.model->transform;
                    if (isBasis(normal)) {
                        continue;
                    }
                    if (glm::length(this->lastCollision) > this->model->furthestExtent) { //collision apparently happened outside the bounds of the model
                        // cout << "illegal collision point = " << this->lastCollision << endl;
                        continue;
                    }
                    // std::cout << "collision detected with " << &other << std::endl;
                    // cout << "normal = " << normal << endl;
                    float combinedElasticity = (this->elasticity + other.elasticity); //the average ratio of energy conserved
                    // cout << "total elasticity = " << combinedElasticity << endl;
                    vec3 force = toVec3(vec4() - combinedElasticity * dot(normalize(normal), normalize(velocity[3])) * normal);
                    // cout << "force = " << force << endl;
                    // cout << "collision point = " << this->lastCollision << endl;
                    applyForce(force, vec3(0,0,0));
                    // applyForce(force, this->lastCollision);
                    // this->positionFixed = true;
                    // this->collisionEnabled = false; //only collide once for testing purposes
                    // for (int i = 0; i < 3; i++) {
                    //     if (velocity[3][i] < 0.01f) {
                    //         velocity[3][i] = 0;
                    //     }
                    // }
                    // cout << endl;
                //}
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
    // vec3 scaledCenter = model->center * model->getScale();
    // cout << "scaledCenter = " << scaledCenter << endl;
    // position += scaledCenter;
    // cout << "adjusted position = " << position << endl;
    if (position == vec3(0,0,0)) {
        velocity[3] += vec4(force, 0);
        return;
    }
    cout << "normalised position = " << normalize(position) << endl;
    vec3 linearForce = dot(force, normalize(position)) * normalize(position);
    // vec3 linearForce = normalize(position) * glm::length(force) / glm::length(position);
    cout << "linear force = " << linearForce << endl;
    //vec3 arbitraryOrthoganolVector = cross(force, position);
    velocity[3] += vec4(linearForce/mass, 0);

    cout << "linearAcc = " << vec4(linearForce/mass, 0) << endl;

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

void Rigidbody::applyForce(vec3 force) {
    applyForce(force, vec3(0,0,0));
}