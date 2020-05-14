#pragma once

#include <glm/glm.hpp>
#include "Model.hpp"
#include "Updatable.hpp"
#include "Temporal.hpp"

class Rigidbody : public Updatable, Temporal {
    public:
        glm::mat4 velocity;
        void update() override;
        bool collide(Rigidbody other);
        bool hasGravity = true;
        bool collisionEnabled = true;
        bool positionFixed = true;
        static glm::vec3 gravity;
        Rigidbody(Model* model, vector<Rigidbody*>& RBList);
        Rigidbody();
        Model* model;
        std::vector<Rigidbody*> collidedWith; //all RBs that this RB has collided with this frame, ensuring that they do not collide twice
        float mass = 1;
        void applyForce(vec3 force);
        void applyForce(vec3 force, vec3 position);
        float elasticity = 0.9; //how much energy is conserved in collisions
        bool suckable = true; //can be moved by magnets
        std::vector<Rigidbody*>* allRBs;
        int collisionLayer = 0;
        int maxCollisions = -1;
    protected:
        static const glm::mat3 collisionLayers;
        //hack for the sake of efficiency
        //set every time there is a collision to the model-relative 3D coordinates of the vertex assumed to have caused it
        glm::vec3 lastCollision;
        bool intersection(ModelTriangle localTri, ModelTriangle otherTri, mat4 localTransform, mat4 otherTransform);
        int collisions = 0;
};