#include <glm/glm.hpp>
#include "Model.hpp"
#include "Updatable.hpp"

class Rigidbody : public Updatable {
    public:
        glm::mat4 velocity;
        void update() override;
        bool collide(Rigidbody other);
        bool hasGravity = true;
        bool collisionEnabled = true;
        bool positionFixed = true;
        static glm::vec3 gravity;
        const static bool realTimeScale = false;
        Rigidbody(Model* model);
        Model* model;
        std::vector<Rigidbody*> collidedWith; //all RBs that this RB has collided with this frame, ensuring that they do not collide twice
        float mass = 1;
        void applyForce(vec3 force, vec3 position);
        float elasticity = 0.9; //how much energy is conserved in collisions
    protected:
        static std::vector<Rigidbody*> allRBs;
};