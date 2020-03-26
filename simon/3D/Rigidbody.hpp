#include <glm/glm.hpp>
#include "Model.hpp"
#include "Updatable.hpp"

class Rigidbody : public Updatable {
    public:
        glm::mat4 velocity;
        void update() override;
        void collide(Rigidbody other);
        bool hasGravity = true;
        bool collisionEnabled = true;
        static glm::vec3 gravity;
        Rigidbody(Model* model);
        Model* model;
    protected:
        static std::vector<Rigidbody*> allRBs;
};