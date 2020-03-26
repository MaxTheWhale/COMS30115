#pragma once

#include <stack>

#include "Movement.hpp"
#include "Transformable.hpp"
#include "Times.hpp"
#include "VectorOutput.hpp"
#include "Updatable.hpp"

using namespace std;

class Animatable : public Transformable, public Updatable {
    public:
        Animatable();
        stack<Movement> moves;
        void update() override;
    protected:
        Movement previous;
        bool firstUpdate = true;
};