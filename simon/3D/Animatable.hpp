#pragma once

#include <stack>

#include "Movement.hpp"
#include "Transformable.hpp"
#include "Times.hpp"
#include "VectorOutput.hpp"

using namespace std;

class Animatable : public Transformable {
    public:
        Animatable();
        stack<Movement> moves;
        void update();
    protected:
        Movement previous;
        bool firstUpdate = true;
};