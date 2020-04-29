#pragma once

#include <stack>

#include "Movement.hpp"
#include "Transformable.hpp"
#include "Times.hpp"
#include "VectorUtil.hpp"
#include "Updatable.hpp"
#include "Temporal.hpp"

using namespace std;

class Animatable : public Transformable, public Updatable, Temporal {
    public:
        Animatable();
        stack<Movement> moves;
        void update() override;
    protected:
        Movement previous;
        bool firstUpdate = true;
};