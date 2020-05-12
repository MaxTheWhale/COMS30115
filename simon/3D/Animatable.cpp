#include "Animatable.hpp"
#include <glm/glm.hpp>
#include <iostream>

using namespace glm;
using namespace std;

Animatable::Animatable()
{
    this->previous.transform = this->transform;
}

void Animatable::update()
{
    //cout << "Animatable update, moves length = " << this->moves.size() << endl;
    if (firstUpdate)
    {
        previous.transform = transform;
        firstUpdate = false;
    }
    if (!moves.empty())
    {
       if (moves.top()->execute(this, this->previous)) {
           moves.top()->repeats--;
           previous = *moves.top();
           if (moves.top()->repeats == 0) {
               moves.pop();
           }
       }
    }
    else
    {
        previous.transform = transform;
    }
    
}