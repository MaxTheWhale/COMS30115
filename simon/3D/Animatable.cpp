#include "Animatable.hpp"
#include <glm/glm.hpp>
#include <iostream>

using namespace glm;
using namespace std;

Animatable::Animatable()
{
    this->previous.transform = this->transform;
}

float mat4Dist(mat4 &a, mat4 &b)
{
    float result = 0;
    for (int i = 0; i < 4; i++)
    {
        result += distance(a[i], b[i]);
    }
    return result;
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
        //cout << "delta = deltaTime (" << Times::deltaTime() << ") * (previous transform (" << this->previous.transform << ") - next transform (" << moves.top().transform << "))" << endl;
        mat4 delta = -timeStep() / moves.top().time * (previous.transform - moves.top().transform);
        cout << "previous = " << previous.transform << endl << "target = " << moves.top().transform << endl << "delta = " << delta << endl;
        mat4 newTransform = mat4();
        for (int i = 0; i < 4; i++)
        {
            newTransform[i] = transform[i] + delta[i];
        }
        //would the move take us further from our goal
        float currentDist = mat4Dist(transform, moves.top().transform);
        float newDist = mat4Dist(newTransform, moves.top().transform);
        cout << "currentDist = " << currentDist << " newDist = " << newDist << endl;
        if (currentDist < newDist)
        {
            transform = moves.top().transform;
            previous = moves.top();
            moves.pop();
        }
        else
        {
            for (int i = 0; i < 4; i++)
            {
                transform[i] += delta[i];
            }
        }
    }
    else
    {
        previous.transform = transform;
    }
    
}