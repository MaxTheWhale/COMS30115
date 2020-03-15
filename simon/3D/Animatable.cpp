#include "Animatable.hpp"
#include <glm/glm.hpp>
#include <iostream>
#include "VectorOutput.hpp"

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
    if (!this->moves.empty())
    {
        //cout << "delta = deltaTime (" << Times::deltaTime() << ") * (previous transform (" << this->previous.transform << ") - next transform (" << moves.top().transform << "))" << endl;
        mat4 delta = -Times::deltaTime() / moves.top().time * (this->previous.transform - moves.top().transform);
        //cout << "delta = " << delta << endl;
        mat4 newTransform = mat4();
        for (int i = 0; i < 4; i++)
        {
            newTransform[i] = this->transform[i] + delta[i];
        }
        //would the move take us further from our goal
        float currentDist = mat4Dist(this->transform, moves.top().transform);
        float newDist = mat4Dist(newTransform, moves.top().transform);
        cout << "currentDist = " << currentDist << " newDist = " << newDist << endl;
        if (currentDist < newDist)
        {
            this->previous = this->moves.top();
            this->moves.pop();
        }
        else
        {
            for (int i = 0; i < 4; i++)
            {
                this->transform[i] += delta[i];
            }
        }
    }
}