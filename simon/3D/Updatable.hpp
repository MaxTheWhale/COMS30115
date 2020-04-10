#pragma once
#include <iostream>

class Updatable {
    public:
        virtual void update() = 0;//{ std::cout << "Updatable update called, this means something is wrong with inheritence" << std::endl; };
};