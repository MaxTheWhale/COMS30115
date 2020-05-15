#pragma once

//provides a simple interface for any class which needs to do something once per frame
class Updatable {
    public:
        virtual void update() = 0;
};