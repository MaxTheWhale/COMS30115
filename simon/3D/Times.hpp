#pragma once

class Times {
    public:
        static unsigned long long getTime();
        static float deltaTime();
        static void update();
        static void init();
    private:
        static float getDeltaTime();
        static float delta;
        static long long frameCount;
        static long long lastFrameTime;
        Times(){}
};