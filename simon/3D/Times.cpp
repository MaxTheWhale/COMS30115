#include "Times.hpp"
#include <sys/time.h>

float getDeltaTime();
float Times::delta;
long long Times::frameCount;
long long Times::lastFrameTime;

void Times::init() {
    lastFrameTime = getTime();
    frameCount = 0;
}

//returns milliseconds since epoch
unsigned long long Times::getTime()
{
    struct timeval tv;

    gettimeofday(&tv, nullptr);

    unsigned long long millisecondsSinceEpoch =
        (unsigned long long)(tv.tv_sec) * 1000 +
        (unsigned long long)(tv.tv_usec) / 1000;
    return millisecondsSinceEpoch;
}
void Times::update() {
    delta = getDeltaTime();
    lastFrameTime = getTime();
    frameCount++;
}
float Times::deltaTime() {
    return delta;
}

float Times::getDeltaTime()
{
    long long ms = getTime() - lastFrameTime;
    return (float)ms / (float)1000;
}