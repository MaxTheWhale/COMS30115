#include "Vertex.hpp"

Vertex Vertex::operator+=(const Vertex& rhs)
{
    pos += rhs.pos;
    pos_3d += rhs.pos_3d;
    normal += rhs.normal;
    brightness += rhs.brightness;
    u += rhs.u;
    v += rhs.v;
    return *this;
}

Vertex operator+(Vertex lhs, const Vertex& rhs)
{
    lhs += rhs;
    return lhs;
}

Vertex Vertex::operator-=(const Vertex& rhs)
{
    pos -= rhs.pos;
    pos_3d -= rhs.pos_3d;
    normal -= rhs.normal;
    brightness -= rhs.brightness;
    u -= rhs.u;
    v -= rhs.v;
    return *this;
}

Vertex operator-(Vertex lhs, const Vertex& rhs)
{
    lhs -= rhs;
    return lhs;
}

Vertex Vertex::operator*=(float rhs)
{
    pos *= rhs;
    pos_3d *= rhs;
    normal *= rhs;
    brightness *= rhs;
    u *= rhs;
    v *= rhs;
    return *this;
}

Vertex operator*(Vertex lhs, float rhs)
{
    lhs *= rhs;
    return lhs;
}

Vertex Vertex::operator/=(float rhs)
{
    pos /= rhs;
    pos_3d /= rhs;
    normal /= rhs;
    brightness /= rhs;
    u /= rhs;
    v /= rhs;
    return *this;
}

Vertex operator/(Vertex lhs, float rhs)
{
    lhs /= rhs;
    return lhs;
}