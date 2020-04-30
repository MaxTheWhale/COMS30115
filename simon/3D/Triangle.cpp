#include "Triangle.hpp"

Triangle::Triangle(ModelTriangle &tri) {
    for (int i = 0; i < 3; i++) {
        vertices[i].pos.x = tri.vertices[i].x;
        vertices[i].pos.y = tri.vertices[i].y;
        vertices[i].pos.z = tri.vertices[i].z;
        vertices[i].pos.w = tri.vertices[i].w;
        vertices[i].u = tri.uvs[i].x;
        vertices[i].v = tri.uvs[i].y;
        vertices[i].brightness = tri.brightness[i];
        vertices[i].normal = tri.normals[i];
    }
    normal = tri.normal;
    tangent = tri.tangent;
    mat = tri.material;
}

Triangle::Triangle(const Vertex &v0, const Vertex &v1, const Vertex &v2, const Material &tMat, const glm::vec4 &tNormal, const glm::mat3 &tTBN) {
    vertices[0] = v0;
    vertices[1] = v1;
    vertices[2] = v2;
    mat = tMat;
    normal = tNormal;
    TBN = tTBN;
}