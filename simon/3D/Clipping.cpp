#include "Clipping.hpp"
#include <vector>
#include <algorithm>

using namespace std;
using namespace glm;

inline Vertex mixVertex(const Vertex &p, const Vertex &q, float a) {
  return p * (1.0f - a) + q * a;
}

// based on https://casual-effects.com/research/McGuire2011Clipping/McGuire-Clipping.pdf
int clipTriangle(list<Triangle>& tris, const vec4& normal) {
  Vertex temp;
  list<Triangle>::iterator tri = tris.begin();
  int n = tris.size();
  for (int i = 0; (i < n) && tri != tris.end(); i++) {
    vector<float> distances;
    distances.push_back(dot((*tri).vertices[0].pos, normal));
    distances.push_back(dot((*tri).vertices[1].pos, normal));
    distances.push_back(dot((*tri).vertices[2].pos, normal));
    if (distances[0] >= 0.0f && distances[1] >= 0.0f && distances[2] >= 0.0f) {
      tri++;
      continue;
    }
    if (distances[0] < 0.0f && distances[1] < 0.0f && distances[2] < 0.0f) {
      tri = tris.erase(tri);
      continue;
    }
    bool nextInside;
    if (distances[1] >= 0.0f && distances[0] < 0.0f) {
      nextInside = (distances[2] >= 0.0f);
      temp = (*tri).vertices[0];
      (*tri).vertices[0] = (*tri).vertices[1];
      (*tri).vertices[1] = (*tri).vertices[2];
      (*tri).vertices[2] = temp;
      rotate(distances.begin(),distances.begin()+1,distances.end());
    }
    else if (distances[2] >= 0.0f && distances[1] < 0.0f) {
      nextInside = (distances[0] >= 0.0f);
      temp = (*tri).vertices[2];
      (*tri).vertices[2] = (*tri).vertices[1];
      (*tri).vertices[1] = (*tri).vertices[0];
      (*tri).vertices[0] = temp;
      rotate(distances.begin(),distances.begin()+2,distances.end());
    }
    else {
      nextInside = (distances[1] >= 0.0f);
    }
    temp = mixVertex((*tri).vertices[0], (*tri).vertices[2], (distances[0] / (distances[0] - distances[2])));
    if (nextInside) {
      (*tri).vertices[2] = mixVertex((*tri).vertices[1], (*tri).vertices[2], (distances[1] / (distances[1] - distances[2])));
      Triangle newTri = Triangle((*tri).vertices[0], (*tri).vertices[2], temp, (*tri).mat, (*tri).normal);
      tris.push_back(newTri);
    }
    else {
      (*tri).vertices[1] = mixVertex((*tri).vertices[0], (*tri).vertices[1], (distances[0] / (distances[0] - distances[1])));
      (*tri).vertices[2] = temp;
    }
    tri++;
  }
  return tris.size();
}

int clipToView(list<Triangle>& tris) {
  const vec4 normals[6] = {vec4(1, 0, 0, 1), vec4(-1, 0, 0, 1), vec4(0, 1, 0, 1), vec4(0, -1, 0, 1), vec4(0, 0, 1, 1), vec4(0, 0, -1, 1)};
  for (auto n : normals) {
    int num = clipTriangle(tris, n);
    if (num == 0) {
      return 0;
    }
  }
  return tris.size();
}