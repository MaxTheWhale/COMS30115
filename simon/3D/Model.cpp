#include "Model.hpp"
#include <Utils.h>
#include <fstream>

using namespace std;

Model::Model(string filename) {
    palette = loadMTL(filename + ".mtl");
    tris = loadOBJ(filename + ".obj", palette);
}

vector<ModelTriangle> Model::loadOBJ(string fileName,
                              unordered_map<string, Colour> palette)
{
  ifstream f;
  string s;
  Colour colour;
  vector<glm::vec3> vertices;
  vector<ModelTriangle> faces;
  f.open(fileName, ios::in);
  while (!f.eof())
  {
    if (f >> s)
    {
      if (s == "mtllib")
      {
        f >> s;
      }
      if (s == "o")
      {
        f >> s;
      }
      if (s == "usemtl")
      {
        f >> s;
        colour = palette[s];
      }
      if (s == "v")
      {
        float x, y, z;
        f >> x >> y >> z;
        vertices.push_back(glm::vec3(x, y, z));
      }
      if (s == "f")
      {
        string a, b, c;
        f >> a >> b >> c;
        faces.push_back(ModelTriangle(vertices[stoi(split(a, '/')[0]) - 1],
                                      vertices[stoi(split(b, '/')[0]) - 1],
                                      vertices[stoi(split(c, '/')[0]) - 1],
                                      colour));
      }
    }
  }
  return faces;
}

unordered_map<string, Colour> Model::loadMTL(string fileName) {
  unordered_map<string, Colour> palette;

  ifstream f;
  string s;
  f.open(fileName, ios::in);
  if (!f.good()) {
    return palette;
  }
  while (!f.eof())
  {
    f >> s;
    if (s == "newmtl")
    {
      string key, r, g, b;
      f >> key;
      f >> s;
      f >> r;
      f >> g;
      f >> b;
      palette[key] = Colour(key, stof(r) * 255, stof(g) * 255, stof(b) * 255);
    }
  }
  return palette;
}