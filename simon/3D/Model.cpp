#include "Model.hpp"
#include <Utils.h>
#include <fstream>

using namespace std;

Model::Model(string filename) {
    palette = loadMTL(filename + ".mtl", texture.data, texture.width, texture.height);
    tris = loadOBJ(filename + ".obj", palette);
}

vector<ModelTriangle> Model::loadOBJ(string fileName,
                              unordered_map<string, Colour> palette)
{
  ifstream f;
  string s;
  string name = "";
  Colour colour = Colour(255, 255, 255);
  vector<glm::vec3> vertices;
  vector<glm::vec2> uvs;
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
        name = s;
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
      if (s == "vt")
      {
        float u, v;
        f >> u >> v;
        uvs.push_back(glm::vec2(u, v));
      }
      if (s == "f")
      {
        string a, b, c;
        f >> a >> b >> c;
        faces.push_back(ModelTriangle(vertices[stoi(split(a, '/')[0]) - 1],
                                      vertices[stoi(split(b, '/')[0]) - 1],
                                      vertices[stoi(split(c, '/')[0]) - 1],
                                      colour, name));
      }
    }
  }
  return faces;
}

void skipHashWS(ifstream &f)
{
  ws(f);
  char current;
  f.read(&current, 1);
  if (current == '#')
  {
    f.ignore(1000, '\n');
    ws(f);
  }
  else
  {
    f.seekg(-1, f.cur);
  }
}

int *loadPPM(string fileName, int &width, int &height)
{
  ifstream f;
  string s;
  f.open(fileName, ios::in | ios::binary);
  f >> s;
  skipHashWS(f);
  f >> s;
  width = stoi(s);
  skipHashWS(f);
  f >> s;
  height = stoi(s);
  skipHashWS(f);
  f >> s;
  f.seekg(1, f.cur);

  int *buff = new int[width * height];
  for (int i = 0; i < width * height; i++)
  {
    buff[i] = 0xff000000;
    f.read((char *)&buff[i], 3);
  }
  return buff;
}

unordered_map<string, Colour> Model::loadMTL(string fileName, int*& data, int& width, int& height) {
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
    if (s == "map_Kd") {
      string texture_file;
      f >> texture_file;
      data = loadPPM(texture_file, width, height);
    }
  }
  return palette;
}