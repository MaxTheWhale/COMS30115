#include "Model.hpp"
#include <Utils.h>
#include <fstream>
#include "Material.h"

using namespace std;

Model::Model(string filename) {
    palette = loadMTL(filename + ".mtl", texture.data, texture.width, texture.height);
    tris = loadOBJ(filename + ".obj", palette);
    if (texture.data != nullptr) {
      texture.dataVec = new glm::vec3[texture.width * texture.height];
      for (int i = 0; i < texture.width * texture.height; i++) {
        texture.dataVec[i].r = ((texture.data[i] & 0xff0000) >> 16) / 255.0f;
        texture.dataVec[i].g = ((texture.data[i] & 0x00ff00) >> 8) / 255.0f;
        texture.dataVec[i].b = (texture.data[i] & 0x0000ff) / 255.0f;
      }
    }
}

vector<ModelTriangle> Model::loadOBJ(string fileName,
                              unordered_map<string, Material> palette)
{
  ifstream f;
  string s;
  string name = "";
  Material material = Material("missing", Colour(0, 0, 0), Colour(255, 255, 255), Colour(0, 0, 0));
  vector<glm::vec3> vertices;
  vector<glm::vec4> normals;
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
        if (s.find(".mtl") == s.npos)
          material = palette[s];
      }
      if (s == "v")
      {
        float x, y, z;
        f >> x >> y >> z;
        vertices.push_back(glm::vec3(x, y, z));
      }
      if (s == "vn")
      {
        float x, y, z;
        f >> x >> y >> z;
        normals.push_back(glm::vec4(x, y, z, 0.0f));
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
        int numFields = count(a.begin(), a.end(), '/') + 1;
        string *as = split(a, '/');
        string *bs = split(b, '/');
        string *cs = split(c, '/');
        ModelTriangle tri = ModelTriangle(vertices[stoi(as[0]) - 1],
                                      vertices[stoi(bs[0]) - 1],
                                      vertices[stoi(cs[0]) - 1],
                                      material, name);
        if (numFields > 1) {
          if (as[1] != "" && bs[1] != "" && cs[1] != "") {
            tri.uvs[0] = uvs[stoi(as[1]) - 1];
            tri.uvs[1] = uvs[stoi(bs[1]) - 1];
            tri.uvs[2] = uvs[stoi(cs[1]) - 1];
          }
        }
        if (numFields > 2) {
          if (as[2] != "" && bs[2] != "" && cs[2] != "") {
            tri.normals[0] = normals[stoi(as[2]) - 1];
            tri.normals[1] = normals[stoi(bs[2]) - 1];
            tri.normals[2] = normals[stoi(cs[2]) - 1];
          }
        }
        delete [] as;
        delete [] bs;
        delete [] cs;
        faces.push_back(tri);
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
    f.read(((char*)&buff[i]) + 2, 1);
    f.read(((char*)&buff[i]) + 1, 1);
    f.read((char*)&buff[i], 1);
  }
  return buff;
}

unordered_map<string, Material> Model::loadMTL(string fileName, int*& data, int& width, int& height) {
  unordered_map<string, Material> palette;

  ifstream f;
  string s;
  string key;
  data = nullptr;
  f.open(fileName, ios::in);
  if (!f.good()) {
    return palette;
  }
  while (!f.eof())
  {
    if (f >> s)
    {
      if (s == "newmtl")
      {
        f >> key;
        palette[key] = Material(key);
      }
      if (s == "Ka") {
        string r, g, b;
        f >> r;
        f >> g;
        f >> b;
        palette[key].ambientVec = glm::vec3(stof(r), stof(g), stof(b));
        palette[key].ambient = Colour(key, stof(r) * 255, stof(g) * 255, stof(b) * 255);
      }
      if (s == "Kd") {
        string r, g, b;
        f >> r;
        f >> g;
        f >> b;
        palette[key].diffuseVec = glm::vec3(stof(r), stof(g), stof(b));
        palette[key].diffuse = Colour(key, stof(r) * 255, stof(g) * 255, stof(b) * 255);
        if (palette[key].ambient.red == -1) {
          palette[key].ambient = palette[key].diffuse;
          palette[key].ambientVec = palette[key].diffuseVec;
        }
      }
      if (s == "Ks") {
        string r, g, b;
        f >> r;
        f >> g;
        f >> b;
        palette[key].specularVec = glm::vec3(stof(r), stof(g), stof(b));
        palette[key].specular = Colour(key, stof(r) * 255, stof(g) * 255, stof(b) * 255);
      }
      if (s == "Ns") {
        float Ns;
        f >> Ns;
        palette[key].highlights = Ns;
      }
      if (s == "illum") {
        int illum;
        f >> illum;
        palette[key].illum = illum;
      }
      if (s == "d") {
        f >> s;
        palette[key].dissolve = stoi(s);
      }
      if (s == "map_Kd") {
        string texture_file;
        f >> texture_file;
        size_t pos = fileName.find('/');
        if (pos != fileName.npos) {
          fileName.erase(pos + 1);
          texture_file = fileName + texture_file;
          cout << texture_file << '\n';
        }
        data = loadPPM(texture_file, width, height);
      }
    }
  }
  return palette;
}