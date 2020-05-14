#include "Model.hpp"
#include <Utils.h>
#include <fstream>
#include "Material.h"
#include "VectorUtil.hpp"
#include <string>

using namespace std;

Model::Model(string filename) {
    palette = loadMTL(filename + ".mtl");
    tris = loadOBJ(filename + ".obj", palette);
    this->furthestExtent = calcExtent();
}

Model::Model(const Model& original) {
    palette = original.palette;
    tris = original.tris;
    furthestExtent = original.furthestExtent;
    center = original.center;
}

vector<ModelTriangle> Model::loadOBJ(string fileName,
                              unordered_map<string, Material> palette)
{
  cout << "loading OBJ " << fileName << endl;
  ifstream f;
  string s;
  string name = "";
  Material material = Material("missing", Colour(255, 255, 255), Colour(255, 255, 255), Colour(0, 0, 0));
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
        if (material.normal_map.dataVec != nullptr) {
          tri.tangent = calcTangent(tri);
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
  cout << "loading PPM " << fileName << endl;
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

unordered_map<string, Material> Model::loadMTL(string fileName) {
  unordered_map<string, Material> palette;

  ifstream f;
  string s;
  string key;
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
        float illum;
        f >> illum;
        palette[key].illum = (int)illum;
      }
      if (s == "d") {
        float dissolve;
        f >> dissolve;
        palette[key].dissolve = dissolve;
      }
      if (s == "map_Kd") {
        string texture_file;
        f >> texture_file;
        size_t pos = fileName.find('/');
        if (pos != fileName.npos) {
          fileName.erase(pos + 1);
          texture_file = fileName + texture_file;
        }
        palette[key].texture.data = loadPPM(texture_file, palette[key].texture.width, palette[key].texture.height);

        palette[key].texture.dataVec = new glm::vec3[palette[key].texture.width * palette[key].texture.height];
        for (int i = 0; i < palette[key].texture.width * palette[key].texture.height; i++) {
          palette[key].texture.dataVec[i].r = ((palette[key].texture.data[i] & 0xff0000) >> 16) / 255.0f;
          palette[key].texture.dataVec[i].g = ((palette[key].texture.data[i] & 0x00ff00) >> 8) / 255.0f;
          palette[key].texture.dataVec[i].b = (palette[key].texture.data[i] & 0x0000ff) / 255.0f;
        }
      }
      if (s == "map_bump" || s == "bump") {
        string texture_file;
        f >> texture_file;
        size_t pos = fileName.find('/');
        if (pos != fileName.npos) {
          fileName.erase(pos + 1);
          texture_file = fileName + texture_file;
        }
        palette[key].normal_map.data = loadPPM(texture_file, palette[key].normal_map.width, palette[key].normal_map.height);

        palette[key].normal_map.dataVec = new glm::vec3[palette[key].normal_map.width * palette[key].normal_map.height];
        for (int i = 0; i < palette[key].normal_map.width * palette[key].normal_map.height; i++) {
          palette[key].normal_map.dataVec[i].x = ((palette[key].normal_map.data[i] & 0xff0000) >> 16) / 255.0f;
          palette[key].normal_map.dataVec[i].y = ((palette[key].normal_map.data[i] & 0x00ff00) >> 8) / 255.0f;
          palette[key].normal_map.dataVec[i].z = (palette[key].normal_map.data[i] & 0x0000ff) / 255.0f;
          palette[key].normal_map.dataVec[i].x *= 2.0f;
          palette[key].normal_map.dataVec[i].x -= 1.0f;
          palette[key].normal_map.dataVec[i].y *= 2.0f;
          palette[key].normal_map.dataVec[i].y -= 1.0f;
          palette[key].normal_map.dataVec[i].z *= 2.0f;
          palette[key].normal_map.dataVec[i].z -= 1.0f;
          palette[key].normal_map.dataVec[i] = glm::normalize(palette[key].normal_map.dataVec[i]);
        }
      }
    }
  }
  return palette;
}

//adapted from https://stackoverflow.com/questions/2083771/a-method-to-calculate-the-centre-of-mass-from-a-stl-stereo-lithography-file
vec3 Model::centerOfMass() {
  float totalVolume = 0;
  vec3 center;
  for (unsigned int i = 0; i < tris.size(); i++) {
    //maaaaaaaaaths
    float currentVolume = (tris[i].vertices[0][0]*tris[i].vertices[1][1]*tris[i].vertices[2][2] -
                           tris[i].vertices[0][0]*tris[i].vertices[2][1]*tris[i].vertices[1][2] -
                           tris[i].vertices[1][0]*tris[i].vertices[0][1]*tris[i].vertices[2][2] +
                           tris[i].vertices[1][0]*tris[i].vertices[2][1]*tris[i].vertices[0][2] +
                           tris[i].vertices[2][0]*tris[i].vertices[0][1]*tris[i].vertices[1][2] -
                           tris[i].vertices[2][0]*tris[i].vertices[1][1]*tris[i].vertices[0][2]) / 6;
    totalVolume += currentVolume;
    center[0] += ((tris[i].vertices[0][0] + tris[i].vertices[1][0] + tris[i].vertices[2][0]) / 4) * currentVolume;
    center[1] += ((tris[i].vertices[0][1] + tris[i].vertices[1][1] + tris[i].vertices[2][1]) / 4) * currentVolume;
    center[2] += ((tris[i].vertices[0][2] + tris[i].vertices[1][2] + tris[i].vertices[2][2]) / 4) * currentVolume;
  }
  return center;
}

float Model::calcExtent() {
  float value = -1.0f;
  vec4 scales;
  for (int i = 0; i < 3; i++) {
    scales[i] = glm::length(this->transform[i]);
  }
  scales[3] = 0;
  for (unsigned int i = 0; i < this->tris.size(); i++) {
    for (unsigned int j = 0; j < 3; j++) {
      float dist = glm::length(this->tris[i].vertices[j] * scales);
      if (dist > value) {
        value = dist;
      }
    }
  }
  return value;
}

vec4 Model::calcTangent(ModelTriangle& tri) {
  vec4 edge1 = tri.vertices[1] - tri.vertices[0];
  vec4 edge2 = tri.vertices[2] - tri.vertices[0];
  vec2 deltaUV1 = tri.uvs[1] - tri.uvs[0];
  vec2 deltaUV2 = tri.uvs[2] - tri.uvs[0];
  float f = 1.0f / (deltaUV1.x * deltaUV2.y - deltaUV2.x * deltaUV1.y);

  vec4 tangent;
  tangent.x = f * (deltaUV2.y * edge1.x - deltaUV1.y * edge2.x);
  tangent.y = f * (deltaUV2.y * edge1.y - deltaUV1.y * edge2.y);
  tangent.z = f * (deltaUV2.y * edge1.z - deltaUV1.y * edge2.z);
  tangent.w = 0;
  return normalize(tangent);
}