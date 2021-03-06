#include <ModelTriangle.h>
#include <CanvasTriangle.h>
#include <DrawingWindow.h>
#include <Utils.h>
#include <fstream>
#include <glm/glm.hpp>
#include <vector>

using namespace std;
using namespace glm;

#define WIDTH 320
#define HEIGHT 240

void draw();
void line(CanvasPoint p, CanvasPoint q, int colour);
void triangle(CanvasTriangle t, int colour, bool filled = false);
int *loadPPM(string fileName, int& width, int& height);
void skipHashWS(ifstream& f);
void update();
void handleEvent(SDL_Event event);
vector<float> Interpolate(float a, float b, int n);
vector<vec3> Interpolate(vec3 a, vec3 b, int n);

class Texture
{
  public:
    int width, height;
    int *buff;
    Texture(string fileName)
    {
      buff = loadPPM(fileName, width, height);
    }

    void draw(DrawingWindow window)
    {
      for (int y = 0; y < height; y++)
      {
        for (int x = 0; x < width; x++)
        {
          window.setPixelColour(x, y, buff[y * width + x]);
        }
      }
    }
};

void texturedTriangle(CanvasTriangle screenTri, CanvasTriangle texTri, Texture tex);

DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
Texture img = Texture("texture.ppm");

int main(int argc, char *argv[]) {
  SDL_Event event;
  CanvasTriangle tTri = CanvasTriangle(CanvasPoint(195, 5), CanvasPoint(395, 380), CanvasPoint(65, 330));
  CanvasTriangle sTri = CanvasTriangle(CanvasPoint(160, 10), CanvasPoint(300, 230), CanvasPoint(10, 150));
  while (true) {
    // We MUST poll for events - otherwise the window will freeze !
    if (window.pollForInputEvents(&event))
      handleEvent(event);
    update();
    draw();
    // Need to render the frame at the end, or nothing actually gets shown on
    // the screen !
    window.renderFrame();
  }
}

void draw() {
}

void update() {
  // Function for performing animation (shifting artifacts or moving the camera)
}

void handleEvent(SDL_Event event) {
  if (event.type == SDL_KEYDOWN) {
    if (event.key.keysym.sym == SDLK_LEFT)
      cout << "LEFT" << endl;
    else if (event.key.keysym.sym == SDLK_RIGHT)
      cout << "RIGHT" << endl;
    else if (event.key.keysym.sym == SDLK_UP)
      cout << "UP" << endl;
    else if (event.key.keysym.sym == SDLK_DOWN)
      cout << "DOWN" << endl;
    else if (event.key.keysym.sym == SDLK_u) {
      cout << "U" << endl;
      triangle(CanvasTriangle(CanvasPoint(rand() % WIDTH, rand() % HEIGHT), CanvasPoint(rand() % WIDTH, rand() % HEIGHT), CanvasPoint(rand() % WIDTH, rand() % HEIGHT)), rand());
    }
    else if (event.key.keysym.sym == SDLK_f) {
      cout << "F" << endl;
      triangle(CanvasTriangle(CanvasPoint(rand() % WIDTH, rand() % HEIGHT), CanvasPoint(rand() % WIDTH, rand() % HEIGHT), CanvasPoint(rand() % WIDTH, rand() % HEIGHT)), rand(), true);
    }
    else if (event.key.keysym.sym == SDLK_t) {
      cout << "T" << endl;
      texturedTriangle(CanvasTriangle(CanvasPoint(rand() % WIDTH, rand() % HEIGHT), CanvasPoint(rand() % WIDTH, rand() % HEIGHT), CanvasPoint(rand() % WIDTH, rand() % HEIGHT)), CanvasTriangle(CanvasPoint(rand() % img.width, rand() % img.height), CanvasPoint(rand() % img.width, rand() % img.height), CanvasPoint(rand() % img.width, rand() % img.height)), img);
    }
  } else if (event.type == SDL_MOUSEBUTTONDOWN)
    cout << "MOUSE CLICKED" << endl;
}

vector<float> Interpolate(float a, float b, int n) {
  vector<float> result;
  if (n == 1) {
    result.push_back((a + b) / 2);
    return result;
  } else if (n > 1) {
    float step = (b - a) / (n - 1);
    for (int i = 0; i < n; ++i) {
      result.push_back(a);
      a += step;
    }
  }
  return result;
}

vector<vec3> Interpolate(vec3 a, vec3 b, int n) {
  vector<vec3> result;
  if (n == 1) {
    result.push_back((a + b) / 2.0f);
    return result;
  } else if (n > 1) {
    vec3 step = (b - a) / (n - 1.0f);
    for (int i = 0; i < n; ++i) {
      result.push_back(a);
      a += step;
    }
  }
  return result;
}

void line(CanvasPoint p, CanvasPoint q, int colour) {
  float x_diff = p.x - q.x;
  float y_diff = p.y - q.y;
  float max_diff = std::max(abs(x_diff), abs(y_diff));
  vector<float> interpolate_x = Interpolate(p.x, q.x, max_diff + 1);
  vector<float> interpolate_y = Interpolate(p.y, q.y, max_diff + 1);
  for (int i = 0; i < max_diff; i++) {
    window.setPixelColour(round(interpolate_x[i]), round(interpolate_y[i]), colour);
  }
}

void triangle(CanvasTriangle t, int colour, bool filled) {
  if (filled) {
    if(t.vertices[0].y > t.vertices[1].y) swap(t.vertices[0], t.vertices[1]);
    if(t.vertices[1].y > t.vertices[2].y) swap(t.vertices[1], t.vertices[2]);
    if(t.vertices[0].y > t.vertices[1].y) swap(t.vertices[0], t.vertices[1]);

    vector<float> line1 = Interpolate(t.vertices[0].x, t.vertices[1].x, abs(t.vertices[0].y - t.vertices[1].y) + 1);
    vector<float> line2 = Interpolate(t.vertices[0].x, t.vertices[2].x, abs(t.vertices[0].y - t.vertices[2].y) + 1);
    vector<float> line3 = Interpolate(t.vertices[1].x, t.vertices[2].x, abs(t.vertices[1].y - t.vertices[2].y) + 1);

    for (int y = t.vertices[0].y; y < t.vertices[1].y; y++) {
      for (int x = std::min(line1[y - t.vertices[0].y], line2[y - t.vertices[0].y]); x < std::max(line1[y - t.vertices[0].y], line2[y - t.vertices[0].y]); x++) {
        window.setPixelColour(x, y, colour);
      }
    }
    for (int y = t.vertices[1].y; y < t.vertices[2].y; y++) {
      for (int x = std::min(line3[y - t.vertices[1].y], line2[y - t.vertices[0].y]); x < std::max(line3[y - t.vertices[1].y], line2[y - t.vertices[0].y]); x++) {
        window.setPixelColour(x, y, colour);
      }
    }
  }

  line(t.vertices[0], t.vertices[1], colour);
  line(t.vertices[1], t.vertices[2], colour);
  line(t.vertices[2], t.vertices[0], colour);
}

void texturedTriangle(CanvasTriangle screenTri, CanvasTriangle texTri, Texture tex) {
  if(screenTri.vertices[0].y > screenTri.vertices[1].y) {
    swap(screenTri.vertices[0], screenTri.vertices[1]);
    swap(texTri.vertices[0], texTri.vertices[1]);
  }
  if(screenTri.vertices[1].y > screenTri.vertices[2].y) {
    swap(screenTri.vertices[1], screenTri.vertices[2]);
    swap(texTri.vertices[1], texTri.vertices[2]);
  }
  if(screenTri.vertices[0].y > screenTri.vertices[1].y) {
    swap(screenTri.vertices[0], screenTri.vertices[1]);
    swap(texTri.vertices[0], texTri.vertices[1]);
  }

  vector<float> sline1 = Interpolate(screenTri.vertices[0].x, screenTri.vertices[1].x, abs(screenTri.vertices[0].y - screenTri.vertices[1].y) + 1);
  vector<float> sline2 = Interpolate(screenTri.vertices[0].x, screenTri.vertices[2].x, abs(screenTri.vertices[0].y - screenTri.vertices[2].y) + 1);
  vector<float> sline3 = Interpolate(screenTri.vertices[1].x, screenTri.vertices[2].x, abs(screenTri.vertices[1].y - screenTri.vertices[2].y) + 1);
  vector<float> txline1 = Interpolate(texTri.vertices[0].x, texTri.vertices[1].x, abs(screenTri.vertices[0].y - screenTri.vertices[1].y) + 1);
  vector<float> txline2 = Interpolate(texTri.vertices[0].x, texTri.vertices[2].x, abs(screenTri.vertices[0].y - screenTri.vertices[2].y) + 1);
  vector<float> txline3 = Interpolate(texTri.vertices[1].x, texTri.vertices[2].x, abs(screenTri.vertices[1].y - screenTri.vertices[2].y) + 1);
  vector<float> tyline1 = Interpolate(texTri.vertices[0].y, texTri.vertices[1].y, abs(screenTri.vertices[0].y - screenTri.vertices[1].y) + 1);
  vector<float> tyline2 = Interpolate(texTri.vertices[0].y, texTri.vertices[2].y, abs(screenTri.vertices[0].y - screenTri.vertices[2].y) + 1);
  vector<float> tyline3 = Interpolate(texTri.vertices[1].y, texTri.vertices[2].y, abs(screenTri.vertices[1].y - screenTri.vertices[2].y) + 1);

  for (int y = screenTri.vertices[0].y; y < screenTri.vertices[1].y; y++) {
    float tx1 = txline1[y - screenTri.vertices[0].y];
    float ty1 = tyline1[y - screenTri.vertices[0].y];
    float tx2 = txline2[y - screenTri.vertices[0].y];
    float ty2 = tyline2[y - screenTri.vertices[0].y];
    if (sline1[y - screenTri.vertices[0].y] > sline2[y - screenTri.vertices[0].y]) {
      swap(tx1, tx2);
      swap(ty1, ty2);
    }
    int xStart = std::min(sline1[y - screenTri.vertices[0].y], sline2[y - screenTri.vertices[0].y]);
    int xEnd = std::max(sline1[y - screenTri.vertices[0].y], sline2[y - screenTri.vertices[0].y]);
    vector<float> txs = Interpolate(tx1, tx2, xEnd - xStart);
    vector<float> tys = Interpolate(ty1, ty2, xEnd - xStart);
    for (int x = xStart; x < xEnd; x++) {
      window.setPixelColour(x, y, tex.buff[(int)tys[x - xStart] * tex.width + (int)txs[x - xStart]]);
    }
  }
  for (int y = screenTri.vertices[1].y; y < screenTri.vertices[2].y; y++) {
    float tx1 = txline2[y - screenTri.vertices[0].y];
    float ty1 = tyline2[y - screenTri.vertices[0].y];
    float tx2 = txline3[y - screenTri.vertices[1].y];
    float ty2 = tyline3[y - screenTri.vertices[1].y];
    if (sline2[y - screenTri.vertices[0].y] > sline3[y - screenTri.vertices[1].y]) {
      swap(tx1, tx2);
      swap(ty1, ty2);
    }
    int xStart = std::min(sline3[y - screenTri.vertices[1].y], sline2[y - screenTri.vertices[0].y]);
    int xEnd = std::max(sline3[y - screenTri.vertices[1].y], sline2[y - screenTri.vertices[0].y]);
    vector<float> txs = Interpolate(tx1, tx2, xEnd - xStart);
    vector<float> tys = Interpolate(ty1, ty2, xEnd - xStart);
    for (int x = xStart; x < xEnd; x++) {
      window.setPixelColour(x, y, tex.buff[(int)tys[x - xStart] * tex.width + (int)txs[x - xStart]]);
    }
  }
}

int *loadPPM(string fileName, int& width, int& height) {
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

  int* buff = new int[width * height];
  for (int i = 0; i < width * height; i++) {
    buff[i] = 0xff000000;
    f.read((char*)&buff[i], 3);
  }
  return buff;
}

void skipHashWS(ifstream& f) {
  ws(f);
  char current;
  f.read(&current, 1);
  if (current == '#') {
    f.ignore(1000, '\n');
    ws(f);
  }
  else {
    f.seekg(-1, f.cur);
  }
}