#include "Lines.hpp"

using namespace glm;

enum CLIP_CODE {TOP = 1, RIGHT = 2, BOTTOM = 4, LEFT = 8};

// Line clipping based on https://www.geeksforgeeks.org/line-clipping-set-1-cohen-sutherland-algorithm/
int clipCode(vec4& p, int width, int height) {
  int code = 0;
  if (p.x < 0.0f) code |= LEFT;
  if (p.y < 0.0f) code |= TOP;
  if (p.x > width - 1.0f) code |= RIGHT;
  if (p.y > height - 1.0f) code |= BOTTOM;
  return code;
}

bool clipLine(vec4& p, vec4& q, int width, int height) {
  int p_code = clipCode(p, width, height);
  int q_code = clipCode(q, width, height);
  while (true) {
    if (p_code == 0 && q_code == 0) return true;
    if (p_code & q_code) return false;
    int current_code = p_code ? p_code : q_code;
    float x = 0.0f, y = 0.0f;
    if (current_code & TOP) {
      x = mix(p.x, q.x, (-p.y) / (q.y - p.y));
      y = 0.0f;
    }
    else if (current_code & RIGHT) {
      y = mix(p.y, q.y, ((width - 1.0f) - p.x) / (q.x - p.x));
      x = width - 1.0f;
    }
    else if (current_code & BOTTOM) {
      x = mix(p.x, q.x, ((height - 1.0f) - p.y) / (q.y - p.y));
      y = height - 1.0f;
    }
    else if (current_code & LEFT) {
      y = mix(p.y, q.y, (-p.x) / (q.x - p.x));
      x = 0.0f;
    }
    if (current_code == p_code) {
      p.x = x;
      p.y = y;
      p_code = clipCode(p, width, height);
    }
    else {
      q.x = x; 
      q.y = y; 
      q_code = clipCode(q, width, height);
    }
  }
}

void line(glm::vec4 p, glm::vec4 q, int colour, uint32_t *buffer, int width, int height, glm::vec2 &offset)
{
  if (!clipLine(p, q, width, height)) return;
  float x_diff = p.x - q.x;
  float y_diff = p.y - q.y;
  float max_diff = max(abs(x_diff), abs(y_diff));
  float step_x = (q.x - p.x) / max_diff;
  float step_y = (q.y - p.y) / max_diff;
  p.x += offset.x;
  p.y += offset.y;
  for (int i = 0; i <= max_diff; i++)
  {
    buffer[(int)(p.y) * width + (int)p.x] = colour;
    p.x += step_x;
    p.y += step_y;
  }
}