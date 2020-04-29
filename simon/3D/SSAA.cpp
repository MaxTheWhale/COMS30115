#include "SSAA.hpp"

enum COLOUR_MASK {ALPHA = 0xff000000, RED = 0x00ff0000, GREEN = 0x0000ff00, BLUE = 0x000000ff};

std::vector<glm::vec2> generateRotatedGrid(int gridSize) {
  std::vector<glm::vec2> result;
  const int numSamples = gridSize * gridSize;
  const float step = 1.0f / numSamples;
  float x = (step * 0.5f) + step * (numSamples - gridSize);
  for (float y = step * 0.5f; y < 1.0f; y += step) {
    result.push_back(glm::vec2(x, y));
    x -= step * gridSize;
    if (x < 0)
      x += step * (numSamples + 1);
  }
  return result;
}

void downsample(uint32_t *source, uint32_t *dest, int width, int height, int num_samples) {
  uint32_t pixel, sample, red, green, blue;
  for (int y = 0; y < height; y++) {
    for (int x = 0; x < width; x++) {
      red = 0;
      green = 0;
      blue = 0;
      for (int s = 0; s < num_samples; s++) {
        sample = source[(s * width * height) + x + y * width];
        red += (sample & RED);
        green += (sample & GREEN);
        blue += (sample & BLUE);
      }
      pixel = ALPHA | ((red / num_samples) & RED) | ((green / num_samples) & GREEN) | ((blue / num_samples) & BLUE);
      dest[x + y * width] = pixel;
    }
  }
}
