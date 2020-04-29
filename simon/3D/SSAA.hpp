#pragma once

#include <glm/glm.hpp>
#include <vector>

std::vector<glm::vec2> generateRotatedGrid(int gridSize);
void downsample(uint32_t *source, uint32_t *dest, int width, int height, int num_samples);