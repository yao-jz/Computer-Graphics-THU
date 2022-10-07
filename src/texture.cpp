#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <Vector3f.h>

#ifndef STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#endif // STB_IMAGE_IMPLEMENTATION
#include <iostream>
using namespace std;
#include "texture.hpp"
image_texture::image_texture(const char *filename)
{

	auto components_per_pixel = bytes_per_pixel;

	data = stbi_load(filename, &width, &height, &components_per_pixel, components_per_pixel);
	if (!data)
	{
		std::cerr << "ERROR: Could not load texture image file '" << filename << "'.\n";
		width = height = 0;
	}
	bytes_per_scanline = bytes_per_pixel * width;
}