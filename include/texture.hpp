#ifndef TEXTURE_H
#define TEXTURE_H

#include <Vector3f.h>
#include <iostream>
#include "perlin.hpp"

using namespace std;
inline double clamp(double x, double min, double max)
{
	if (x < min)
		return min;
	if (x > max)
		return max;
	return x;
}
class texture
{
public:
	virtual Vector3f value(double u, double v, const Vector3f &p) const = 0;
};

class solid_color : public texture
{
public:
	solid_color() {}
	solid_color(Vector3f c) : color_value(c) {}
	solid_color(double red, double green, double blue) : solid_color(Vector3f(red, green, blue)) {}
	virtual Vector3f value(double u, double v, const Vector3f &p) const override
	{
		return color_value;
	}
	Vector3f color_value;
};

class checker_texture : public texture
{
public:
	checker_texture() {}
	checker_texture(texture *_even, texture *_odd) : even(_even), odd(_odd) {}
	checker_texture(Vector3f c1, Vector3f c2)
	{
		even = new solid_color(c1);
		odd = new solid_color(c2);
	}
	virtual Vector3f value(double u, double v, const Vector3f &p) const override
	{
		auto sines = sin(10 * p.x()) * sin(10 * p.y()) * sin(10 * p.z());
		if (sines < 0)
		{
			return odd->value(u, v, p);
		}
		else
		{
			return even->value(u, v, p);
		}
	}
	texture *odd;
	texture *even;
};

class image_texture : public texture
{
public:
	const static int bytes_per_pixel = 3;

	image_texture() : data(nullptr), width(0), height(0), bytes_per_scanline(0) {}
	image_texture(const char *filename);
	~image_texture() { delete data; }
	virtual Vector3f value(double u, double v, const Vector3f &p) const override
	{
		//如果没有颜色，返回青色
		if (data == nullptr)
			return Vector3f(0., 1., 1.);

		if(u<0) u=-u;
		if(v<0) v=-v;

		while(u>1.0) u-=1.0;
		while(v>1.0) v-=1.0;
		v = 1.0-v;

		auto i = static_cast<int>(u * width);
		auto j = static_cast<int>(v * height);
		if (i >= width)
			i = width - 1;
		if (j >= height)
			j = height - 1;

		const auto color_scale = 1.0 / 255.0;
		auto pixel = data + j * bytes_per_scanline + i * bytes_per_pixel;
		return Vector3f(color_scale * pixel[0], color_scale * pixel[1], color_scale * pixel[2]);
	}
	unsigned char *data;
	int width, height;
	int bytes_per_scanline;
};

class noise_texture : public texture
{
public:
	noise_texture(double sc) : scale(sc) {}
	noise_texture() {}
	Vector3f value(double u, double v, const Vector3f &vec) const override
	{
		return Vector3f(1, 1, 1) * noise.noise(scale * vec);
	}

	perlin noise;
	double scale;
};

#endif //TEXTURE_H