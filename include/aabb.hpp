#ifndef AABB_HPP
#define AABB_HPP

#include "ray.hpp"
#include <cmath>
#include <vecmath.h>
#include <iostream>
using namespace std;

// Axis-aligned bounding box
class aabb
{
public:
	aabb() {}
	aabb(const Vector3f &a, const Vector3f &b)
	{
		min = a;
		max = b;
	}

	inline bool hit(const Ray &r) const
	{
		double t_min = -1; // 初值如何取？
		double t_max = 1e6; 
		for (int a = 0; a < 3; a++)
		{
			auto t0 = fmin((min[a] - r.getOrigin()[a]) / r.getDirection()[a], (max[a] - r.getOrigin()[a]) / r.getDirection()[a]);
			auto t1 = fmax((min[a] - r.getOrigin()[a]) / r.getDirection()[a], (max[a] - r.getOrigin()[a]) / r.getDirection()[a]);
			t_min = fmax(t0, t_min);
			t_max = fmin(t1, t_max);
			if (t_max <= t_min)
				return false;
		}
		return true;
	}
	// inline bool hit(const Ray &r, double t_min, double t_max) const
	// {
	// 	for (int a = 0; a < 3; a++)
	// 	{
	// 		auto t0 = fmin((minimum[a] - r.origin()[a]) / r.direction()[a], (maximum[a] - r.origin()[a]) / r.direction()[a]);
	// 		auto t1 = fmax((minimum[a] - r.origin()[a]) / r.direction()[a], (maximum[a] - r.origin()[a]) / r.direction()[a]);
	// 		t_min = fmax(t0, t_min);
	// 		t_max = fmin(t1, t_max);
	// 		if (t_max <= t_min)
	// 			return false;
	// 	}

	// 	return true;
	// }
	Vector3f min, max;
};

#endif