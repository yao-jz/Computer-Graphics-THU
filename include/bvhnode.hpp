#ifndef BVH_HPP
#define BVH_HPP

#include "object3d.hpp"
#include "ray.hpp"
#include "group.hpp"
#include <algorithm>
#include "hit.hpp"
#include <iostream>
#include <vector>

using namespace std;
class bvh_node : public Object3D
{
public:
	bvh_node() {}
	// bvh_node(const Group &g, double time0, double time1) : bvh_node(g.group, 0, g.group.size(), time0, time1) {}
	bvh_node(const vector<Object3D *> &src_objects, size_t start, size_t end, double time0, double time1)
	{
		auto objects = src_objects;
		// cout << start << " " << end << endl;
		const int axis = random_int(0, 2);
		size_t num = end - start;
		if (num == 1)
		{
			left = right = objects[start];
		}
		else if (num == 2)
		{
			if (box_compare(objects[start], objects[start + 1], axis))
			{
				left = objects[start];
				right = objects[start + 1];
			}
			else
			{
				left = objects[start + 1];
				right = objects[start];
			}
		}
		else
		{
			if (axis == 0)
				sort(objects.begin() + start, objects.begin() + end, comparer<0>());
			else if (axis == 1)
				sort(objects.begin() + start, objects.begin() + end, comparer<1>());
			else if (axis == 2)
				sort(objects.begin() + start, objects.begin() + end, comparer<2>());
			auto mid = start + num / 2;
			left = new bvh_node(objects, start, mid, time0, time1);
			right = new bvh_node(objects, mid, end, time0, time1);
		}
		aabb box_left, box_right;
		if (!left->bounding_box(time0, time1, box_left) || !right->bounding_box(time0, time1, box_right))
			std::cerr << "no bounding box in constructor\n ";
		box = surrounding_box(box_left, box_right);
	}
	bool intersect(const Ray &r, Hit &h, float tmin) override
	{
		if (!box.hit(r))
			return false;
		bool hit_left = left->intersect(r, h, 0);
		bool hit_right = right->intersect(r, h, 0);
		return hit_left || hit_right;
	}
	virtual bool bounding_box(double time0, double time1, aabb &output_box) const override
	{
		output_box = box;
		return true;
	}
	Object3D *left;
	Object3D *right;
	aabb box;
	template <int dim>
	struct comparer
	{
		bool operator()(const Object3D *a, const Object3D *b)
		{
			aabb box_a;
			aabb box_b;
			if (!a->bounding_box(0, 0, box_a) || !b->bounding_box(0, 0, box_b))
				std::cerr << "No bounding box in bvh_node constructor.\n";
			return box_a.min[dim] < box_b.min[dim];
		}
	};
	inline bool box_compare(const Object3D *a, const Object3D *b, int axis)
	{
		aabb box_a;
		aabb box_b;
		if (!a->bounding_box(0, 0, box_a) || !b->bounding_box(0, 0, box_b))
			std::cerr << "No bounding box in bvh_node constructor.\n";
		return box_a.min[axis] < box_b.min[axis];
	}
	inline double random_double() { return rand() / (RAND_MAX + 1.0); }
	inline double random_double(double min, double max) { return min + (max - min) * random_double(); }
	inline int random_int(int min, int max) { return static_cast<int>(random_double(min, max + 1)); }
};

#endif