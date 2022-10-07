#ifndef PPM_HPP
#define PPM_HPP

#include <vecmath.h>
#include <iostream>
#include <vector>
#include "camera.hpp"
#include "light.hpp"
#include "group.hpp"
#include "scene_parser.hpp"

using namespace std;

class ViewPoint
{
	friend class PPM;

public:
	Vector3f c;
	Vector3f n;
	Vector3f color;
	double strength;
	int x, y;
	ViewPoint() {}
	ViewPoint(Vector3f pos, Vector3f N, Vector3f color, double s, int x, int y);
};

class PPM
{
public:
	struct KDTreeNode
	{
		ViewPoint value;
		KDTreeNode *l;
		KDTreeNode *r;
		int split_dimension;
		Vector3f bd_max, bd_min;
	};
	template <int dim>
	struct comparer
	{
		bool operator()(const ViewPoint &p1, const ViewPoint &p2)
		{
			if (dim == 0)
				return p1.c.x() < p2.c.x();
			else if (dim == 1)
				return p1.c.y() < p2.c.y();
			else if (dim == 2)
				return p1.c.z() < p2.c.z();
		}
	};
	KDTreeNode *root;

	Camera *camera;
	int x, y;
	int max_depth;
	int start_row, start_col;
	int total_round;
	int photon_num;
	double total_brightness;
	double round_decay; //0.7 or 2/3?
	double initial_r;
	int max_jump;
	vector<ViewPoint> view_points;
	Vector3f* bg_pic;
	Vector3f* board; // 最终绘制的颜色数组
	Vector3f bg_color;
	SceneParser* parser;
	Group* baseGroup;

	PPM();
	~PPM();
	void run();
	void RayTracing(const Ray &ray, Vector3f rayC, int xx, int yy, int depth, double lambda, double dist);
	void photonTracing(const Ray &ray, const Vector3f &rayC, int depth, double r, double lambda, double dist, int diff_count);

	void buildKDTree(KDTreeNode *&now, vector<ViewPoint> &lst, int l = -1, int r = -1, int dim = 0);
	void rmKDTree(KDTreeNode *&node);
	void queryKDTree(KDTreeNode *node, vector<const ViewPoint *> &result, const Vector3f &pos, double r);
};

#endif