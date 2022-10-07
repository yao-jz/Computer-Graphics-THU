#include "ppm.hpp"
#include <string>
#include <vecmath.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include "hit.hpp"
#include "material.hpp"
#include "ray.hpp"
#include <cstdlib>
#include <omp.h>
using namespace std;

inline Vector3f get_max(Vector3f &v1, Vector3f &v2) { return Vector3f(max(v1.x(), v2.x()), max(v1.y(), v2.y()), max(v1.z(), v2.z())); }
inline Vector3f get_min(Vector3f &v1, Vector3f &v2) { return Vector3f(min(v1.x(), v2.x()), min(v1.y(), v2.y()), min(v1.z(), v2.z())); }

ViewPoint::ViewPoint(Vector3f pos, Vector3f N, Vector3f color, double s, int x, int y) : c(pos), n(N), color(color), strength(s), x(x), y(y) {}
PPM::PPM() { root = nullptr; }
PPM::~PPM() {}

void PPM::run()
{
	bg_pic = new Vector3f[x * y];
	int *samples_count = new int[x * y];
	for (int i = 0; i < x * y; i++)
		samples_count[i] = 0;

	double current_r = initial_r;
	double energy = 1.0 / log(total_round);
	double current_e = energy;
	for (int iter = 0; iter < total_round; iter++)
	{
		for (int i = 0; i < x * y; i++)
			bg_pic[i] = Vector3f(0, 0, 0);
		view_points.clear();
		srand(time(0));
#pragma omp parallel for num_threads(12)
		for (int i = start_row; i < x; i++)
		{
			for (int j = start_col; j < y; j++)
			{
				for (int ww = 0; ww < 8; ww++)
				{
					double _i = i + (rand() * 1.0 / RAND_MAX - 0.5) * 1;
					double _j = j + (rand() * 1.0 / RAND_MAX - 0.5) * 1;
					// Ray camRay = camera->generateRay(Vector2f(i, j));
					Ray camRay = camera->generateRay(Vector2f(_i, _j));
					RayTracing(camRay, Vector3f(1, 1, 1), i, j, 0, current_e, 0);
					samples_count[i * y + j]++;
				}
			}
		}
		buildKDTree(root, view_points);
		cout << "kdtree build over" << endl;
		double brightness = 0;
		// 增加每个光源的brightness
		for (int li = 0; li < parser->getNumLights(); li++)
		{
			Light *light = parser->getLight(li);
			brightness += light->brightness;
		}
		for (int li = 0; li < parser->getNumLights(); li++)
		{
			Light *light = parser->getLight(li);
			int light_emit_count = photon_num * light->brightness / brightness;
#pragma omp parallel for num_threads(12)
			for (int i = 0; i < light_emit_count; i++)
			{
				fprintf(stderr, "\rRendering %d of %d iteration %.2f% \%", iter + 1, total_round, 100 * double(i) / double(light_emit_count));
				Ray ray = light->randomlyEmit();
				photonTracing(ray, light->getColor(), 0, current_r, total_brightness / photon_num, 0, 0);
			}
		}
		rmKDTree(root);
		for (int i = 0; i < x; i++)
		{
			for (int j = 0; j < y; j++)
			{
				board[i * y + j] += bg_pic[i * y + j] * (1.0 / samples_count[i * y + j]);
			}
		}
		current_r *= round_decay; // 迭代半径
		current_e /= round_decay;
	}
}

void PPM::photonTracing(const Ray &ray, const Vector3f &rayC, int depth, double r, double lambda, double dist, int diff_count)
{
	// cout << "enter photontracing" << endl;
	if (depth > max_jump)
		return;
	Hit hit;
	bool isIntersect = baseGroup->intersect(ray, hit, 0);
	if (!isIntersect)
		return;
	dist += hit.getT();
	if (hit.getMaterial()->refl > 0.000001)
	{
		//处理反射
		// cout << "handle refl" << endl;
		Vector3f I = ray.getDirection().normalized();
		Vector3f N = hit.getNormal().normalized();
		Vector3f reflect_direction = (I - 2.0 * (Vector3f::dot(I, N)) * N);
		Vector3f origin = ray.pointAtParameter(hit.getT());
		origin += 0.001 * N;
		Ray reflect_ray(origin, reflect_direction.normalized());
		photonTracing(reflect_ray, rayC, depth + 1, r, lambda * hit.getMaterial()->refl, dist, diff_count);
		// cout << "refl end" << endl;
	}
	if (hit.getMaterial()->refr > 0.000001)
	{
		//处理折射
		// cout << " handle refr" << endl;
		Vector3f I = ray.getDirection().normalized();
		Vector3f N = hit.getNormal().normalized();
		Vector3f reflect_direction = (I - 2.0 * (Vector3f::dot(I, N)) * N);
		Vector3f origin = ray.pointAtParameter(hit.getT());
		origin += 0.001 * N;
		Ray reflect_ray(origin, reflect_direction);
		Vector3f nl = (Vector3f::dot(N, I) < 0) ? N : N * -1;
		bool into = Vector3f::dot(nl, N) > 0;
		float nc = 1, nt = 1.5, nnt = (into) ? (nc / nt) : (nt / nc), ddn = Vector3f::dot(I, nl), cos2t;
		if ((cos2t = 1 - nnt * nnt * (1 - ddn * ddn)) < 0)
		{
			photonTracing(reflect_ray, rayC, depth + 1, r, lambda * hit.getMaterial()->refr, dist, diff_count);
		}
		else
		{
			Vector3f t;
			t = nnt * I + (nnt * Vector3f::dot(-1 * I, nl) - sqrt(cos2t)) * nl;
			t = t.normalized();
			origin = ray.pointAtParameter(hit.getT());
			origin += 0.01 * t;
			Ray refract_ray(origin, t);
			photonTracing(refract_ray, rayC, depth + 1, r, lambda * hit.getMaterial()->refr, dist, diff_count);
		}
		// cout << "refr end" << endl;
	}
	if (hit.getMaterial()->diff > 0.00001)
	{
		//处理漫反射
		// cout << "handle diff" << endl;

		vector<const ViewPoint *> vps;
		Vector3f p = ray.pointAtParameter(hit.getT()) + hit.getNormal() * 0.00002;
		queryKDTree(root, vps, p, r);
		for (auto &w : vps)
		{
			if ((Vector3f::dot(w->n, ray.getDirection()) < -0.0001))
			{
				Vector3f res = rayC * w->color * pow(((r - (w->c - p).length()) / r), 2) * lambda * (1.0 / (diff_count + 1)) * w->strength;
#pragma omp critical
				//res 就是color
				bg_pic[w->x * y + w->y] += res;
			}
		}
		// 处理漫反射
		Vector3f DD = (hit.getNormal() * Vector3f::randomVectorOnSphere()).normalized();
		double phi = acos(rand() * 1.0 / RAND_MAX);
		p = ray.pointAtParameter(hit.getT()) + hit.getNormal() * 0.0001;
		Ray diffuse_ray(p, DD * sin(phi) + hit.getNormal() * cos(phi));
		photonTracing(diffuse_ray, rayC * hit.getMaterial()->getColor(hit.u, hit.v, ray.pointAtParameter(hit.getT())), depth + 1, r, lambda * hit.getMaterial()->diff, dist + hit.getT(), diff_count + 1);
		// cout << "diff end"<<endl;
	}
}

void PPM::RayTracing(const Ray &ray, Vector3f rayC, int xx, int yy, int depth, double lambda, double dist)
{
	if (lambda < 1e-9 || depth > max_depth)
		//如果没有能量了或者迭代次数超过上限
		return;
	Hit hit;
	bool isIntersect = baseGroup->intersect(ray, hit, 0);
	if (isIntersect)
	{
		if (hit.getMaterial()->emissionColor.x() != 0)
		{
			// 光源，直接光照
			bg_pic[xx * y + yy] += hit.getMaterial()->emissionColor * lambda;
		}
		else
		{
			if (hit.getMaterial()->refr > 0.000001)
			{
				//折射
				// cout << "enter refr" << endl;
				Vector3f I = ray.getDirection().normalized();
				Vector3f N = hit.getNormal().normalized();
				Vector3f reflect_direction = (I - 2.0 * (Vector3f::dot(I, N)) * N);
				Vector3f origin = ray.pointAtParameter(hit.getT());
				origin += 0.001 * N;
				Ray reflect_ray(origin, reflect_direction);
				Vector3f nl = (Vector3f::dot(N, I) < 0) ? N : N * -1;
				bool into = Vector3f::dot(nl, N) > 0;
				float nc = 1, nt = 1.5, nnt = (into) ? (nc / nt) : (nt / nc), ddn = Vector3f::dot(I, nl), cos2t;
				if ((cos2t = 1 - nnt * nnt * (1 - ddn * ddn)) < 0)
				{
					RayTracing(reflect_ray, rayC, xx, yy, depth + 1, lambda * hit.getMaterial()->refr, dist + hit.getT());
				}
				else
				{
					Vector3f t;
					t = nnt * I + (nnt * Vector3f::dot(-1 * I, nl) - sqrt(cos2t)) * nl;
					t = t.normalized();
					origin = ray.pointAtParameter(hit.getT());
					origin += 0.01 * t;
					Ray refract_ray(origin, t);
					RayTracing(refract_ray, rayC, xx, yy, depth + 1, lambda * hit.getMaterial()->refr, dist + hit.getT());
				}
				// cout << "out refr " << endl;
			}
			if (hit.getMaterial()->refl > 0.000001)
			{
				//反射
				// cout << "Enter refl" << endl;
				Vector3f I = ray.getDirection().normalized();
				Vector3f N = hit.getNormal().normalized();
				Vector3f reflect_direction = (I - 2.0 * (Vector3f::dot(I, N)) * N);
				Vector3f origin = ray.pointAtParameter(hit.getT());
				origin += 0.001 * N;
				Ray reflect_ray(origin, reflect_direction.normalized());
				RayTracing(reflect_ray, rayC, xx, yy, depth + 1, lambda * hit.getMaterial()->refl, dist + hit.getT());
				// cout << "out refl" << endl;
			}
			if (hit.getMaterial()->diff > 0.000001)
			{
				//漫反射
				// cout << "enter diff" << endl;
#pragma omp critical
				view_points.push_back(
					ViewPoint(
						ray.pointAtParameter(hit.getT()) + hit.getNormal() * 0.00002,
						hit.getNormal(),
						rayC * hit.getMaterial()->getColor(hit.u, hit.v, ray.pointAtParameter(hit.getT())),
						lambda * hit.getMaterial()->diff,
						xx,
						yy));
				// cout << "out diff" << endl;
			}
		}
	}
	else
	{
		bg_pic[xx * y + yy] += bg_color * lambda;
	}
}

void PPM::buildKDTree(KDTreeNode *&node, vector<ViewPoint> &lst, int l, int r, int dim)
{
	// lst储存这个节点内保存的所有点
	if (l == -1 && r == -1)
		l = 0, r = lst.size();
	if (l >= r)
		return;
	int mid = (l + r) >> 1;
	// 根据维度进行排序
	if (dim == 0)
		nth_element(lst.begin() + l, lst.begin() + mid, lst.begin() + r, comparer<0>());
	else if (dim == 1)
		nth_element(lst.begin() + l, lst.begin() + mid, lst.begin() + r, comparer<1>());
	else if (dim == 2)
		nth_element(lst.begin() + l, lst.begin() + mid, lst.begin() + r, comparer<2>());

	//初始化当前节点
	node = new KDTreeNode();
	node->value = lst[mid];
	node->l = node->r = nullptr;
	node->split_dimension = dim;
	node->bd_max = node->value.c;
	node->bd_min = node->value.c;

	// 建树并更新区域
	buildKDTree(node->l, lst, l, mid, (dim + 1) % 3);
	if (node->l)
	{
		node->bd_max = get_max(node->bd_max, node->l->bd_max);
		node->bd_min = get_min(node->bd_min, node->l->bd_min);
	}
	buildKDTree(node->r, lst, mid + 1, r, (dim + 1) % 3);
	if (node->r)
	{
		node->bd_max = get_max(node->bd_max, node->r->bd_max);
		node->bd_min = get_min(node->bd_min, node->r->bd_min);
	}
}

void PPM::queryKDTree(KDTreeNode *node, vector<const ViewPoint *> &result, const Vector3f &pos, double r)
{
	// 在给定半径的范围内进行搜索
	double x, y, z;

	if (pos.x() <= node->bd_max.x() && pos.x() >= node->bd_min.x())
		x = 0;
	else
		x = min(abs(pos.x() - node->bd_max.x()), abs(pos.x() - node->bd_min.x()));
	if (pos.y() <= node->bd_max.y() && pos.y() >= node->bd_min.y())
		y = 0;
	else
		y = min(abs(pos.y() - node->bd_max.y()), abs(pos.y() - node->bd_min.y()));
	if (pos.z() <= node->bd_max.z() && pos.z() >= node->bd_min.z())
		z = 0;
	else
		z = min(abs(pos.z() - node->bd_max.z()), abs(pos.z() - node->bd_min.z()));

	// 外面，停止递归
	if (x * x + y * y + z * z > r * r)
		return;

	// 判断节点内viewpoint是否在r内
	if ((node->value.c - pos).length() <= r)
		result.push_back(&(node->value));
	// 递归查询左子树右子树
	if (node->l)
		queryKDTree(node->l, result, pos, r);
	if (node->r)
		queryKDTree(node->r, result, pos, r);
}

void PPM::rmKDTree(KDTreeNode *&node)
{
	if (!node)
		return;
	rmKDTree(node->l);
	rmKDTree(node->r);
	delete node;
}