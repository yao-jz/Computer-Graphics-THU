#ifndef MESH_H
#define MESH_H

#include <vector>
#include "object3d.hpp"
#include "triangle.hpp"
#include "Vector2f.h"
#include "Vector3f.h"
#include "group.hpp"
#include <iostream>
using namespace std;
class Mesh : public Object3D
{

public:
    Mesh(const char *filename, Material *m);

    struct TriangleIndex
    {
        TriangleIndex()
        {
            x[0] = 0;
            x[1] = 0;
            x[2] = 0;
        }
        int &operator[](const int i) { return x[i]; }
        // By Computer Graphics convention, counterclockwise winding is front face
        int x[3]{}; //逆时针缠绕是正面
    };

    bvh_node *bvh_root = nullptr;
    vector<Object3D *> triangles;
    std::vector<Vector3f> v;
    std::vector<TriangleIndex> t;
    std::vector<Vector3f> n;
    bool intersect(const Ray &r, Hit &h, float tmin) override;
    virtual bool bounding_box(double time0, double time1, aabb &output_box) const
    {
        float minx = 1e7, miny = 1e7, minz = 1e7, maxx = -10000, maxy = -10000, maxz = -10000;
        for (auto p : v)
        {
            minx = min(minx, p.x());
            miny = min(miny, p.y());
            minz = min(minz, p.z());
            maxx = max(maxx, p.x());
            maxy = max(maxy, p.y());
            maxz = max(maxz, p.z());
        }
        output_box = aabb(Vector3f(minx, miny, minz), Vector3f(maxx, maxy, maxz));
        return true;
    }

    void init()
    {
        if(bvh_root) return;
        // 初始化mesh节点
        for (int triId = 0; triId < (int)t.size(); ++triId)
        {
            auto triIndex = t[triId];
            Triangle* ttt = new Triangle(v[triIndex[0]],
                                        v[triIndex[1]], v[triIndex[2]], material);
            ttt->normal = n[triId];
            triangles.push_back(ttt);
        }
        // cout << triangles.size() << endl;
        if (!bvh_root)
            bvh_root = new bvh_node(triangles, 0, triangles.size(), 0, 0);
    }
    void drawGL() override
    {
        // TODO (PA2): Call drawGL for each individual triangle.
        for (int triId = 0; triId < (int)t.size(); ++triId)
        {
            Object3D::drawGL();
            TriangleIndex &triIndex = t[triId];
            glBegin(GL_TRIANGLES);
            glNormal3fv(n[triId]);
            glVertex3fv(v[triIndex[0]]);
            glVertex3fv(v[triIndex[1]]);
            glVertex3fv(v[triIndex[2]]);
            glEnd();
        }
    }

private:
    // Normal can be used for light estimation
    void computeNormal();
};

#endif
