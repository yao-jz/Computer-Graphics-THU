#ifndef REVSURFACE_HPP
#define REVSURFACE_HPP

#include "object3d.hpp"
#include "curve.hpp"
#include <iostream>
#include <cmath>
#include <tuple>
using namespace std;
typedef std::tuple<unsigned, unsigned, unsigned> Tup3u;
class RevSurface : public Object3D
{

    Curve *pCurve;
    Vector3f pos;

public:
    RevSurface(Curve *pCurve, Material *material, Vector3f position = Vector3f::ZERO) : pCurve(pCurve), Object3D(material), pos(position)
    {
        // Check flat.
        for (const auto &cp : pCurve->getControls())
        {
            if (cp.z() != 0.0)
            {
                printf("Profile of revSurface must be flat on xy plane.\n");
                exit(0);
            }
        }
        init_surface();
    }

    inline Vector3f R(const Ray &ray, double t)
    {
        return ray.getOrigin() + t * ray.getDirection();
    }
    inline Vector3f S(double u, double theta)
    {
        return pos + Vector3f(cos(theta) * pCurve->px(u), pCurve->py(u), sin(theta) * pCurve->px(u));
    }
    inline Vector3f F(const Ray &ray, double t, double u, double theta)
    {
        return R(ray, t) - S(u, theta);
    }
    inline Matrix3f dF(const Ray &ray, double t, double u, double theta)
    {
        return Matrix3f(ray.getDirection(), Vector3f(-cos(theta) * pCurve->dpx(u), -pCurve->dpy(u), -sin(theta) * pCurve->dpx(u)), Vector3f(sin(theta) * pCurve->px(u), 0, -cos(theta) * pCurve->px(u)));
    }

    struct Surface
    {
        std::vector<Vector3f> VV; // 储存顶点
        std::vector<double> Vt;   // 储存顶点在Bezier曲线对应的t值
        std::vector<Vector3f> VN; // 储存顶点法向
        std::vector<Tup3u> VF;
        std::vector<double> angle; //角度
    } surface;
    void init_surface()
    {
        // Definition for drawable surface.
        // Surface is just a struct that contains vertices, normals, and
        // faces.  VV[i] is the position of vertex i, and VN[i] is the normal
        // of vertex i.  A face is a triple i,j,k corresponding to a triangle
        // with (vertex i, normal i), (vertex j, normal j), ...
        std::vector<CurvePoint> curvePoints;
        pCurve->discretize(50, curvePoints);
        const int steps = 20;
        for (unsigned int ci = 0; ci < curvePoints.size(); ++ci)
        {
            const CurvePoint &cp = curvePoints[ci];
            for (unsigned int i = 0; i < steps; ++i)
            {
                float t = (float)i / steps; // 弧度
                Quat4f rot;
                rot.setAxisAngle(t * 2 * 3.14159, Vector3f::UP);
                Vector3f pnew = Matrix3f::rotation(rot) * cp.V;
                Vector3f pNormal = Vector3f::cross(cp.T, -Vector3f::FORWARD);
                Vector3f nnew = Matrix3f::rotation(rot) * pNormal;
                surface.VV.push_back(pnew + pos);           // 这里实现位移
                surface.VN.push_back(nnew);                 // 储存顶点法向
                surface.Vt.push_back(cp.t);                 // 储存t
                surface.angle.push_back(t * 2 * 3.1415926); // 储存弧度值
                int i1 = (i + 1 == steps) ? 0 : i + 1;      // 是否已经反过来了
                if (ci != curvePoints.size() - 1)
                {
                    surface.VF.emplace_back((ci + 1) * steps + i, ci * steps + i1, ci * steps + i);
                    surface.VF.emplace_back((ci + 1) * steps + i, (ci + 1) * steps + i1, ci * steps + i1);
                }
            }
        }
        cout << surface.angle.size() << endl;
    }

    ~RevSurface() override
    {
        delete pCurve;
    }

    virtual bool bounding_box(double time0, double time1, aabb& output_box)const {
        float minx = 1e7, miny = 1e7, minz = 1e7, maxx = -10000, maxy=-10000, maxz=-10000;
        for(auto p : surface.VV){
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

    bool intersect(const Ray &r, Hit &h, float tmin) override
    {
        // (PA3 optional TODO): implement this for the ray-tracing routine using G-N iteration.
        bool result = false;
        float t, u, theta;
        // 获得估计值
        for (int triId = 0; triId < surface.VF.size(); ++triId)
        {
            // cout << "开始寻找 ";
            unsigned i, j, k;
            i = get<0>(surface.VF[triId]);
            j = get<1>(surface.VF[triId]);
            k = get<2>(surface.VF[triId]);
            Triangle triangle(surface.VV[i], surface.VV[j], surface.VV[k], material);
            if (triangle.intersect(r, h, tmin))
            {
                t = h.getT();
                result = true;
                Vector3f c = triangle.getCenterAxis(r); // 中心坐标
                u = c.x() * surface.Vt[i] + c.y() * surface.Vt[j] + c.z() * surface.Vt[k];
                theta = c.x() * surface.angle[i] + c.y() * surface.angle[j] + c.z() * surface.angle[k];
            }
        }
        if (result)
        {
            // 交点坐标估计值
            Vector3f x(t, u, theta);
            double lr = 0.7;
            // 牛顿迭代
            for (int kk = 0; kk < 20; kk++)
            {
                t = x.x();
                u = x.y();
                theta = x.z();
                Vector3f f = F(r, t, u, theta);
                Matrix3f df = dF(r, t, u, theta);
                x = x - (df.inverse() * f);
            }
            // 更新hit的t和normal，以及贴图使用到的u和v
            h.u = x.y();
            h.v = x.z();
            h.t = x.x();
            Vector3f v1(cos(x.z()) * pCurve->dpx(x.y()), pCurve->dpy(x.y()), sin(x.z()) * pCurve->dpx(x.y()));
            Vector3f v2(-sin(x.z()) * pCurve->px(x.y()), 0, cos(x.z()) * pCurve->px(x.y()));
            h.normal = Vector3f::cross(v1, v2).normalized();
            if (Vector3f::dot(r.getDirection(), h.normal) > 0)
                h.normal = -h.normal;
        }
        return result;
    }
};

#endif //REVSURFACE_HPP
