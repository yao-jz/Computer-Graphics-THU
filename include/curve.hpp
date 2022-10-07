#ifndef CURVE_HPP
#define CURVE_HPP

#include "object3d.hpp"
#include <vecmath.h>
#include <vector>
#include <utility>

#include <algorithm>
#include <iostream>

using namespace std;

double C(int n, int m) //组合数
{
    if (m == 0)
        return 1;
    double ans = 1;
    if (m > n - m)
        m = n - m;
    for (int i = n - m + 1, j = 1; i <= n; ++i, ++j)
        ans = ans * i / j;
    return ans;
}
// TODO (PA3): Implement Bernstein class to compute spline basis function.
//       You may refer to the python-script for implementation.

// The CurvePoint object stores information about a point on a curve
// after it has been tesselated: the vertex (V) and the tangent (T)
// It is the responsiblility of functions that create these objects to fill in all the data.
struct CurvePoint
{

    Vector3f V; // Vertex
    Vector3f T; // Tangent  (unit)单位切向量
    double t;   // 对应的t值
};

class Curve : public Object3D
{
protected:
    std::vector<Vector3f> controls;

public:
    explicit Curve(std::vector<Vector3f> points) : controls(std::move(points)) {}

    bool intersect(const Ray &r, Hit &h, float tmin) override
    {
        return false;
    }

    virtual bool bounding_box(double time0, double time1, aabb& output_box)const {
        return false;
    }

    std::vector<Vector3f> &getControls()
    {
        return controls;
    }

    virtual void discretize(int resolution, std::vector<CurvePoint> &data) = 0;
    virtual double dpx(double t) = 0;
    virtual double dpy(double t) = 0;
    virtual double px(double t) = 0;
    virtual double py(double t) = 0;

    void drawGL() override
    {
        Object3D::drawGL();
        glPushAttrib(GL_ALL_ATTRIB_BITS);
        glDisable(GL_LIGHTING);
        glColor3f(1, 1, 0);
        glBegin(GL_LINE_STRIP);
        for (auto &control : controls)
        {
            glVertex3fv(control);
        }
        glEnd();
        glPointSize(4);
        glBegin(GL_POINTS);
        for (auto &control : controls)
        {
            glVertex3fv(control);
        }
        glEnd();
        std::vector<CurvePoint> sampledPoints;
        discretize(30, sampledPoints);
        glColor3f(1, 1, 1);
        glBegin(GL_LINE_STRIP);
        for (auto &cp : sampledPoints)
        {
            glVertex3fv(cp.V);
        }
        glEnd();
        glPopAttrib();
    }
};

class BezierCurve : public Curve
{
public:
    explicit BezierCurve(const std::vector<Vector3f> &points) : Curve(points)
    {
        if (points.size() < 4 || points.size() % 3 != 1)
        {
            printf("Number of control points of BezierCurve must be 3n+1!\n");
            exit(0);
        }
    }
    virtual bool bounding_box(double time0, double time1, aabb& output_box)const {
        return false;
    }

    double b(int i, double t, int n)
    {

        if (i == 0)
            return (1 - t) * (1 - t) * (1 - t);
        else if (i == 1)
            return 3 * t * (1 - t) * (1 - t);
        else if (i == 2)
            return 3 * t * t * (1 - t);
        else if (i == 3)
            return t * t * t;
    }
    double d_b(int i, double t, int n)
    {

        if (i == 0)
            return -3 * (1 - t) * (1 - t);
        else if (i == 1)
            return 9 * t * t - 12 * t + 3;
        else if (i == 2)
            return -9 * t * t + 6 * t;
        else if (i == 3)
            return 3 * t * t;
    }

    Vector3f f(double t)
    {
        Vector3f result(0, 0, 0);
        for (int i = 0; i <= n; i++)
        {
            result.x() += b(i, t, 3) * controls[i].x();
            result.y() += b(i, t, 3) * controls[i].y();
            result.z() += b(i, t, 3) * controls[i].z();
        }
        return result;
    }

    Vector3f de_f(double t)
    {
        Vector3f result(0, 0, 0);
        for (int i = 0; i <= n; i++)
        {
            result.x() += d_b(i, t, 3) * controls[i].x();
            result.y() += d_b(i, t, 3) * controls[i].y();
            result.z() += d_b(i, t, 3) * controls[i].z();
        }
        return result;
    }

    double dpx(double t) override
    {
        return controls[0].x() * d_b(0, t, 3) + controls[1].x() * d_b(1, t, 3) + controls[2].x() * d_b(2, t, 3) + controls[3].x() * d_b(3, t, 3);
    }
    double dpy(double t) override
    {
        return controls[0].y() * d_b(0, t, 3) + controls[1].y() * d_b(1, t, 3) + controls[2].y() * d_b(2, t, 3) + controls[3].y() * d_b(3, t, 3);
    }
    double px(double t) override
    {
        return controls[0].x() * b(0, t, 3) + controls[1].x() * b(1, t, 3) + controls[2].x() * b(2, t, 3) + controls[3].x() * b(3, t, 3);
    }
    double py(double t) override
    {
        return controls[0].y() * b(0, t, 3) + controls[1].y() * b(1, t, 3) + controls[2].y() * b(2, t, 3) + controls[3].y() * b(3, t, 3);
    }

    void discretize(int resolution, std::vector<CurvePoint> &data) override
    { //离散化
        data.clear();
        // TODO (PA3): fill in data vector
        n = int(controls.size() - 1);
        double dt = 1.0 / resolution;
        double now_t = 0;
        while (now_t < 1)
        {
            CurvePoint temp;
            temp.V = f(now_t);
            temp.T = de_f(now_t);
            temp.t = now_t;
            data.push_back(temp);
            now_t += dt;
        }
    }

protected:
    int n;
};

class BsplineCurve : public Curve
{
public:
    BsplineCurve(const std::vector<Vector3f> &points) : Curve(points)
    {
        if (points.size() < 4)
        {
            printf("Number of control points of BspineCurve must be more than 4!\n");
            exit(0);
        }
    }

    void init_knots(int n, int k)
    {
        knots.clear();
        for (int i = 0; i <= n + k + 1; i++)
        {
            knots.push_back(double(i) / double(n + k + 1));
        }
    }
    virtual bool bounding_box(double time0, double time1, aabb& output_box)const {
        return false;
    }

    float b(int i, int k, float t)
    {
        for (int j = 0; j < 100; j++)
        {
            for (int w = 0; w < 100; w++)
            {
                store[j][w] = 0;
            }
        }
        for (int j = 0; j < n + k + 1; j++)
        {
            if (knots[j] <= t && t < knots[j + 1])
                store[j][0] = 1;
            else
                store[j][0] = 0;
        }
        for (int p = 1; p <= k; p++)
        {
            for (int j = 0; j < n + 1; j++)
            {
                store[j][p] = ((t - knots[j]) / (knots[j + p] - knots[j])) * store[j][p - 1] + ((knots[j + p + 1] - t) / (knots[j + p + 1] - knots[j + 1])) * store[j + 1][p - 1];
            }
        }
        return store[i][k];
    }

    float d_b(int i, int k, float t)
    {
        return k * (b(i, k - 1, t) /
                        (knots[i + k] - knots[i]) -
                    b(i + 1, k - 1, t) /
                        (knots[i + k + 1] - knots[i + 1]));
    }

    Vector3f f(float t)
    {
        Vector3f result(0, 0, 0);
        for (int i = 0; i <= n; i++)
        {
            result.x() += b(i, 3, t) * controls[i].x();
            result.y() += b(i, 3, t) * controls[i].y();
            result.z() += b(i, 3, t) * controls[i].z();
        }
        return result;
    }

    Vector3f de_f(float t)
    {
        Vector3f result(0, 0, 0);
        for (int i = 0; i <= n; i++)
        {
            result.x() += d_b(i, 3, t) * controls[i].x();
            result.y() += d_b(i, 3, t) * controls[i].y();
            result.z() += d_b(i, 3, t) * controls[i].z();
        }
        return result;
    }

    void discretize(int resolution, std::vector<CurvePoint> &data) override
    {
        data.clear();
        n = int(controls.size() - 1);
        init_knots(n, 3);
        // TODO (PA3): fill in data vector
        float dt = (knots[controls.size()] - knots[3]) / float(resolution * (controls.size()));
        float now_t = knots[3];
        while (now_t < knots[controls.size()])
        {
            CurvePoint temp;
            temp.V = f(now_t);
            temp.T = de_f(now_t);
            data.push_back(temp);
            now_t += dt;
        }
    }
    double dpx(double t){};
    double dpy(double t){};
    double px(double t){};
    double py(double t){};

protected:
    vector<float> knots;
    float store[100][100];
    int n;
};

#endif // CURVE_HPP
