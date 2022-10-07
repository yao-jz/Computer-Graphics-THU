#ifndef SPHERE_H
#define SPHERE_H

#include "object3d.hpp"
#include <vecmath.h>
#include <cmath>
#include <glut.h>


// TODO (PA2): Copy from PA1  finish

class Sphere : public Object3D
{
public:
    Sphere()
    {
        // unit ball at the center
    }

    Sphere(const Vector3f &center, float radius, Material *material) : Object3D(material)
    {
        //
        this->center = center;
        this->radius = radius;
    }

    ~Sphere() override = default;

    inline double random_double() { return rand() / (RAND_MAX + 1.0); }
    inline double random_double(double min, double max) { return min + (max - min) * random_double(); }
    inline Vector3f random_vec() { return Vector3f(random_double(), random_double(), random_double()); }
    inline Vector3f random_vec(double min, double max) { return Vector3f(random_double(min, max), random_double(min, max), random_double(min, max)); }
    Vector3f random_in_unit_sphere()
    { //生成单位球内的向量
        while (true)
        {
            Vector3f p = random_vec(-1, 1);
            if (p.length() >= 1)
                continue;
            return p;
        }
    }

    // void generatePhoton(Vector3f &origin, Vector3f &direction, double &powerscale, Hit &h)
    // {
    //     origin = center;
    //     direction = random_in_unit_sphere();
    //     powerscale = Vector3f::dot(direction, h.getNormal()); //不同角度光强不同
    // }

    virtual bool bounding_box(double time0, double time1, aabb& output_box) const {
        output_box = aabb(center - Vector3f(radius, radius, radius), center + Vector3f(radius, radius, radius));
        return true;
    }

    bool intersect(const Ray &r, Hit &h, float tmin) override
    {
        float t = 1e38;
        //计算t
        bool position = true; //true: 在球体外部，false: 在球体内部
        Vector3f l = this->center - r.getOrigin();
        if (l.squaredLength() <= this->radius * this->radius)
            position = false;
        else
            position = true;
        float tp = Vector3f::dot(l, r.getDirection().normalized());
        if (tp < 0 && position)
            return false;
        float d = sqrt(l.squaredLength() - tp * tp);
        if (d > this->radius)
            return false;
        float tt = sqrt(this->radius * this->radius - d * d);
        //检查t
        if (position) //球体外部
            t = tp - tt;
        else //球体内部
            t = tp + tt;
        //更新t
        if (t <= h.getT() && t >= tmin)
        {
            //计算法向量
            Vector3f jiaodian = r.getOrigin() + t * r.getDirection().normalized();
            Vector3f newN = (jiaodian - this->center).normalized();
            h.set(t, this->material, newN);
            Sphere::get_sphere_uv(newN, u, v);
            h.u = u;
            h.v = v;
            return true;
        }
        else if (t > h.getT() && t >= tmin)
        {
            return true;
        }
        else
            return false;
    }

    void drawGL() override
    {
        Object3D::drawGL();
        glMatrixMode(GL_MODELVIEW);
        glPushMatrix();
        glTranslatef(center.x(), center.y(), center.z());
        glutSolidSphere(radius, 80, 80);
        glPopMatrix();
    }

    inline void roll_angle(float &angle, double a, double max)
    {
        if (angle > max - a)
            angle = angle + a - max;
        else
            angle += a;
    }
    void get_sphere_uv(const Vector3f &p, double &u, double &v)
    {
        float pi = 3.14159265;
        auto theta = acos(-p.y()); // [0, pi]
        auto phi = atan2(-p.z(), p.x()) + pi;
        //下面的东西都是图像在正中央的情况
        roll_angle(theta, -pi / 14, pi);      //减小时向下
        roll_angle(phi, 5 * pi / 14, 2 * pi); // 减小是向右
        u = phi / (2 * pi);
        v = theta / pi;
    }

protected:
    Vector3f center;
    float radius;
};

class moving_sphere : public Sphere
{
public:
    moving_sphere() {}
    moving_sphere(Vector3f c0, Vector3f c1, double _time0, double _time1, double r, Material *mat)
    {
        radius = r;
        material = mat;
        center0 = c0;
        center1 = c1;
        time0 = _time0;
        time1 = _time1;
    }
    virtual bool bounding_box(double time0, double time1, aabb& output_box) const {
        aabb box0(center(time0) - Vector3f(radius, radius, radius), center(time0)+Vector3f(radius, radius, radius));
        aabb box1(center(time1) - Vector3f(radius, radius, radius), center(time1)+Vector3f(radius, radius, radius));
        output_box = surrounding_box(box0, box1);
        return true;
    }

    bool intersect(const Ray &r, Hit &h, float tmin) override
    {
        float t = 1e38;
        //计算t
        bool position = true; //true: 在球体外部，false: 在球体内部
        Vector3f l = center(r.time) - r.getOrigin();
        if (l.squaredLength() <= this->radius * this->radius)
            position = false;
        else
            position = true;
        float tp = Vector3f::dot(l, r.getDirection().normalized());
        if (tp < 0 && position)
            return false;
        float d = sqrt(l.squaredLength() - tp * tp);
        if (d > this->radius)
            return false;
        float tt = sqrt(this->radius * this->radius - d * d);
        //检查t
        if (position) //球体外部
            t = tp - tt;
        else //球体内部
            t = tp + tt;
        //更新t
        if (t <= h.getT() && t >= tmin)
        {
            //计算法向量
            Vector3f jiaodian = r.getOrigin() + t * r.getDirection().normalized();
            Vector3f newN = (jiaodian - center(r.time)).normalized();
            h.set(t, this->material, newN);
            Sphere::get_sphere_uv(newN, u, v);
            h.u = u;
            h.v = v;
            return true;
        }
        else if (t > h.getT() && t >= tmin)
        {
            return true;
        }
        else
            return false;
    }

    Vector3f center(double time) const
    {
        return center0 + ((time - time0) / (time1 - time0)) * (center1 - center0);
    }

    Vector3f center0, center1;
    double time0, time1;
};

#endif
