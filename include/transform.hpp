#ifndef TRANSFORM_H
#define TRANSFORM_H

#include <vecmath.h>
#include <iostream>
#include "object3d.hpp"
using namespace std;
static Vector3f transformPoint(const Matrix4f &mat, const Vector3f &point)
{
    return (mat * Vector4f(point, 1)).xyz();
}

// transform a 3D directino using a matrix, returning a direction
// This function *does not* take the inverse tranpose for you.
static Vector3f transformDirection(const Matrix4f &mat, const Vector3f &dir)
{
    return (mat * Vector4f(dir, 0)).xyz();
}

class Transform : public Object3D
{
public:
    Transform() {}

    Transform(const Matrix4f &m, Object3D *obj, Vector3f s, Vector3f trans) : o(obj),scale(s),translate(trans)
    {
        transform = m.inverse();
    }

    ~Transform()
    {
    }

    virtual bool intersect(const Ray &r, Hit &h, float tmin)
    {
        Vector3f trSource = transformPoint(transform, r.getOrigin());
        Vector3f trDirection = transformDirection(transform, r.getDirection());
        Ray tr(trSource, trDirection);
        bool inter = o->intersect(tr, h, tmin);
        if (inter)
        {
            h.set(h.getT(), h.getMaterial(), transformDirection(transform.transposed(), h.getNormal()).normalized());
        }
        return inter;
    }

    virtual bool bounding_box(double time0, double time1, aabb &output_box) const
    {
        aabb temp_box;
        o->bounding_box(time0, time1, temp_box);
        // 输出信息
        cout << "transform " << (translate+temp_box.min)*scale << " " << (translate+temp_box.max)*scale << endl;
        output_box = aabb((translate+temp_box.min)*scale, (translate+temp_box.max)*scale);
        return true;
    }

    void drawGL() override
    {
        Object3D::drawGL();
        glMatrixMode(GL_MODELVIEW);
        glPushMatrix();
        glMultMatrixf(transform.inverse());
        o->drawGL();
        glPopMatrix();
    }

protected:
    Object3D *o; //un-transformed object
    Matrix4f transform;
    Vector3f scale;
    Vector3f translate;
};

#endif //TRANSFORM_H
