#ifndef PLANE_H
#define PLANE_H

#include "object3d.hpp"
#include <vecmath.h>
#include <cmath>

// TODO (PA2): Copy from PA1 finish

class Plane : public Object3D
{
public:
    Plane()
    {
    }

    Plane(const Vector3f &normal, float d, Material *m, float s = 10) : Object3D(m)
    {
        this->norm = normal;
        this->d = -d;
        float length = this->norm.length();
        this->norm = this->norm.normalized();
        this->d /= length;
        this->scale = int(s);
    }

    ~Plane() override = default;

    virtual bool bounding_box(double time0, double time1, aabb& output_box) const {
        return false;
    }

    bool intersect(const Ray &r, Hit &h, float tmin) override
    {
        float t = -1.0 * (this->d + Vector3f::dot(this->norm, r.getOrigin())) / Vector3f::dot(this->norm, r.getDirection());
        if (t > 0)
        {
            if (t >= tmin && t <= h.getT())
            {
                h.set(t, this->material, this->norm.normalized());
                if (norm.z() == 1)
                {
                    h.u = r.pointAtParameter(h.getT()).x() / scale;
                    h.v = r.pointAtParameter(h.getT()).y() / scale;
                }
                else if (norm.y() == 1)
                {
                    h.u = r.pointAtParameter(h.getT()).z() / scale;
                    h.v = r.pointAtParameter(h.getT()).x() / scale;
                }
                else if (norm.x() == 1)
                {
                    h.u = r.pointAtParameter(h.getT()).z() / scale;
                    h.v = r.pointAtParameter(h.getT()).y() / scale;
                }
                return true;
            }
            else
            {
                return false;
            }
        }
        else
        {
            return false;
        }
    }

    void drawGL() override
    {
        Object3D::drawGL();
        Vector3f xAxis = Vector3f::RIGHT;
        Vector3f yAxis = Vector3f::cross(norm, xAxis);
        xAxis = Vector3f::cross(yAxis, norm);
        const float planeSize = 10.0;
        glBegin(GL_TRIANGLES);
        glNormal3fv(norm);
        glVertex3fv(d * norm + planeSize * xAxis + planeSize * yAxis);
        glVertex3fv(d * norm - planeSize * xAxis - planeSize * yAxis);
        glVertex3fv(d * norm + planeSize * xAxis - planeSize * yAxis);
        glNormal3fv(norm);
        glVertex3fv(d * norm + planeSize * xAxis + planeSize * yAxis);
        glVertex3fv(d * norm - planeSize * xAxis + planeSize * yAxis);
        glVertex3fv(d * norm - planeSize * xAxis - planeSize * yAxis);
        glEnd();
    }

public:
    int scale = 10;
    Vector3f norm;
    float d;
};

#endif //PLANE_H
