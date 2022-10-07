#ifndef OBJECT3D_H
#define OBJECT3D_H

#include "ray.hpp"
#include "hit.hpp"
#include "material.hpp"
#include <glut.h>
#include "aabb.hpp"

// Base class for all 3d entities.
class Object3D
{
public:
    Object3D() : material(nullptr) {}

    virtual ~Object3D() = default;

    explicit Object3D(Material *material)
    {
        this->material = material;
    }

    // PA1: Intersect Ray with this object. If hit, store information in hit structure.
    // This will not be used in PA2.
    virtual bool intersect(const Ray &r, Hit &h, float tmin) = 0;
    virtual bool bounding_box(double time0, double time1, aabb &output_box) const = 0;
    aabb surrounding_box(aabb box0, aabb box1) const
    {
        Vector3f small(fmin(box0.min.x(), box1.min.x()), fmin(box0.min.y(), box1.min.y()),
                       fmin(box0.min.z(), box1.min.z()));
        Vector3f big(fmax(box0.max.x(), box1.max.x()), fmax(box0.max.y(), box1.max.y()),
                     fmax(box0.max.z(), box1.max.z()));
        return aabb(small, big);
    }
    // virtual void generatePhoton(Vector3f &origin, Vector3f &direction, double &powerscale, Hit& h);
    // PA2: draw using OpenGL pipeline.
    virtual void
    drawGL()
    {
        if (material)
            material->Use();
    }
    double u;
    double v;

protected:
    Material *material;
};

#endif
