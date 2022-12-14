#ifndef LIGHT_H
#define LIGHT_H

#include <Vector3f.h>
#include "object3d.hpp"
#include "ray.hpp"

class Light
{
public:
    Light() = default;
    double brightness;
    virtual ~Light() = default;

    virtual void getIllumination(const Vector3f &p, Vector3f &dir, Vector3f &col) const = 0;
    virtual void turnOn(int idx) const = 0;
    virtual Ray randomlyEmit() const = 0;
    virtual Vector3f getColor() const = 0;
};

class DirectionalLight : public Light
{
public:
    DirectionalLight() = delete;

    DirectionalLight(const Vector3f &d, const Vector3f &c)
    {
        direction = d.normalized();
        color = c;
    }

    ~DirectionalLight() override = default;

    ///@param p unsed in this function
    ///@param distanceToLight not well defined because it's not a point light
    void getIllumination(const Vector3f &p, Vector3f &dir, Vector3f &col) const override
    {
        // the direction to the light is the opposite of the
        // direction of the directional light source
        dir = -direction;
        col = color;
    }

    void turnOn(int idx) const override
    {
        glEnable(GL_LIGHT0 + idx);
        glLightfv(GL_LIGHT0 + idx, GL_DIFFUSE, Vector4f(color, 1.0));
        glLightfv(GL_LIGHT0 + idx, GL_SPECULAR, Vector4f(color, 1.0));
        // Last component is 0.0, indicating directional light.
        glLightfv(GL_LIGHT0 + idx, GL_POSITION, Vector4f(-direction, 0.0));
    }
    Ray randomlyEmit() const { return Ray(Vector3f::ZERO, Vector3f::ZERO); }
    Vector3f getColor() const
    {
        return color;
    }

private:
    Vector3f direction;
    Vector3f color;
};

class PointLight : public Light
{
public:
    PointLight() = delete;

    PointLight(const Vector3f &p, const Vector3f &c, double bright = 1)
    {
        position = p;
        color = c;
        brightness = bright;
    }

    ~PointLight() override = default;

    void getIllumination(const Vector3f &p, Vector3f &dir, Vector3f &col) const override
    {
        // the direction to the light is the opposite of the
        // direction of the directional light source
        dir = (position - p);
        dir = dir / dir.length();
        col = color;
    }
    Vector3f getColor() const
    {
        return color;
    }

    void turnOn(int idx) const override
    {
        glEnable(GL_LIGHT0 + idx);
        glLightfv(GL_LIGHT0 + idx, GL_DIFFUSE, Vector4f(color, 1.0));
        glLightfv(GL_LIGHT0 + idx, GL_SPECULAR, Vector4f(color, 1.0));
        glLightfv(GL_LIGHT0 + idx, GL_POSITION, Vector4f(position, 1.0));
    }

    Ray randomlyEmit() const
    {
        return Ray(position, Vector3f::randomVectorOnSphere());
    }

private:
    Vector3f position;
    Vector3f color;
};

#endif // LIGHT_H
