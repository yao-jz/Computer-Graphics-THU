#ifndef CAMERA_H
#define CAMERA_H

#include "ray.hpp"
#include <vecmath.h>
#include <vecio.h>
#include <float.h>
#include <cmath>
#include <glut.h>

class Camera
{
public:
    double time0, time1;
    Camera(const Vector3f &center, const Vector3f &direction, const Vector3f &up, int imgW, int imgH, double _time0, double _time1)
    {
        this->center = center;
        this->direction = direction.normalized();
        this->horizontal = Vector3f::cross(this->direction, up);
        this->horizontal.normalize();
        this->up = Vector3f::cross(this->horizontal, this->direction);
        this->width = imgW;
        this->height = imgH;
        this->time0 = _time0;
        this->time1 = _time1;
    }
    Camera(const Vector3f &center, const Vector3f &direction, const Vector3f &up, int imgW, int imgH, float fdepth)
    {
        this->center = center;
        this->direction = direction.normalized();
        this->horizontal = Vector3f::cross(this->direction, up);
        this->horizontal.normalize();
        this->up = Vector3f::cross(this->horizontal, this->direction);
        this->width = imgW;
        this->height = imgH;
        this->time0 = 0;
        this->time1 = 0;
        this->fdepth = fdepth;
    }

    // Generate rays for each screen-space coordinate
    virtual Ray generateRay(const Vector2f &point) = 0;
    virtual void setupGLMatrix()
    {
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
        gluLookAt(center.x(), center.y(), center.z(),                                                 // Position
                  center.x() + direction.x(), center.y() + direction.y(), center.z() + direction.z(), // LookAt
                  up.x(), up.y(), up.z());                                                            // Up direction
    }

    virtual ~Camera() = default;

    int getWidth() const { return width; }
    int getHeight() const { return height; }

    void setCenter(const Vector3f &pos)
    {
        this->center = pos;
    }
    Vector3f getCenter() const
    {
        return this->center;
    }

    void setRotation(const Matrix3f &mat)
    {
        this->horizontal = mat.getCol(0);
        this->up = -mat.getCol(1);
        this->direction = mat.getCol(2);
    }
    Matrix3f getRotation() const
    {
        return Matrix3f(this->horizontal, -this->up, this->direction);
    }

    virtual void resize(int w, int h)
    {
        width = w;
        height = h;
    }

public:
    // Extrinsic parameters
    Vector3f center;
    Vector3f direction;
    Vector3f up;
    Vector3f horizontal;
    // Intrinsic parameters
    int width;
    int height;
    Vector3f dx, dy;
    Vector3f origin;
    float fdepth;
    virtual Ray addDepth(const Ray &ray) = 0;
};

class PerspectiveCamera : public Camera
{

    // Perspective intrinsics
    float fx; //焦距
    float fy;
    float cx;
    float cy;
    float fovyd; //度数

public:
    Ray addDepth(const Ray &ray) override
    {
        Vector3f fpoint = ray.getOrigin() + fdepth * ray.getDirection();
        return Ray(ray.getOrigin() + Vector3f::randomVectorOnSphere() * fdepth, (fpoint - ray.getOrigin()).normalized());
    }
    float getFovy() const { return fovyd; }

    PerspectiveCamera(const Vector3f &center, const Vector3f &direction,
                      const Vector3f &up, int imgW, int imgH, float angle, double _time0 = 0.0, double _time1 = 0.0) : Camera(center, direction, up, imgW, imgH, _time0, _time1)
    {
        // angle is fovy in radian.
        fovyd = angle / 3.1415 * 180.0;
        fx = fy = (float)height / (2 * tanf(angle / 2));
        cx = width / 2.0f;
        cy = height / 2.0f;
    }
    PerspectiveCamera(const Vector3f &center, const Vector3f &direction,
                      const Vector3f &up, int imgW, int imgH, float angle, float fdepth) : Camera(center, direction, up, imgW, imgH, fdepth)
    {
        // angle is fovy in radian.
        fovyd = angle / 3.1415 * 180.0;
        fx = fy = (float)height / (2 * tanf(angle / 2));
        cx = width / 2.0f;
        cy = height / 2.0f;
    }
    void resize(int w, int h) override
    {
        fx *= (float)h / height;
        fy = fx;
        Camera::resize(w, h);
        cx = width / 2.0f;
        cy = height / 2.0f;
    }
    inline double random_double() { return rand() / (RAND_MAX + 1.0); }
    inline double random_double(double min, double max) { return min + (max - min) * random_double(); }
    Ray generateRay(const Vector2f &point) override
    {
        float csx = (point.x() - cx) / fx;
        float csy = (point.y() - cy) / fy;
        Vector3f dir(csx, -csy, 1.0f);
        Matrix3f R(this->horizontal, -this->up, this->direction);
        dir = R * dir;
        dir = dir / dir.length();
        Ray ray(this->center, dir, random_double(time0, time1));
        if (fdepth != 0)
            return addDepth(ray);
        else
            return ray;
    }

    void setupGLMatrix() override
    {
        // Extrinsic.
        Camera::setupGLMatrix();
        // Perspective Intrinsic.
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        // field of view in Y, aspect ratio, near crop and far crop.
        gluPerspective(fovyd, cx / cy, 0.01, 100.0);
    }
};

#endif //CAMERA_H
