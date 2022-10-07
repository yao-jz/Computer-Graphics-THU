#ifndef MATERIAL_H
#define MATERIAL_H

#include <cassert>
#include <vecmath.h>

#include "ray.hpp"
#include "texture.hpp"
#include "hit.hpp"
#include <iostream>
#include <glut.h>
using namespace std;
class Material
{
public:
    explicit Material(const Vector3f &d_color, const Vector3f &s_color, Vector3f &e_color, float s = 0, int t = 2, int tt = 0, char *filename = nullptr, double d1 = 0, double d2 = 0, double d3 = 0, double d4 = 0, Vector3f co = Vector3f::ZERO) : diffuseColor(d_color), specularColor(s_color), emissionColor(e_color), shininess(s), type(t), tt(tt), refl(d1), diff(d2), spec(d3), refr(d4)
    {
        // if(tt == 1) albedo = new checker_texture(Vector3f(0.2, 0.3, 0.1), Vector3f(0.9, 0.9, 0.9));
        if (tt == 0)
        {
            // 单色
            albedo = new solid_color(co);
        }
        else if (tt == 1)
        {
            // 图片
            albedo = new image_texture(filename);
        }
    }

    virtual ~Material() = default;

    virtual Vector3f getDiffuseColor() const
    {
        return diffuseColor;
    }

    virtual Vector3f getSpecularColor() const
    {
        return specularColor;
    }

    virtual Vector3f getEmissionColor() const
    {
        return emissionColor;
    }

    Vector3f Shade(const Ray &ray, const Hit &hit,
                   const Vector3f &dirToLight, const Vector3f &lightColor)
    {
        // 材质颜色
        Vector3f this_color = diffuseColor;
        if (tt == 1)
            this_color = albedo->value(hit.u, hit.v, ray.pointAtParameter(hit.getT()));

        Vector3f shaded = Vector3f::ZERO;
        //
        Vector3f normalDirToLight = dirToLight.normalized();
        Vector3f hitNormal = hit.getNormal().normalized();
        Vector3f fanshe = 2 * (Vector3f::dot(hitNormal, normalDirToLight)) * hitNormal - normalDirToLight;
        float num1 = Vector3f::dot(normalDirToLight, hitNormal);
        if (num1 < 0)
            num1 = 0;
        float num2 = Vector3f::dot((-1 * ray.getDirection().normalized()), fanshe.normalized());
        //这里，应该是先判断是否小于零，再计算s次方
        if (num2 < 0)
            num2 = 0;
        num2 = pow(num2, this->shininess);
        for (int i = 0; i < 3; i++)
        {
            shaded[i] = lightColor[i] * (this_color[i] * num1 + specularColor[i] * num2);
        }
        return shaded;
    }

    // For OpenGL, this is fully implemented
    void Use()
    {
        glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, Vector4f(diffuseColor, 1.0f));
        glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, Vector4f(specularColor, 1.0f));
        glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, Vector2f(shininess * 4.0, 1.0f));
    }

    Vector3f getColor(double u, double v, Vector3f pos)
    {
        return albedo->value(u, v, pos);
    }

    int type; //1反射，2漫反射，3折射

    double refl; // reflection ratio
    double diff; // diffusion ratio
    double spec; // high light diffusion
    double refr; // refraction ratio

    Vector3f diffuseColor;  //漫反射
    Vector3f specularColor; //镜面反射
    Vector3f emissionColor;
    float shininess;
    // 未初始化albedo的solid color
    texture *albedo;
    texture *absorb; // for ppm
    int tt = 0;      //1有图片，0是固定颜色
};

#endif // MATERIAL_H
