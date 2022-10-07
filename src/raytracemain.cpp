#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <iostream>
#include "scene_parser.hpp"
#include "image.hpp"
#include "camera.hpp"
#include <stdio.h>
#include "group.hpp"
#include "light.hpp"
#include "photon.hpp"
#include <cmath>
#include <cstdlib>
#include <string>

using namespace std;
inline double clamp(double x) { return x < 0 ? 0 : x > 1 ? 1
                                                         : x; }
const int samples_per_pixel = 4000; //抗锯齿程度
const int MAX_PHOTON = 1000;
Group *baseGroup;
Object3D *lightSphere;
Vector3f BackgroundColor;
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
Vector3f random_in_half_unit_sphere(const Vector3f &n)
{
    while (true)
    {
        Vector3f p = random_in_unit_sphere();
        if (Vector3f::dot(p, n) <= 0)
            continue;
        return p;
    }
}

Vector3f raytrace(Ray &ray, Hit &hit, int depth, SceneParser &parser)
{
    if (++depth > 5)
        return Vector3f(0, 0, 0);
    bool isIntersect = baseGroup->intersect(ray, hit, 0);
    if (!isIntersect)
        return Vector3f(0, 0, 0);
    Vector3f color(0, 0, 0);
    color = hit.getMaterial()->getEmissionColor();
    Vector3f origin = ray.pointAtParameter(hit.getT());
    if (hit.getMaterial()->type != 3)
    {
        // for (int li = 0; li < parser.getNumLights(); li++)
        // {
        //     Light *light = parser.getLight(li);
        //     Vector3f L, lightColor;
        //     light->getIllumination(ray.pointAtParameter(hit.getT()), L, lightColor);
        //     Vector3f direction = L.normalized();
        //     Vector3f origin = ray.pointAtParameter(hit.getT());
        //     origin = origin + 0.001 * direction;
        //     origin = ray.pointAtParameter(hit.getT());
        //     color += hit.getMaterial()->Shade(ray, hit, L, lightColor);
        // }
    }
    //漫反射，需要处理texture

    Vector3f direction = hit.getNormal() + random_in_unit_sphere();
    Ray diffuse_ray(origin, direction);
    Hit h1;
    if (hit.getMaterial()->type != 3)
        color += hit.getMaterial()->getDiffuseColor() * raytrace(diffuse_ray, h1, depth, parser);
    else
        color += 0.2 * hit.getMaterial()->getDiffuseColor() * raytrace(diffuse_ray, h1, depth, parser);
    origin = ray.pointAtParameter(hit.getT());
    
    if (hit.getMaterial()->type == 1)
    {
        // 反射
        Vector3f I = ray.getDirection().normalized();
        Vector3f N = hit.getNormal().normalized();
        Vector3f reflect_direction = (I - 2.0 * (Vector3f::dot(I, N)) * N);
        origin += 0.01 * reflect_direction;
        Ray reflect_ray(origin, reflect_direction.normalized());
        Hit h1;
        color += hit.getMaterial()->getDiffuseColor() * raytrace(reflect_ray, h1, depth, parser);
    }
    origin = ray.pointAtParameter(hit.getT());
    if (hit.getMaterial()->type == 3)
    {
        //折射
        Vector3f I = ray.getDirection().normalized();
        Vector3f N = hit.getNormal().normalized();
        Vector3f reflect_direction = (I - 2.0 * (Vector3f::dot(I, N)) * N);
        origin += 0.01 * reflect_direction;
        Ray reflect_ray(origin, reflect_direction);
        Hit h2;
        Vector3f nl = (Vector3f::dot(N, I) < 0) ? N : N * -1;
        bool into = Vector3f::dot(nl, N) > 0;
        float nc = 1, nt = 1.5, nnt = (into) ? (nc / nt) : (nt / nc), ddn = Vector3f::dot(I, nl), cos2t;
        if ((cos2t = 1 - nnt * nnt * (1 - ddn * ddn)) < 0)
        {
            color += hit.getMaterial()->getDiffuseColor() * raytrace(reflect_ray, h2, depth, parser);
        }
        else
        {
            Vector3f t;
            t = nnt * I + (nnt * Vector3f::dot(-1 * I, nl) - sqrt(cos2t)) * nl;
            t = t.normalized();
            origin = ray.pointAtParameter(hit.getT());
            origin += 0.01 * t;
            Ray refract_ray(origin, t);
            Hit h3;
            color += hit.getMaterial()->getDiffuseColor() * raytrace(refract_ray, h3, depth, parser);
        }
    }
    return color;
}

int main(int argc, char *argv[])
{
    for (int argNum = 1; argNum < argc; ++argNum)
    {
        std::cout << "Argument " << argNum << " is: " << argv[argNum] << std::endl;
    }

    if (argc != 3)
    {
        cout << "Usage: ./bin/PA3 <input scene file> <output bmp/ppm file>" << endl;
        return 1;
    }
    string inputFile = argv[1];
    string outputFile = argv[2]; // only bmp is allowed.

    SceneParser parser(argv[1]);
    BackgroundColor = parser.getBackgroundColor();
    Image image(parser.getCamera()->getWidth(), parser.getCamera()->getHeight());
    auto camera = parser.getCamera();
    camera->time0 = 0.0;
    camera->time1 = 8.0;
#pragma omp parallel for schedule(dynamic, 1)
    for (int x = 0; x < camera->getWidth(); ++x)
    {
        fprintf(stderr, "\rRendering %5.2f%%", samples_per_pixel, 100. * x / (camera->getWidth() - 1));
        for (int y = 0; y < camera->getHeight(); ++y)
        {
            Vector3f color = Vector3f(0, 0, 0);
            for (int j = 0; j < samples_per_pixel; j++)
            {
                // Ray camRay = parser.getCamera()->generateRay(Vector2f(x, y));
                Ray camRay = parser.getCamera()->generateRay(Vector2f(x + random_double(), y + random_double()));
                baseGroup = parser.getGroup();
                Hit hit;
                bool isIntersect = baseGroup->intersect(camRay, hit, 0);
                Hit h;
                if (isIntersect)
                {
                    color += raytrace(camRay, h, 0, parser);
                }
                else
                {
                    color += Vector3f(0, 0, 0);
                }
            }
            color = color / samples_per_pixel;
            // color = Vector3f(sqrt(color[0]), sqrt(color[1]), sqrt(color[2]));
            image.SetPixel(x, y, color);
        }
    }
    image.SaveImage(argv[2]);
    cout << "Hello! Computer Graphics!" << endl;
    return 0;
}
