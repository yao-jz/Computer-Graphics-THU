#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <iostream>
#include <ctime>
#include "scene_parser.hpp"
#include "image.hpp"
#include "camera.hpp"
#include <stdio.h>
#include "ppm.hpp"
#include "group.hpp"
#include "light.hpp"
#include "photon.hpp"
#include <cmath>
#include <cstdlib>
#include <string>

using namespace std;
inline double clamp(double x) { return x < 0 ? 0 : x > 1 ? 1
                                                         : x; }
Group *baseGroup;
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

int main(int argc, char *argv[])
{
    cout << "programme start" << endl;
    clock_t a,b;
    a = clock();
    string inputFile = argv[1];
    string outputFile = argv[2]; // only bmp is allowed.
    SceneParser parser(argv[1]);
    cout << "parse over " << endl;
    BackgroundColor = parser.getBackgroundColor();
    Image image(parser.getCamera()->getWidth(), parser.getCamera()->getHeight());
    auto camera = parser.getCamera();
    camera->time0 = 0.0;
    camera->time1 = 8.0;
    cout << "ppm start" << endl;
    PPM *ppm = new PPM(); ppm->camera = camera;
    ppm->x = camera->getWidth(); ppm->y = camera->getHeight();
    ppm->max_depth = 5; ppm->max_jump = 5;
    ppm->start_col = 0; ppm->start_row = 0;
    // ppm->total_round = 5; ppm->photon_num = 1000000;
    // ppm->total_brightness = 2500; ppm->round_decay = 0.95;
    ppm->total_round = 5; ppm->photon_num = 100000;
    ppm->total_brightness = 1800; ppm->round_decay = 0.95;
    ppm->initial_r = 5; ppm->bg_color = Vector3f(0,0,0);
    ppm->parser = &parser; ppm->baseGroup = parser.getGroup();
    ppm->board = new Vector3f[ppm->x*ppm->y];
    cout << ppm->total_brightness << " " << ppm->total_round << " " << ppm->photon_num << endl;
    for(int i = 0; i < ppm->x*ppm->y; i++)
        ppm->board[i] = Vector3f(0,0,0);
    ppm->run();   
    for(int i = 0; i < ppm->x; i++) {
        for(int j = 0; j < ppm->y; j++){
            Vector3f& c = (ppm->board)[i*ppm->y+j];
            // cout << c << endl;
            image.SetPixel(i, j, c);
        }
    }
    image.SaveImage(argv[2]);
    b = clock();
    cout << endl;
    cout << double(b-a) / CLOCKS_PER_SEC << endl;
    return 0;
}
