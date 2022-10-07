#ifndef TRIANGLE_H
#define TRIANGLE_H

#include "object3d.hpp"
#include <vecmath.h>
#include <cmath>
#include <iostream>

using namespace std;

// TODO (PA2): Copy from PA1   finish
class Triangle: public Object3D
{

public:
	Triangle() = delete;
        ///@param a b c are three vertex positions of the triangle

	Triangle( const Vector3f& a, const Vector3f& b, const Vector3f& c, Material* m) : Object3D(m) {
	    vertices[0] = a;
		vertices[1] = b;
		vertices[2] = c;
		Vector3f E1= vertices[0] - vertices[1];
		Vector3f E2 = vertices[0] - vertices[2];
		normal = Vector3f::cross(E1, E2).normalized();
	}

	Vector3f getCenterAxis(const Ray& ray) {
		Vector3f E1= vertices[0] - vertices[1];
		Vector3f E2 = vertices[0] - vertices[2];
		Vector3f S = vertices[0] - ray.getOrigin();
		Vector3f Rd = ray.getDirection();
		float t = Matrix3f(S, E1, E2).determinant() / Matrix3f(Rd, E1, E2).determinant();
		float beta = Matrix3f(Rd, S, E2).determinant() / Matrix3f(Rd, E1, E2).determinant();
		float gamma = Matrix3f(Rd, E1, S).determinant() / Matrix3f(Rd, E1, E2).determinant();
		return Vector3f(1-beta-gamma, beta, gamma);//对应0，1，2的顺序
	}

	virtual bool bounding_box(double time0, double time1, aabb& output_box) const {
		float minx = min(vertices[0].x(), min(vertices[1].x(), vertices[2].x()));
		float miny = min(vertices[0].y(), min(vertices[1].y(), vertices[2].y()));
		float minz = min(vertices[0].z(), min(vertices[1].z(), vertices[2].z()));
		float maxx = max(vertices[0].x(), max(vertices[1].x(), vertices[2].x()));
		float maxy = max(vertices[0].y(), max(vertices[1].y(), vertices[2].y()));
		float maxz = max(vertices[0].z(), max(vertices[1].z(), vertices[2].z()));
		output_box = aabb(Vector3f(minx, miny, minz), Vector3f(maxx, maxy, maxz));
		return true;
	}

	bool intersect( const Ray& ray,  Hit& hit , float tmin) override {
        Vector3f E1= vertices[0] - vertices[1];
		Vector3f E2 = vertices[0] - vertices[2];
		Vector3f S = vertices[0] - ray.getOrigin();
		Vector3f Rd = ray.getDirection();
		float t = Matrix3f(S, E1, E2).determinant() / Matrix3f(Rd, E1, E2).determinant();
		float beta = Matrix3f(Rd, S, E2).determinant() / Matrix3f(Rd, E1, E2).determinant();
		float gamma = Matrix3f(Rd, E1, S).determinant() / Matrix3f(Rd, E1, E2).determinant();
		if(t > 0){
			if(beta >= 0 && beta <= 1 && gamma <= 1 && gamma >= 0 && beta + gamma <= 1){
				//在三角形内部
				
				if(t >= tmin && t <= hit.getT()){
					hit.set(t, this->material, this->normal);
					return true;
				} else{
					return false;//应该返回false还是true？
				}
			} else {
				return false;
			}
		} else {
			return false;
		}
	}
	Vector3f normal;
	Vector3f vertices[3];

    void drawGL() override {
        Object3D::drawGL();
        glBegin(GL_TRIANGLES);
        glNormal3fv(normal);
        glVertex3fv(vertices[0]); glVertex3fv(vertices[1]); glVertex3fv(vertices[2]);
        glEnd();
    }

protected:
};

#endif //TRIANGLE_H
