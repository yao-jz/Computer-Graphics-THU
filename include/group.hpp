
#ifndef GROUP_H
#define GROUP_H

#include "object3d.hpp"
#include "ray.hpp"
#include "hit.hpp"
#include <iostream>
#include "bvhnode.hpp"
#include <vector>

using namespace std;

// TODO (PA2): Implement Group - copy from PA1 finish
class Group : public Object3D
{

public:
    Group()
    {
        group.resize(10);
    }

    explicit Group(int num_objects)
    {
        group.resize(num_objects);
        // 构建的时候用这个
    }

    ~Group() override
    {
    }

    bool intersect(const Ray &r, Hit &h, float tmin) override
    {
        if (bvh_root)
            return bvh_root->intersect(r, h, tmin);

        bool result = false;
        for (int i = 0; i < group.size(); i++)
        {
            //如果result是true的话，反过来就不会调用intersect这个函数了
            result = group[i]->intersect(r, h, tmin) || result;
        }
        return result;
    }

    virtual bool bounding_box(double time0, double time1, aabb &output_box) const override
    {
        // 处理叶结点也是group的情况
        if (group.empty())
            return false;
        aabb temp_box;
        bool first_box = true;
        for (const auto &object : group)
        {
            if (!object->bounding_box(time0, time1, temp_box))
                return false;
            output_box = first_box ? temp_box : surrounding_box(output_box, temp_box);
            first_box = false;
        }
        return true;
    }

    void init()
    {
        aabb temp_box;
        if (!bvh_root and this->bounding_box(0, 0, temp_box))
            bvh_root = new bvh_node(group, 0, group.size(), 0, 0);
    }

    void drawGL() override
    {
        for (int i = 0; i < group.size(); i++)
        {
            group[i]->drawGL();
        }
    }

    void addObject(int index, Object3D *obj)
    {
        group[index] = obj;
    }

    int getGroupSize()
    {
        return group.size();
    }

    bvh_node *bvh_root = nullptr;
    vector<Object3D *> group;
};

#endif
