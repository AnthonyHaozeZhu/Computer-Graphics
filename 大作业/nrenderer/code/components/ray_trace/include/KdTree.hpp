#ifndef __KDTREE_HPP__
#define __KDTREE_HPP__

#define MIN_BOX_SIZE 10
#define MAX_OBJ_BOX 100

#include <vector>
#include <list>

#include "scene/Scene.hpp"
#include "Ray.hpp"

#define SPHERE 100
#define TRIANGLE 101

struct Box {
    Box() {}
    Box(NRenderer::Sphere* sp) {
        type = SPHERE;
        this->sp = sp;
        float r = sp->radius;
        min = std::vector<float>{ sp->position.x - r, sp->position.y - r, sp->position.z - r };
        max = std::vector<float>{ sp->position.x + r, sp->position.y + r, sp->position.z + r };
    }
    Box(NRenderer::Triangle* tr) {
        type = TRIANGLE;
        this->tr = tr;
        min = std::vector<float>{ std::min(tr->v1.x, std::min(tr->v2.x, tr->v3.x)),
            std::min(tr->v1.y, std::min(tr->v2.y, tr->v3.y)),
            std::min(tr->v1.z, std::min(tr->v2.z, tr->v3.z)) };

        max = std::vector<float>{ std::max(tr->v1.x, std::max(tr->v2.x, tr->v3.x)),
            std::max(tr->v1.y, std::max(tr->v2.y, tr->v3.y)),
            std::max(tr->v1.z, std::max(tr->v2.z, tr->v3.z)) };
    }
    union
    {
        NRenderer::Sphere* sp;
        NRenderer::Triangle* tr;
    };
    int type = -1;
    std::vector<float> min;
    std::vector<float> max;
};

struct KDNode {
    KDNode() {
        min = std::vector<float>(3);
        max = std::vector<float>(3);
        left = right = nullptr;
    }
    KDNode(Box* box) {
        min = std::vector<float>(3, INFINITY);
        max = std::vector<float>(3, -INFINITY);
        Update(box);
        left = right = nullptr;
        Box_list.push_back(box);
    }

    KDNode(KDNode* n, int index, bool flag, int v) {
        min = std::vector<float>(3);
        max = std::vector<float>(3);
        for (int i = 0; i < min.size(); i++) {
            min[i] = n->min[i];
        }
        for (int i = 0; i < max.size(); i++) {
            max[i] = n->max[i];
        }
        if (flag) {
            max[index] = v;
        }
        else {
            min[index] = v;
        }
        left = right = nullptr;
    }

    void Update(Box* box) {
        for (int i = 0; i < box->min.size(); i++) {
            min[i] = std::min(box->min[i], min[i]);
        }
        for (int i = 0; i < box->max.size(); i++) {
            max[i] = std::max(box->max[i], max[i]);
        }
    }

    bool IsInNode(const Box* box) {
        for (int i = 0; i < box->min.size(); i++) {
            if (box->min[i] < min[i]) return false;
        }
        for (int i = 0; i < box->max.size(); i++) {
            if (box->max[i] > max[i]) return false;
        }
        return true;
    }

    void Insert(Box* box) {
        Update(box);
        Box_list.push_back(box);
    }

    std::vector<float> min;
    std::vector<float> max;
    std::list<Box*> Box_list;
    KDNode* left, * right;
    bool is_leaf = true;
    float r;
};

struct KDTree {
    KDTree() { root = new KDNode; }
    void Insert(std::vector<Box*>& boxs, int s, int e, KDNode* n) {
        if (s == e) return;
        else if (e - s <= MAX_SIZE) {
            for (int i = s; i < e; i++) {
                n->Insert(boxs[i]);
            }
        }
        else {
            for (int i = s; i < e; i++) {
                n->Update(boxs[i]);
            }
            n->is_leaf = false;
            n->left = new KDNode;
            n->right = new KDNode;
            std::sort(boxs.begin() + s, boxs.begin() + e, [&](const Box* A, const Box* B) {
                if (index == 1) {
                    return A->min[0] < B->min[0];
                }
                else if (index == 2) {
                    return A->min[1] < B->min[1];
                }
                else {
                    return A->min[2] < B->min[2];
                }
                });
            index = (index + 1) % 3;
            int mid = (s + e) / 2;
            Insert(boxs, s, mid, n->left);
            Insert(boxs, mid, e, n->right);
        }
    }
    void Insert(std::vector<NRenderer::Sphere>& sps, std::vector<NRenderer::Triangle>& trs) {
        std::vector<Box*> boxs(sps.size() + trs.size());
        int index = 0;
        for (auto& sp : sps) {
            boxs[index++] = new Box(&sp);
        }
        for (auto& tr : trs) {
            boxs[index++] = new Box(&tr);
        }
        Insert(boxs, 0, boxs.size(), root);
    }

    bool IsHit(const RayTrace::Ray& r, KDNode* n) {
        auto inv_dir = 1.0f / r.direction;
        NRenderer::Vec3 cube_min(n->min[0], n->min[1], n->min[2]), cube_max(n->max[0], n->max[1], n->max[2]);

        auto tMin = (cube_min - r.origin) * inv_dir;
        auto tMax = (cube_max - r.origin) * inv_dir;
        auto t1 = min(tMin, tMax);
        auto t2 = max(tMin, tMax);
        float tNear = std::max(std::max(t1.x, t1.y), t1.z);
        float tFar = std::min(std::min(t2.x, t2.y), t2.z);

        return  tNear < tFar;
    }

    std::list<KDNode*> find_Node(const RayTrace::Ray& r){
        std::list<KDNode*> result;
        find_Node(r, root, result);
        return result;
    }

    void find_Node(const RayTrace::Ray& r, KDNode* n, std::list<KDNode*>& result) {
        if (IsHit(r, n)) {
            if (n->is_leaf) {
                result.push_back(n);
            }
            else {
                find_Node(r, n->left, result);
                find_Node(r, n->right, result);
            }
        }
    }

    KDNode* root = nullptr;
    int MAX_SIZE = 3;
    int index = 0;
};

#endif