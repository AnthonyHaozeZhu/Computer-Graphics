#ifndef KDTREE_H
#define KDTREE_H

#define MIN_BOX_SIZE 10
#define MAX_OBJ_BOX 100

#include <vector>
#include <list>
#include <algorithm>
#include <queue>
#include <numeric>
#include <eigen3/Eigen/Eigen>

#include "scene/Scene.hpp"
#include "Ray.hpp"
#include "PhotonMapping.hpp"

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

struct KdNode {
    KdNode() {
        min = std::vector<float>(3);
        max = std::vector<float>(3);
        left = right = nullptr;
    }
    KdNode(Box* box) {
        min = std::vector<float>(3, INFINITY);
        max = std::vector<float>(3, -INFINITY);
        Update(box);
        left = right = nullptr;
        Box_list.push_back(box);
    }

    KdNode(KdNode* n, int index, bool flag, int v) {
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
    KdNode* left, * right;
    bool is_leaf = true;
    float r;
};

struct KdTree {
    KdTree() { root = new KdNode; }
    void Insert(std::vector<Box*>& boxs, int s, int e, KdNode* n) {
        for (int i = s; i < e; i++) {
            n->Update(boxs[i]);
        }
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
            n->left = new KdNode;
            n->right = new KdNode;
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

    bool IsHit(const PhotoMapping::Ray& r, KdNode* n) {
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

    std::list<KdNode*> find_Node(const PhotoMapping::Ray& r) {
        std::list<KdNode*> result;
        find_Node(r, root, result);
        return result;
    }

    void find_Node(const PhotoMapping::Ray& r, KdNode* n, std::list<KdNode*>& result) {
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

    KdNode* root = nullptr;
    int MAX_SIZE = 3;
    int index = 0;
};

struct photon {
    PhotoMapping::Vec3 position;
    PhotoMapping::Vec3 r;
    PhotoMapping::Ray in, out;
    PhotoMapping::Vec3 norm;
    float operator[](int index) const {
        return position[index];
    }
    float& operator[](int index) {
        return position[index];
    }
};

template<const size_t DIM>
using VecXf = Eigen::Vector<float, DIM>;

template<const size_t DIM>
struct KDNode {
    int axis = -1;
    int index = -1;
    VecXf<DIM> pos = {};
    KDNode* left = nullptr, * right = nullptr;
};



template<const size_t DIM>
struct KDTree {
    void buildTree(const std::vector<VecXf<DIM>>& points) {
        std::vector<int> indexs(points.size());
        for (int i = 0; i < indexs.size(); i++) indexs[i] = i;
        this->points = points;
        root = Insert(0, indexs.size(), indexs, 0);
    }

    void buildTree(std::vector<VecXf<DIM>>&& points) {
        std::vector<int> indexs(points.size());
        for (int i = 0; i < indexs.size(); i++) indexs[i] = i;
        this->points = points;
        root = Insert(0, indexs.size(), indexs, 0);
    }

    KDNode<DIM>* Insert(int s, int e, std::vector<int>& indexs, int depth) {
        if (s >= e) return nullptr;
        if (e - s == 1) {
            return new KDNode<DIM>{ depth % DIM , indexs[s], points[indexs[s]] };
        }
        std::sort(indexs.begin() + s, indexs.begin() + e, [&](const int& i, const int& j) {
            return points[i][depth % DIM] < points[j][depth % DIM];
            });
        int mid = (e + s) / 2;
        KDNode<DIM>* n = new KDNode<DIM>{ depth % DIM , indexs[mid], points[indexs[mid]] };
        n->left = Insert(s, mid, indexs, depth + 1);
        n->right = Insert(mid + 1, e, indexs, depth + 1);
        return n;
    }

    using Queue = std::priority_queue<std::pair<float, int>>;

    template<typename Point>
    std::vector<int> nearstSearch(const Point& points, float& max_r, int k) {
        VecXf<DIM> p;
        for (int i = 0; i < DIM; i++) {
            p[i] = points[i];
        }
        Queue queue;
        nearstSearch(root, p, queue, k);
        std::vector<int> result_indexs(queue.size());
        max_r = 0.f;
        for (int i = 0; i < result_indexs.size(); ++i) {
            int index = queue.top().second;
            float dis = queue.top().first;
            result_indexs[i] = index;
            max_r = std::max(max_r, std::sqrt(dis));
            queue.pop();
        }
        return result_indexs;
    }

    template<typename Point>
    int nearstSearch(const Point& p, float& max_r) {
        return nearstSearch(p, max_r, 1)[0];
    }

    static inline float distance(const VecXf<DIM>& p1, const VecXf<DIM>& p2) {
        float dis2 = 0.f;
        for (int i = 0; i < DIM; i++) {
            dis2 += std::pow(p1[i] - p2[i], 2);
        }
        return dis2;
    }

    void nearstSearch(KDNode<DIM>* n, const VecXf<DIM>& p, Queue& queue, int k) {
        if (n == nullptr) return;
        float dis = distance(n->pos, p);
        queue.push(std::make_pair(dis, n->index));
        if (queue.size() > k) {
            queue.pop();
        }
        dis = std::pow(p[n->axis] - n->pos[n->axis], 2);
        if (p[n->axis] < n->pos[n->axis]) {
            nearstSearch(n->left, p, queue, k);
            if (n->right && (dis < queue.top().first)) {
                nearstSearch(n->right, p, queue, k);
            }
        }
        else {
            nearstSearch(n->right, p, queue, k);
            if (n->left && (dis < queue.top().first)) {
                nearstSearch(n->left, p, queue, k);
            }
        }

    }
    KDNode<DIM>* root = nullptr;
    std::vector<VecXf<DIM>> points;
};
#endif