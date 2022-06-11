#ifndef KDTREE_H
#define KDTREE_H

#define MIN_BOX_SIZE 10
#define MAX_OBJ_BOX 100

#include <vector>
#include <list>
#include <algorithm>
#include <queue>
#include <numeric>

#include "scene/Scene.hpp"
#include "Ray.hpp"
#include "PathPhotonMapping.hpp"

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

    bool IsHit(const PathPhotonMapping::Ray& r, KdNode* n) {
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

    std::list<KdNode*> find_Node(const PathPhotonMapping::Ray& r) {
        std::list<KdNode*> result;
        find_Node(r, root, result);
        return result;
    }

    void find_Node(const PathPhotonMapping::Ray& r, KdNode* n, std::list<KdNode*>& result) {
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
    PathPhotonMapping::Vec3 position;
    PathPhotonMapping::Vec3 r;
    PathPhotonMapping::Ray in, out;
    PathPhotonMapping::Vec3 norm;
    float operator[](int index) const {
        return position[index];
    }
    float& operator[](int index) {
        return position[index];
    }
};

struct KDNode {
    int axis = -1;
    int idx = -1;
    int leftChildIdx = -1;
    int rightChildIdx = -1;

    KDNode() {}
};

class KDTree {
private:
    std::vector<KDNode> nodes;
    std::vector<photon> photons;
    int nPoints;

    void buildNode(int* indices, int n_points, int depth) {
        if (n_points <= 0) return;

        int axis = depth % 3;

        std::sort(indices, indices + n_points, [&](const int idx1, const int idx2) {
            return photons[idx1][axis] < photons[idx2][axis];
            });

        int mid = (n_points - 1) / 2;

        int parentIdx = nodes.size();
        KDNode node;
        node.axis = axis;
        node.idx = indices[mid];
        nodes.push_back(node);

        int leftChildIdx = nodes.size();
        buildNode(indices, mid, depth + 1);

        if (leftChildIdx == nodes.size()) {
            nodes[parentIdx].leftChildIdx = -1;
        }
        else {
            nodes[parentIdx].leftChildIdx = leftChildIdx;
        }

        const int rightChildIdx = nodes.size();
        buildNode(indices + mid + 1, n_points - mid - 1, depth + 1);

        if (rightChildIdx == nodes.size()) {
            nodes[parentIdx].rightChildIdx = -1;
        }
        else {
            nodes[parentIdx].rightChildIdx = rightChildIdx;
        }
    }

    static inline float distance(const photon& p1, const photon& p2) {
        float x = p1[0] - p2[0];
        float y = p1[1] - p2[1];
        float z = p1[2] - p2[2];
        return std::pow(x, 2) + std::pow(y, 2) + std::pow(z, 2);
    }

    void searchKNearestNode(int nodeIdx, const photon& p, int k, std::priority_queue<std::pair<float, int>>& queue) {
        if (nodeIdx == -1 || nodeIdx >= nodes.size()) return;

        KDNode node = nodes[nodeIdx];

        photon median = photons[node.idx];

        float dist2 = distance(p, median);
        queue.push(std::make_pair(dist2, node.idx));

        if (queue.size() > k) {
            queue.pop();
        }

        bool isLower = p[node.axis] < median[node.axis];
        if (isLower) {
            searchKNearestNode(node.leftChildIdx, p, k, queue);
        }
        else {
            searchKNearestNode(node.rightChildIdx, p, k, queue);
        }

        float dist_to_siblings = median[node.axis] - p[node.axis];
        if (queue.top().first > dist_to_siblings * dist_to_siblings) {
            if (isLower) {
                searchKNearestNode(node.rightChildIdx, p, k, queue);
            }
            else {
                searchKNearestNode(node.leftChildIdx, p, k, queue);
            }
        }
    }

    void buildTree() {
        std::vector<int> indices(nPoints);
        std::iota(indices.begin(), indices.end(), 0);

        buildNode(indices.data(), nPoints, 0);
    }
public:
    KDTree() {}

    void setPoints(const std::vector<photon> points) {
        this->photons = points;
        this->nPoints = points.size();
        buildTree();
    }

    std::vector<int> searchKNearest(const photon& p, int k, float& max) {
        std::priority_queue<std::pair<float, int>> queue;
        searchKNearestNode(0, p, k, queue);

        std::vector<int> ret(queue.size());
        max = 0;
        for (int i = 0; i < ret.size(); ++i) {
            auto p = queue.top().second;
            auto dis = queue.top().first;
            ret[i] = p;
            max = std::max(max, std::sqrt(dis));
            queue.pop();
        }
        return ret;
    }
};

#endif