#include <vector>
#include <algorithm>
using namespace std;

class AnchorUnionFind {
public:
    vector<int> parent_;  // 父节点
    vector<int> rank_;    // 树的秩（近似高度）
    vector<int> anchor_; // anchor_ 节点

public:
    // 构造函数：初始化 n 个元素，每个元素的父节点是自己
    AnchorUnionFind(int m) {
        parent_.resize(m);
        rank_.resize(m, 0);
        anchor_.resize(m);
        for (int i = 0; i < m; i++) {
            parent_[i] = i; anchor_[i] = i;
        }
    }

    // 查找操作（Find with Path Compression）
    int find(int x) {
        if (parent_[x] != x) {
            parent_[x] = find(parent_[x]); // 路径压缩
        }
        return parent_[x];
    }

    // 合并操作（Union by rank_）
    void update(int x, int y, vector<int>& sp) {
        int px = find(x);
        int py = find(y);

        if (px == py) return; // 本来就同一个集合
        int new_root = -1;
        // 按秩合并，让低秩树挂到高秩树下
        if (rank_[px] < rank_[py]) {
            parent_[px] = py;
            new_root = py;
        } else if (rank_[py] < rank_[px]) {
            parent_[py] = px;
            new_root = px;
        } else {
            parent_[py] = px;
            new_root = px;
            rank_[px]++;  // 合并后树高 +1
        }
        int ax = anchor_[px], ay = anchor_[py];
        if (sp[ax] <= sp[ay]) {
            anchor_[new_root] = ay;
        } else {
            anchor_[new_root] = ax;
        }
    }

    void init() {
        int m = parent_.size();
        for (int i = 0; i < m; i++) {
            parent_[i] = i; anchor_[i] = i; 
        }
        fill(rank_.begin(), rank_.end(), 0);
    }

};
