#ifndef MBA_DECOMP_H_
#define MBA_DECOMP_H_

#include <iostream>
#include <list>
#include <unordered_map>
#include<fstream>
#include <map>
#include <cstdint>
#include <string>
#include <vector>
#include <algorithm>
#include "tgraph.h"

using namespace std;

// struct IteratorHash {
// 	std::size_t operator()(const std::list<uint32_t>::iterator& it) const {
// 		return std::hash<const void*>()(static_cast<const void*>(&(*it)));
// 	}
// };

// // 定义相等函数
// struct IteratorEqual {
// 	bool operator()(const std::list<uint32_t>::iterator& a, const std::list<uint32_t>::iterator& b) const {
// 		return &(*a) == &(*b);
// 	}
// };

typedef struct final {
    uint32_t vid;
    uint32_t eid;
  } ArrayEntry;
// edge type
typedef std::pair<uint32_t, uint32_t> EdgT;
// triangle type
typedef std::tuple<uint32_t, uint32_t, uint32_t> Triangle;
typedef std::vector<uint32_t> VI;

class Mba final {
public:
    
    Mba(const uint32_t n, const uint32_t l, const std::string& fn);
    Mba(const Mba&) = delete;
    Mba& operator=(const Mba&) = delete;
    ~Mba() {}
    void loadGraph(const std::string& dataset_path);
    void trussDecomp();
    void KdeltaTrussDecomp();
    vector<vector<int>> findkdCommunityForQuery(int query, int k, uint32_t delta);
    void writeToFile(const std::string& file_name);
    void get_KSpan_dict(int k, map<uint32_t, vector<int>>& kspan_dict){
        for (int i = 0; i < tg_.m_; i++) {
            uint32_t sp = kspan_[k][i];
            if (sp != UINT32_MAX) {
                kspan_dict[sp].push_back(i);
            }
            
        }
    }
    inline void getminmax(uint32_t a, uint32_t b, uint32_t c){
        if (a > b) a ^= b ^= a ^= b;
        if (a > c) a ^= c ^= a ^= c;
        if (b > c){maxe = b; mine = a;}
        else {maxe = c; mine = a;}
    }

    inline Triangle sortTri(uint32_t a, uint32_t b, uint32_t c){
        if (a > b) a ^= b ^= a ^= b;
        if (a > c) a ^= c ^= a ^= c;
        if (b > c) b ^= c ^= b ^= c;
        return {a, b, c};
    }

    void print_info(uint32_t t_, uint32_t n_, uint32_t m_) {
        cout << "-----------info--------------" << endl;
        cout << "Number of timestamps: " << t_ << endl;
        cout << "Number of vertices: " << n_ << endl;
        cout << "Number of temporal edges: " << m_ << endl;
        cout << "-----------------------------" << endl;
    }
    TGraph tg_;
    vector<vector<uint32_t>> kspan_;
    int kmax() { return k_max; }

private:
    const uint32_t l_;
    const uint32_t n_;
    uint32_t ttshreld_;
    uint32_t m_;
    vector<int> trn_;
    // data members
    int k_max;	//the # maximum k value
    uint32_t t_max; //the #maximum timestamps
    // the k-support
    VI ks_;
    //the support
    VI s_;
    // the edge peeling order
    VI ord_;
    //the position of edges in order 
    VI pos;
    uint32_t mine, maxe;
};
#endif
