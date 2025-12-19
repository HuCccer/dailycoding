#pragma once
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <list>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <cstdint>
#include "Tgraph.h"
#include "anchor_union_find.h"

using namespace std;

typedef vector<int> SGN;

struct pair_hash {
    size_t operator()(const pair<int, int>& p) const {
        return std::hash<int>()(p.first) ^ (std::hash<int>()(p.second) << 1);
    }
};

struct TN {
    int id_; uint32_t delta_;
    int father_;              // father TN id, -1 if none
    vector<int> children_;   // children TN ids
    vector<int> SGNids_;
    // TN() {};
    TN(uint32_t delta, int id) : delta_(delta), id_(id), father_(-1) {};
};


struct Forest {
    vector<int> etoTreeNID;  // edge ID → Tree node ID
    vector<SGN> idSGN; // Super node ID → SGN object
    unordered_map<int, TN> idTN; // Super node ID → Tree node object
};


/***************************************************
 * Abstract Base Class: KDEquiForests
 ***************************************************/
class KDEquiForests{
public:
    int num_; // number of forest 
    int m_; // number of edges
    TGraph &tg_;
    // Original vertex → set of super nodes
    vector<Forest> KDEForests;
    AnchorUnionFind auf;
    KDEquiForests(int num, int m, TGraph& tg) : num_(num), m_(m), auf(m), tg_(tg), KDEForests(num_){};
    KDEquiForests(const KDEquiForests&) = delete;
    KDEquiForests& operator=(const KDEquiForests&) = delete;
    ~KDEquiForests() {}


    void constructIndexForK(map<uint32_t, vector<int>>& deltaEdgeLists, vector<uint32_t>& kspan, TGraph& graph, Forest& kdforest, int trussness);
    void mergeNode(int nx, int ny, Forest& forest, unordered_map<int,int>& SIDtoTID);
    vector<vector<int>> findkdCommunityForQuery(int query, int k, uint32_t delta);

};
