#include <queue>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <list>
#include "KDEquiForests.h"


using namespace std;


/***************************************************
 * constructIndex
 ***************************************************/

void KDEquiForests::constructIndexForK(map<uint32_t, vector<int>>& deltaEdgeList, vector<uint32_t>& ksp, TGraph& tg, Forest& forest, int k) {
    unordered_map<int, vector<pair<int, uint32_t>>> edgeigd; // edge id to edge object
    forest.etoTreeNID.resize(tg.m(), -1);
    list<tuple<int, int, uint32_t>> W;
    vector<bool> visited(tg.m());
    vector<int> SGNtodelta(tg.m());
    unordered_map<pair<int, int>, uint32_t, pair_hash> Buffer;
    unordered_map<int,int> SIDtoTID; 
    queue<int> Qk;
    int tnid = 0; // kdequitruss node ID
    int last_tnid = 0;
    for (auto& tpair : deltaEdgeList) {
        uint32_t delta= tpair.first;
        vector<int>& edgelist = tpair.second;
        while (!edgelist.empty()) {
            int e0 = edgelist.back();
            edgelist.pop_back();  
            if (visited[e0]) continue;
            Qk.push(e0);
            forest.idSGN.emplace_back(SGN());
            SGN& Ed = forest.idSGN.back();
            SGNtodelta[tnid] = delta;
            forest.idTN.emplace(tnid, TN(delta, tnid));
            forest.idTN.at(tnid).SGNids_.push_back(tnid);
            SIDtoTID[tnid] = tnid;
            visited[e0] = true;
            SIDtoTID[tnid] = tnid;
            while (!Qk.empty()) {
                int e = Qk.front(); Qk.pop();
                Ed.emplace_back(e);
                auto etolist = edgeigd.find(e);
                if (etolist != edgeigd.end()) {
                    for (const auto& p : etolist->second) {
                        // a kdequitruss set id
                        int sid = p.first;
                        if (sid == tnid) continue;
                        // get delta of the common triangle
                        uint32_t tdelta = p.second;
                        // postpone parent-child relationship setting
                        auto [it, inserted] = Buffer.try_emplace({sid, tnid}, tdelta);
                        if (!inserted) {
                            if (it->second <= delta) continue;
                            it->second = min(it->second, tdelta);
                        } 
                        if (it->second > delta) continue;
                        int rnid = auf.find(sid);
                        int anchorid = auf.anchor_[rnid];
                        int preid = SIDtoTID[anchorid];
                        int curid = SIDtoTID[tnid];
                        TN& preTN = forest.idTN.at(preid);
                        TN& currTN = forest.idTN.at(curid);
                        if (preTN.father_ == curid) continue;
                        // set parent-child relationship    
                        if (preTN.father_ == -1) {
                            preTN.father_ = currTN.id_;
                            currTN.children_.push_back(preTN.id_);
                        } else {
                            mergeNode(preTN.father_, curid, forest, SIDtoTID);
                        }  
                    }
                }
                std::vector<std::pair<uint32_t, uint32_t>> tris = tg.GetTriangles(e, k);
                // enumerate triangle neighbors
                for (const auto& tri : tris) {
                    const uint32_t e1 = tri.first; const uint32_t e2 = tri.second;
                    uint32_t sp1 = ksp[e1], sp2 = ksp[e2], deltatri = tg.GetMst(tg.tau_[e], tg.tau_[e1], tg.tau_[e2]);
                    uint32_t maxsp = max(sp1, sp2);
                    if (maxsp < delta) continue;
                    // function to process_edge
                    // process edge e1
                    if (sp1 == maxsp && !visited[e1]) {
                        if (sp1 == delta && deltatri <= delta) {
                            Qk.push(e1);
                            visited[e1] = true;
                        } else {
                            edgeigd[e1].push_back({tnid, deltatri});
                        }
                    }
                    // process edge e2
                    if (sp2 == maxsp && !visited[e2]) {
                        if (sp2 == delta && deltatri <= delta) {
                            Qk.push(e2);
                            visited[e2] = true;
                        } else {
                            edgeigd[e2].push_back({tnid, deltatri});
                        }
                    }
                }
                edgeigd.erase(e);
            }
            tnid++;
        }
        list<tuple<int, int, uint32_t>> delay_re; {
            for (auto& [re, t] : Buffer) {
                auto& [x, y]  = re;
                if (t > delta) { 
                    delay_re.push_back({x, y, t});
                } else {
                    auf.update(x, y, SGNtodelta);
                }
            }
        }
        Buffer.clear();
        for (auto it = W.begin(); it != W.end();) {
            auto [x, y, t] = *it;
            int rx = auf.find(x), ry = auf.find(y);
            if (rx == ry) {
                it = W.erase(it);
                continue;
            }
            if (t > delta) {
                ++it;
                continue;
            }
            int ax = auf.anchor_[rx], ay = auf.anchor_[ry];
            int tidx = SIDtoTID.at(ax), tidy = SIDtoTID.at(ay);
            uint32_t deltax = forest.idTN.at(tidx).delta_, deltay = forest.idTN.at(tidy).delta_;
            uint32_t maxdelta = max(deltax, deltay);
            if (maxdelta >= t) {
                if (maxdelta == deltax && maxdelta == deltay) {
                    mergeNode(tidx, tidy, forest, SIDtoTID);
                } else if (maxdelta == deltax) {
                    forest.idTN.at(tidy).father_ = tidx;
                    forest.idTN.at(tidx).children_.push_back(tidy);
                } else if (maxdelta == deltay) {
                    forest.idTN.at(tidx).father_ = tidy;
                    forest.idTN.at(tidy).children_.push_back(tidx);
                }
                auf.update(x, y, SGNtodelta);
                it = W.erase(it);
            }  
        }
        W.splice(W.end(), delay_re);
        for (int i = last_tnid; i < tnid; i++) {
            for (const int& e : forest.idSGN[i]) {
                forest.etoTreeNID[e] = SIDtoTID[i];
            }
        }
        last_tnid = tnid;
    }
    unordered_map<uint64_t, uint32_t> reduced_W; {
        for (auto& [x, y, t] : W) {
            int rx = auf.find(x), ry = auf.find(y);
            if (rx > ry) swap(rx, ry);
            uint64_t key = (uint64_t(rx) << 32) | uint32_t(ry);
            auto [it, inserted] = reduced_W.try_emplace(key, t);
            if (!inserted) {
                it->second = min(it->second, t);
            }    
        }
    }
    vector<tuple<int, int, uint32_t>> sorted_reduced_W; {
        for (auto& [key, t] : reduced_W) {
            int rx = key >> 32; int ry = key & 0xffffffff;
            sorted_reduced_W.push_back({t, rx, ry});
        }
    }
    sort(sorted_reduced_W.begin(), sorted_reduced_W.end(),
         [](const std::tuple<int, int, uint32_t>& a, const std::tuple<int, int, uint32_t>& b) {
            return std::get<0>(a) > std::get<0>(b);  // z descending
        });

    for (auto& [t, rx, ry] : sorted_reduced_W) {
        int nrx = auf.find(rx), nry = auf.find(ry);
        int ax = auf.anchor_[nrx], ay = auf.anchor_[nry];
        int tidx = SIDtoTID.at(ax), tidy = SIDtoTID.at(ay);
        uint32_t deltax = forest.idTN.at(tidx).delta_, deltay = forest.idTN.at(tidy).delta_;
        uint32_t maxdelta = max(deltax, deltay);
        if (maxdelta < t) {
            forest.idSGN.emplace_back(SGN());
            SGNtodelta[tnid] = t;
            forest.idTN.emplace(tnid, TN(t, tnid));
            forest.idTN.at(tnid).SGNids_.push_back(tnid);
            SIDtoTID[tnid] = tnid;
            forest.idTN.at(tnid).children_.push_back(tidx);
            forest.idTN.at(tidx).father_ = tnid;
            forest.idTN.at(tnid).children_.push_back(tidy);
            forest.idTN.at(tidy).father_ = tnid;
            auf.update(nrx, tnid, SGNtodelta);
            auf.update(nry, tnid, SGNtodelta);
        } else {
            if (maxdelta == deltax && maxdelta == deltay) {
                mergeNode(tidx, tidy, forest, SIDtoTID);
            } else if (maxdelta == deltax) {
                forest.idTN.at(tidy).father_ = tidx;
                forest.idTN.at(tidx).children_.push_back(tidy);
            } else if (maxdelta == deltay) {
                forest.idTN.at(tidx).father_ = tidy;
                forest.idTN.at(tidy).children_.push_back(tidx);
            }
            auf.update(rx, ry, SGNtodelta);
        }
    }
    
    deltaEdgeList.clear();
    auf.init();
}

void KDEquiForests::mergeNode(int xid, int yid, Forest& forest, unordered_map<int,int>& SIDtoTID) {
    TN& nodex = forest.idTN.at(xid), nodey = forest.idTN.at(yid);
    for (int sid : nodey.SGNids_) {
        nodex.SGNids_.push_back(sid); SIDtoTID[sid] = xid;    
    }
    for (int cid : nodey.children_) {
        forest.idTN.at(cid).father_ = xid;
        nodex.children_.push_back(cid);
    }
    forest.idTN.erase(yid); 
}


vector<vector<int>> KDEquiForests::findkdCommunityForQuery(int query, int k, uint32_t delta) {
    unordered_set<int> visited;
    Forest& Forest = KDEForests[k];
    queue<int> q;
    vector<vector<int>> result;
    for (const auto& ae : tg_.adj_[query]) {
        int e = ae.eid;
        int treeid = Forest.etoTreeNID[e];
        if (visited.count(treeid)) continue;
        const TN* node = &Forest.idTN.at(treeid);
        if (node->delta_ > delta) continue;
        while (node->father_ != -1 && Forest.idTN.at(node->father_).delta_ <= delta) {
            node = &Forest.idTN.at(node->father_);
        }
        vector<int> Ai;
        q.push(node->id_);
        while (!q.empty()) {
            int treeid = q.front(); q.pop();
            visited.insert(treeid);
            const TN* node = &Forest.idTN.at(treeid);
            for (int sgnid : node->SGNids_) {
                Ai.insert(Ai.end(), Forest.idSGN[sgnid].begin(), Forest.idSGN[sgnid].end());
            }
            for (int cid : node->children_) {
                q.push(cid);
            }
        }
        result.emplace_back(std::move(Ai));
    }
    for (auto& community : result) {
		for (int& e : community) {
			cout<<e<<" ";
		}
		cout<<endl;
	}
    return result;
}

