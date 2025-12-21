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
    vector<unordered_map<int, uint32_t>> edgeigd(m_); // edge id to edge object
    forest.etoTreeNID.resize(m_, -1);
    map<int, unordered_set<pair<int, int>, pair_hash>> waitlist;
    vector<bool> visited(m_);
    vector<int> SGNtodelta(m_);
    unordered_set<pair<int, int>, pair_hash>Buffer;
    vector<int> SIDtoTID(m_); 
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
            forest.idTN[tnid].SGNids_.push_back(tnid);
            SIDtoTID[tnid] = tnid;
            visited[e0] = true;
            while (!Qk.empty()) {
                int e = Qk.front(); Qk.pop();
                Ed.emplace_back(e);
                if (!edgeigd[e].empty()) {
                    for (const auto& [sid, t] : edgeigd[e]) {
                        if (sid == tnid) continue; 
                        // postpone parent-child relationship setting
                        if (t > delta) { waitlist[t].insert({sid, tnid}); continue;} 
                        // set parent-child relationship    
                        int rnid = auf.find(sid); int anchorid = auf.anchor_[rnid];
                        int lowLevelId = SIDtoTID[anchorid]; int highLevelId = SIDtoTID[tnid];
                        TN& lowLevelTN = forest.idTN[lowLevelId], &highLevelTN = forest.idTN[highLevelId];
                        if (lowLevelTN.father_ == highLevelId) continue;
                        if (lowLevelTN.father_ == -1) {
                            lowLevelTN.father_ = highLevelId;
                            highLevelTN.children_.push_back(lowLevelId);
                        } else {
                            // cout<<"cur---"<<endl;
                            mergeNode(lowLevelTN.father_, highLevelId, forest, SIDtoTID);
                        }  
                        Buffer.insert({rnid, tnid});
                    }
                    edgeigd[e].clear();
                }  
                // find triangles containing the edge with ID eid
                std::vector<std::pair<uint32_t, uint32_t>> tris = tg.GetTriangles(e, k);
                // enumerate triangle neighbors
                for (const auto& [e1, e2] : tris) {
                    if (visited[e1] && visited[e2]) continue;
                    uint32_t sp1 = ksp[e1], sp2 = ksp[e2];
                    uint32_t maxsp = max(sp1, sp2), deltatri = tg.GetMst(tg.tau_[e], tg.tau_[e1], tg.tau_[e2]);
                    // function to process_edge
                    // process edge e1
                    if (sp1 == maxsp && !visited[e1]) {
                        if (sp1 == delta && deltatri <= delta) {
                            Qk.push(e1);
                            visited[e1] = true;
                        } else {
                            auto [it, inserted] = edgeigd[e1].try_emplace(tnid, deltatri);
                            if (!inserted) {
                                if (it->second <= delta) continue;
                                it->second = min(it->second, deltatri);
                            } 
                        }
                    }
                    // process edge e2
                    if (sp2 == maxsp && !visited[e2]) {
                        if (sp2 == delta && deltatri <= delta) {
                            Qk.push(e2);
                            visited[e2] = true;
                        } else {
                            auto [it, inserted] = edgeigd[e2].try_emplace(tnid, deltatri);
                            if (!inserted) {
                                if (it->second <= delta) continue;
                                it->second = min(it->second, deltatri);
                            } 
                        }
                    }
                }
            }
            tnid++;
        }
        for (auto& [x, y] : Buffer) {
            auf.update(x, y, SGNtodelta);
        }
        Buffer.clear();
        // cout<<"waitlist: "<<waitlist.size()<<endl;
        for (auto it = waitlist.begin(); it != waitlist.end();) {
            auto& [t, relationship] = *it;
            if (t > delta) break;
            for (auto re = relationship.begin(); re != relationship.end();) {
                int x = re->first, y = re->second;
                int rx = auf.find(x), ry = auf.find(y);
                if (rx == ry) {re = relationship.erase(re); continue; }
                int ax = auf.anchor_[rx], ay = auf.anchor_[ry];
                int tidx = SIDtoTID[ax], tidy = SIDtoTID[ay];
                uint32_t deltax = forest.idTN[tidx].delta_, deltay = forest.idTN[tidy].delta_;
                uint32_t maxdelta = max(deltax, deltay);
                if (maxdelta >= t) {
                    if (maxdelta == deltax && maxdelta == deltay) {
                        // cout<<"waitlist"<<endl;
                        mergeNode(tidx, tidy, forest, SIDtoTID);
                    } else if (maxdelta == deltax) {
                        forest.idTN[tidy].father_ = tidx;
                        forest.idTN[tidx].children_.push_back(tidy);
                    } else if (maxdelta == deltay) {
                        forest.idTN[tidx].father_ = tidy;
                        forest.idTN[tidy].children_.push_back(tidx);
                    }
                    auf.update(x, y, SGNtodelta);
                    re = relationship.erase(re); 
                    continue;
                }
                ++re;
            }
             if (relationship.size() == 0) it = waitlist.erase(it);
             else ++it;    
        }
        for (int i = last_tnid; i < tnid; i++) {
            for (const int& e : forest.idSGN[i]) {
                forest.etoTreeNID[e] = SIDtoTID[i];
            }
        }
        last_tnid = tnid;
    }
    // cout<<"waitlist: "<<waitlist.size()<<endl;
    for (auto& [t, relationships] : waitlist) {
        for (auto& [x, y] : relationships) {
            int rx = auf.find(x), ry = auf.find(y);
            if (rx == ry) continue;
            int ax = auf.anchor_[rx], ay = auf.anchor_[ry];
            // cout<<"t: "<<t<<", anchor:"<<SGNtodelta[ax]<<" "<<SGNtodelta[ay]<<endl;
            int tidx = SIDtoTID[ax], tidy = SIDtoTID[ay];
            uint32_t deltax = forest.idTN[tidx].delta_, deltay = forest.idTN[tidy].delta_;
            uint32_t maxdelta = max(deltax, deltay);
            if (maxdelta < t) {
                // cout<<"new_node"<<endl;
                forest.idSGN.emplace_back(SGN());
                SGNtodelta[tnid] = t;
                forest.idTN.emplace(tnid, TN(t, tnid));
                forest.idTN[tnid].SGNids_.push_back(tnid);
                SIDtoTID[tnid] = tnid;
                forest.idTN[tnid].children_.push_back(tidx);
                forest.idTN[tidx].father_ = tnid;
                forest.idTN[tnid].children_.push_back(tidy);
                forest.idTN[tidy].father_ = tnid;
                auf.update(rx, tnid, SGNtodelta);
                auf.update(ry, tnid, SGNtodelta);
                tnid++;
            } else {
                if (maxdelta == deltax && maxdelta == deltay) {
                    mergeNode(tidx, tidy, forest, SIDtoTID);
                } else if (maxdelta == deltax) {
                    forest.idTN[tidy].father_ = tidx;
                    forest.idTN[tidx].children_.push_back(tidy);
                } else if (maxdelta == deltay) {
                    forest.idTN[tidx].father_ = tidy;
                    forest.idTN[tidy].children_.push_back(tidx);
                }
                auf.update(rx, ry, SGNtodelta);
            }
        }
       
    }
   
    deltaEdgeList.clear();
    auf.init();
}

void KDEquiForests::mergeNode(int xid, int yid, Forest& forest, vector<int>& SIDtoTID) {
    TN& nodex = forest.idTN[xid], nodey = forest.idTN[yid];
    for (int sid : nodey.SGNids_) {
        nodex.SGNids_.push_back(sid); SIDtoTID[sid] = xid;    
    }
    for (int cid : nodey.children_) {
        forest.idTN[cid].father_ = xid;
        nodex.children_.push_back(cid);
    }

    // for (int sid : nodex.SGNids_) {
    //     cout<< sid<<" ";
    // }
    // cout<<endl;
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
        const TN* node = &Forest.idTN[treeid];
        if (node->delta_ > delta) continue;
        while (node->father_ != -1 && Forest.idTN[node->father_].delta_ <= delta) {
            node = &Forest.idTN[node->father_];
        }
        vector<int> Ai;
        q.push(node->id_);
        while (!q.empty()) {
            int treeid = q.front(); q.pop();
            visited.insert(treeid);
            const TN* node = &Forest.idTN[treeid];
            for (int sgnid : node->SGNids_) {
                Ai.insert(Ai.end(), Forest.idSGN[sgnid].begin(), Forest.idSGN[sgnid].end());
            }
            for (int cid : node->children_) {
                q.push(cid);
            }
        }
        result.emplace_back(std::move(Ai));
    }
    int size = 0;
    for (auto& community : result) {
		size += community.size();
		// for (int& e : community) {
		// 	cout<<e<<" ";
		// }
	}
	cout<<size<<endl;
    return result;
}

