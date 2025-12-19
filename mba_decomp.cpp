#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <numeric>
#include <utility>
#include <queue>
#include <unordered_set>
#include "mba_decomp.h"

#define ASSERT(truth) \
    if (!(truth)) { \
      std::cerr << "\x1b[1;31mASSERT\x1b[0m: " \
                << "LINE " << __LINE__ \
                << ", " << __FILE__ \
                << std::endl; \
      std::exit(EXIT_FAILURE); \
    } else

#define ASSERT_MSG(truth, msg) \
    if (!(truth)) { \
      std::cerr << "\x1b[1;31mASSERT\x1b[0m: " \
                << "LINE " << __LINE__ \
                << ", " << __FILE__ << '\n' \
                << "\x1b[1;32mINFO\x1b[0m: " << msg \
                << std::endl; \
      std::exit(EXIT_FAILURE); \
    } else

vector<unordered_map<uint32_t, bool>> Mc;
vector<vector<Triangle>> time2triangle;

Mba::Mba(const uint32_t n, const uint32_t l, const std::string& fn)
    : n_(n), l_(l), tg_(n, l, trn_) {
    ASSERT_MSG(1 <= l_ && l_ < (static_cast<uint32_t>(1) << 29),
                "it is required 64 <= l <= 2^29 for the ease of implementation");
    trn_   = std::vector<int>(l_, 0);
    loadGraph(fn);
}


// truss decomposition and the corresponding order
void Mba::loadGraph(const std::string& file_name) {
    std::ifstream infile(file_name, std::ios::in);
    ASSERT_MSG(infile.is_open(), "cannot open the file");
    // read the size of the graph
    uint32_t num_of_timestamp = 0, num_of_vertex = 0, temporal_edges = 0;
    infile >> num_of_timestamp >> num_of_vertex >> temporal_edges;
	ttshreld_ = num_of_timestamp; time2triangle.resize(ttshreld_);
    ASSERT(temporal_edges <= l_ && n_ == num_of_vertex);
    print_info(num_of_timestamp, num_of_vertex, temporal_edges);
    cout << endl << "Load Graph.." << endl;
    ASSERT_MSG(!infile.eof(), "invalid graph file");
    // read the edges
    map<EdgT,vector<uint32_t>> ET;
    for (int te = 0; te < temporal_edges; te++){
        uint32_t t, v1, v2;
        infile >> t >> v1 >> v2;
        if (v1 > v2) std::swap(v1, v2);
        ET[{v1,v2}].emplace_back(t);
    }
    infile.close();
    uint32_t eid = 0;
    for(auto& [edge, times] : ET){
        ASSERT(tg_.LazyInsert(edge.first, edge.second, times) == eid++);
    }
	m_ = tg_.m();
    Mc.resize(m_);
     // rectify the graph
    tg_.Rectify();
    // clear ET
    decltype(ET)().swap(ET);
    cout << "Load graph successfully!" << endl;
}

void Mba::trussDecomp() {
 // truss decomposition
    // 1. compute the support of each edge by triangle listing
    // 1.1. define a total order over the vertices
    const auto pred = [this](const uint32_t v1, const uint32_t v2) {
    const size_t deg1 = tg_.adj_[v1].size();
    const size_t deg2 = tg_.adj_[v2].size();
        if (deg1 != deg2) return deg1 > deg2;
        else return v1 > v2;
    };
    VI verts(n_);
  	iota(verts.begin(), verts.end(), 0);
  	sort(verts.begin(), verts.end(), pred);
  	// 1.2. call the "forward" algorithm to list triangles
  	trn_.resize(m_); s_.resize(m_);
  	vector<vector<ArrayEntry>> A(n_);
    //delta-triangle list
	for (const uint32_t v : verts) {
    	for (const auto& ae : tg_.adj_[v]) {
			uint32_t u = ae.vid;
			uint32_t e = ae.eid;
			if (!pred(v, u)) continue;
			size_t pv = 0, pu = 0;
			while (pv < A[v].size() && pu < A[u].size()) {
				if (A[v][pv].vid == A[u][pu].vid) {
					++trn_[A[v][pv].eid]; ++trn_[A[u][pu].eid];
					++trn_[e];
					uint32_t interval = tg_.GetMst(tg_.tau_[e],tg_.tau_[A[v][pv].eid],tg_.tau_[A[u][pu].eid]);
					if(interval > t_max) t_max = interval;
					time2triangle[interval].emplace_back(sortTri(e,A[v][pv].eid,A[u][pu].eid));
					++pv; ++pu;
				} else if (pred(A[v][pv].vid, A[u][pu].vid)) {
				++pv;
				} else {
				++pu;
				}
			}
      		A[u].emplace_back(ArrayEntry{v, e});
		}
	}
	// 2. decomposition
	// 2.1. sort the edges according to their supports
	memcpy(s_.data(),trn_.data(),s_.size() * sizeof(uint32_t));
	const uint32_t max_sup = *max_element(trn_.cbegin(), trn_.cend());
	VI bin(max_sup + 1, 0);
  	for (uint32_t eid = 0; eid < m_; ++eid) ++bin[trn_[eid]];
	for (uint32_t i = 0, start = 0; i <= max_sup; ++i) {
		start += bin[i];
		bin[i] = start - bin[i];
	}

  	ord_.resize(m_);
	pos.resize(m_);
  	for (uint32_t eid = 0; eid < m_; ++eid) {
    	pos[eid] = bin[trn_[eid]];
    	ord_[pos[eid]] = eid;
    	++bin[trn_[eid]];
  	}
	rotate(bin.rbegin(), bin.rbegin() + 1, bin.rend());
  	bin[0] = 0;
  	// 2.2. peeling
  	ks_.resize(m_, 0);
	vector<bool> removed(m_, false);
  	int k = 0;
  	for (uint32_t i = 0; i < m_; ++i) {
    	k = max(k, trn_[ord_[i]]);
    	// ASSERT(bin[k] == i);
    	const uint32_t eid = ord_[i];
    	++bin[trn_[eid]];
    	removed[eid] = true;
    	// find triangles containing the edge with ID eid
    	vector<std::pair<uint32_t, uint32_t>> tris = tg_.GetTriangles(eid);
		// update ks_[eid]
		for (const auto& tri : tris) {
			const uint32_t e1 = tri.first;
			const uint32_t e2 = tri.second;
			if (trn_[e1] >= k && trn_[e2] >= k) ++ks_[eid];
			if (removed[e1] || removed[e2]) continue;
			for (const uint32_t e : {e1, e2}) {
				if (trn_[e] > k) {
				const uint32_t pe3 = bin[trn_[e]];
				const uint32_t pe = pos[e];
				if (pe3 != pe) {
					const uint32_t e3 = ord_[pe3];
					ord_[pe] = e3;
					pos[e3] = pe;
					ord_[pe3] = e;
					pos[e] = pe3;
				}
				++bin[trn_[e]];
				--trn_[e];
				}
			}
		}
  	}
	k_max = k;
}


void Mba::KdeltaTrussDecomp() {
    trussDecomp();
	vector<int> trussness(m_);
	memcpy(trussness.data(),trn_.data(),trussness.size() * sizeof(int));
    queue<uint32_t> q;
	vector<bool> Ins(m_,false);
    kspan_.resize(k_max + 1, vector<uint32_t>(tg_.m_, UINT32_MAX));
	for(uint32_t t = t_max ; t > 0 ; --t){
		if(time2triangle[t].size() == 0) continue;
		for(auto tri:time2triangle[t]){
			uint32_t e1 = std::get<0>(tri), e2 = std::get<1>(tri), e3 = std::get<2>(tri);
			Mc[e1][e3] = true;
			//update support
			s_[e1]--; s_[e2]--; s_[e3]--;
            uint32_t k1 = trussness[e1], k2 = trussness[e2], k3 = trussness[e3];
			//update ks
			if(k1<=k2 && k1<=k3) {
				if(ks_[e1]-- == k1) q.push(e1);
			} 
			if(k2<=k1 && k2<=k3) {
				if(ks_[e2]-- == k2) q.push(e2);
			} 
			if(k3<=k1 && k3<=k2) {
				if(ks_[e3]-- == k3) q.push(e3);
			}
			//update the trussness of edge
			while(!q.empty()){
				uint32_t e=q.front();
				q.pop();
                Ins[e] = false;
				uint32_t k = trussness[e];
				uint32_t supe = s_[e];
				if(supe == 0){
					kspan_[1][e] = t;
                    trussness[e] = 0;
					continue;
				}
				uint32_t nfound=0;
				uint32_t kse = 0;
                uint32_t v1 = tg_.edge_info_[e].first, v2 = tg_.edge_info_[e].second;
      			size_t p1 = 0, p2 = 0;
      			while (p1 < tg_.adj_[v1].size() && p2 < tg_.adj_[v2].size()) {
					if(nfound == supe) break;
					if (tg_.adj_[v1][p1].vid == tg_.adj_[v2][p2].vid) {
						getminmax(e, tg_.adj_[v1][p1].eid, tg_.adj_[v2][p2].eid);
						if(Mc[mine].find(maxe) != Mc[mine].end()) {
							++p1; ++p2;
							continue;
						}
						nfound++;
						uint32_t e1 = tg_.adj_[v1][p1].eid; uint32_t e2 = tg_.adj_[v2][p2].eid;
						++p1; ++p2;
						//update ks
						if(trussness[e1] < k-1 || trussness[e2] < k-1) continue;
						kse++;
						if(trussness[e1] >= k && trussness[e2] >= k) {
							if(trussness[e1] == k && Ins[e1] == false && ks_[e1]-- == trussness[e1]) {
								q.push(e1); Ins[e1] = true;
							}
							if(trussness[e2] == k && Ins[e2] == false && ks_[e2]-- == trussness[e2]) {
								q.push(e2); Ins[e2] = true;
							}
						}			
					}
					else if (tg_.adj_[v1][p1].vid < tg_.adj_[v2][p2].vid) {
						++p1;
					} else {
						++p2;
					}
      			}
				ks_[e] = kse;
				trussness[e]--;
				kspan_[k][e] = t;
			}
		}//for time2triangle[t]
	}
    for(uint32_t e = 0; e < m_; ++e){
		while(trussness[e] > 0){
			kspan_[trussness[e]][e] = 0;
			trussness[e]--; 
		}  
    }
}

vector<vector<int>> Mba::findkdCommunityForQuery(int query, int k, uint32_t delta) {
    unordered_set<int> visited;
    queue<int> q;
    vector<vector<int>> result;
    for (const auto& ae : tg_.adj_[query]) {
        int e = ae.eid;
		if (kspan_[k][e] <= delta && !visited.count(e)) {
			q.push(e); visited.insert(e); 
			vector<int> Ai;
			while (!q.empty()) {
            	int eid = q.front(); q.pop(); Ai.push_back(eid);
				vector<std::pair<uint32_t, uint32_t>> tris = tg_.GetTriangles(eid, k);
				for (const auto& tri : tris) {
					uint32_t e1 = tri.first, e2 = tri.second, deltatri = tg_.GetMst(tg_.tau_[e], tg_.tau_[e1], tg_.tau_[e2]);;
					if (deltatri > delta) continue;
					if (kspan_[k][e1] <= delta && !visited.count(e1)) {
						q.push(e1); visited.insert(e1);
					}
					if (kspan_[k][e2] <= delta && !visited.count(e2)) {
						q.push(e2); visited.insert(e2);
					}
				}
        	}
			result.emplace_back(std::move(Ai));
		}
    }
	for (auto& community : result) {
		for (int& e : community) {
			cout<<e<<" ";
		}
		cout<<endl;
	}
    return result;
 }

// void Mba::writeToFile(const std::string& filename) {
//     std::ofstream outfile(filename, std::ios::binary);
//     if (!outfile) {
//         std::cerr << "Failed to open file for writing: " << filename << std::endl;
//         return;
//     }
//     outfile.write(reinterpret_cast<const char*>(&n_), sizeof n_)
//          .write(reinterpret_cast<const char*>(&m_), sizeof m_);

//     for (uint32_t i = 0; i < m_; i++) {
//         // uint32_t tsize = edge2times[i].size();
//         // uint32_t header[3] = {edges_[i].first, edges_[i].second, tsize};
//         outfile.write(reinterpret_cast<const char*>(tg_.edge_info_[i].first), sizeof(tg_.edge_info_[i].first));
//         outfile.write(reinterpret_cast<const char*>(tg_.edge_info_[i].second), sizeof(tg_.edge_info_[i].second));
//         // outfile.write(reinterpret_cast<const char*>(edge2times[i].data()), tsize * sizeof(uint32_t));
//     }
//     outfile.write(reinterpret_cast<const char*>(&k_max), sizeof(k_max));
//     for (uint32_t k = 1; k <= k_max; ++k) {
//         for (auto& [eid, delta] : kspan_dict[k]) {
//             outfile.write(reinterpret_cast<const char*>(eid), sizeof(eid));
//             outfile.write(reinterpret_cast<const char*>(delta), sizeof(delta));
//         }
//     }

//     outfile.close();
//     std::cout << "Index successfully written to: " << filename << std::endl;
// }