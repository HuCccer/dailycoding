#include <unistd.h>
#include <chrono>
#include <cinttypes>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>

#include "mba_decomp.h"
#include "kdequiforests.h"
// #include "mba_decomp.cpp"
// #include "KDEquiForests.cpp"

// a sample program
int main(int argc, char** argv) {
  const std::string graph_file = "D:\\file\\program\\KDEquiTruss\\data\\superuser.txt";
  // freopen(".\\edges.txt","w",stdout);  
  std::ifstream infile(graph_file);
  uint32_t t = 0, n = 0, m = 0; // # of timestamps, vertices and temporal edges
  infile >> t >> n >> m;
  infile.close();
  // read the graph and the index
  // We set the param "l" to "m * 2" here, where "l" is the maximum number of
  // edges that a graph can hold. We impose this constraint only for ease of
  // implementation. TODO: remove this constraint in the future.
//   __debugbreak();
  Mba mba(n, m * 2, graph_file);
  auto beg = std::chrono::steady_clock::now();
  mba.KdeltaTrussDecomp();
  // for(int i = 0; i < mba.tg_.edge_info_.size(); i++) {
  //   cout<<i<<": ("<< mba.tg_.edge_info_[i].first<<","<<mba.tg_.edge_info_[i].second<<")"<<endl;
  // }
  auto end = std::chrono::steady_clock::now();
  auto dif = end - beg;
  printf("K-delta truss decomposition costs ");
  printf("%.3fms\n", std::chrono::duration<double, std::milli>(dif).count());
  KDEquiForests KDForsts(mba.kmax() + 1, mba.tg_.m(), mba.tg_);
  map<uint32_t, vector<int>> kspan_dict;
  beg = std::chrono::steady_clock::now();
  for (int k = 1; k <= mba.kmax(); k++) {
      mba.get_KSpan_dict(k, kspan_dict);
      KDForsts.constructIndexForK(kspan_dict, mba.kspan_[k], mba.tg_, KDForsts.KDEForests[k], k);
  }
  end = std::chrono::steady_clock::now();
  dif = end - beg;
  printf("KDEquiForest construction costs ");
  printf("%.3fms\n", std::chrono::duration<double, std::milli>(dif).count());
  // beg = std::chrono::steady_clock::now();
  // KDForsts.findkdCommunityForQuery(3, 2, 8);
  // end = std::chrono::steady_clock::now();
  // dif = end - beg;
  // printf("KDEquiForest-community query costs ");
  // printf("%.3fms\n", std::chrono::duration<double, std::milli>(dif * 10000).count());
  // beg = std::chrono::steady_clock::now();
  // mba.findkdCommunityForQuery(3, 2, 8);
  // end = std::chrono::steady_clock::now();
  // dif = end - beg;
  // printf("KDEquiForest-community query costs ");
  // printf("%.3fms\n", std::chrono::duration<double, std::milli>(dif * 10000).count());

  return 0;
}

