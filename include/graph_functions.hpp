#ifndef GRAPH_FUNCTIONS_HPP_
#define GRAPH_FUNCTIONS_HPP_

#include <vector>
#include <string>
#include <cstdlib>
#include <cstring>
#include <iostream>

#include "mkl.h"
#include "math_functions.hpp"

namespace Graph {

  struct Edge {
    Edge(int v1, int v2) : V1(v1), V2(v2) {};
    Edge(const Edge& e) : V1(e.V1), V2(e.V2) {};
    int V1, V2;
  };

  struct Graph {
    Graph() : E(0), N(0) {};
    Graph(int n, double p=0.2);
    std::vector<Edge> edges;
    int E, N;
  };

  double* incidence(Graph& G);

  int random_connected_graph(int N, Graph& G, double p=0.2);

  int adj_list(Graph& G, std::vector<std::vector<int> >& lst);

  int depth_first_search(std::vector<std::vector<int> >& lst,
    int n0, std::vector<int>& nodes);

  int depth_first_search(Graph& G, int n0, std::vector<int>& nodes);

  int connect(Graph& G);

  int add_edge(Graph& G, bool force=false);

  int add_node(Graph& G, bool force=false);

  int rm_edge(Graph& G, int* idx, bool reconnect=true);

  int rm_node(Graph& G, int* nidx, int* eidx, bool reconnect=true);

  std::ostream& operator<< (std::ostream& os, const Graph& G);

}

#endif
