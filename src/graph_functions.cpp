
#include "graph_functions.hpp"


namespace Graph {

  Graph::Graph(int n, double p) {
    random_connected_graph(n, *this, p);
  }

  bool eq(Edge& e, int n1, int n2) {
    return (e.V1 == n1 && e.V2 == n2)
      || (e.V1 == n2 && e.V2 == n1);
  }

  double* incidence(Graph& G) {
    double* ic = (double*) calloc(2*G.E*G.N, sizeof(double));
    for(int i=0; i<G.E; i++) {
      Edge& e = G.edges[i];
      ic[2*G.E*e.V1 + 2*i] = -1;
      ic[2*G.E*e.V2 + 2*i] = 1;
      ic[2*G.E*e.V1 + 2*i+1] = 1;
      ic[2*G.E*e.V2 + 2*i+1] = -1;
    }
    return ic;
  }

  int random_connected_graph(int N, Graph& G, double p)
  {

    int *found = (int*) calloc(N, sizeof(int));
    int *adj = (int*) calloc(N*N, sizeof(int));

    G.N = N; G.E=0; G.edges.clear();
    int node = Math::rng_int(0, N);
    int next, n_found=1; found[node]=1;

    while (n_found < N) {
      int next = Math::rng_int(0, N);
      if (!found[next]) {
        found[next]=1; n_found++;
        G.edges.push_back(Edge(node, next));
        adj[node*N+next] = 1; adj[next*N+node] = 1;
      }
      node = next;
    }
    G.E = N - 1;

    for (int i=0; i<N; i++) {
      for (int j=i+1; j<N; j++) {
        if (!adj[i*N+j] && Math::rng_uniform()<p) {
          G.edges.push_back(Edge(i, j)); G.E++;
        }
      }
    }
    free(found); free(adj); return 1;
  }

  int adj_list(Graph& G, std::vector<std::vector<int> >& lst)
  {
    for (int i=0; i < G.N; i++) {
      lst.push_back(std::vector<int>());
    }

    for(int i=0; i<G.E; i++) {
      int v1=G.edges[i].V1, v2 = G.edges[i].V2;
      lst[v1].push_back(v2);
      lst[v2].push_back(v1);
    }
  }

  int depth_first_search(std::vector<std::vector<int> >& lst,
    int n0, std::vector<int>& nodes)
  {

    int N = lst.size();
    std::vector<int> to_search;
    to_search.push_back(n0);

    int* visited = (int*) calloc(N, sizeof(int));
    nodes.clear(); visited[n0] = 1;

    while (to_search.size() > 0) {
      int k, node = to_search.back();
      to_search.pop_back();

      for (int i=0; i<lst[node].size(); i++) {
        k = lst[node][i];

        if (!visited[k]) {
          to_search.push_back(k);
          visited[k] = 1;
        }
      }
    }

    for (int i=0; i<N; i++) {
      if (visited[i]) nodes.push_back(i);
    }
    free(visited); return 1;

  }

  int depth_first_search(Graph& G, int n0, std::vector<int>& nodes)
  {
    std::vector<std::vector<int> > adjlst;
    adj_list(G, adjlst);
    depth_first_search(adjlst, n0, nodes);
  }

  int connect(Graph& G)
  {
    std::vector<std::vector<int> > adjlst;
    std::vector<std::vector<int> > clusters;

    std::vector<int> nodes; adj_list(G, adjlst);
    int* found = (int*) calloc(G.N, sizeof(int));

    for (int i=0; i<G.N; i++) {
      if (!found[i]) {
        depth_first_search(adjlst, i, nodes);
        clusters.push_back(nodes);
        for (int j=0; j<nodes.size(); j++) {
          found[nodes[j]] = 1;
        }
      }
    }

    // graph is already connected
    if (clusters.size() == 1) {
      free(found); return 1;
    }

    free(found); int n_clusters = clusters.size();
    int *c_found = (int*) calloc(G.N, sizeof(int));
    int node = Math::rng_int(0, n_clusters);
    int next, n_found=1; c_found[node] = 1;

    while (n_found < n_clusters) {
      next = Math::rng_int(0, n_clusters);
      if (!c_found[next]) {
        int r1 = Math::rng_int(0, clusters[node].size());
        int r2 = Math::rng_int(0, clusters[next].size());
        int v1=clusters[node][r1], v2=clusters[next][r2];
        G.edges.push_back(Edge(v1, v2));
        c_found[next] = 1; n_found++; G.E++;
      }
    }
    free(c_found); return 1;
  }

  int add_edge(Graph& G, bool force)
  {
    if (!force) {
      int r1 = Math::rng_int(0, G.N);
      int r2 = Math::rng_int(0, G.N-1);
      r2 = r2 >= r1 ? r2+1 : r2;

      for (int i = 0; i < G.E; i++) {
        if (eq(G.edges[i], r1, r2)) return 0;
      }
      G.edges.push_back(Edge(r1, r2));
      G.E++; return 1;
    }

    else {

      int N = G.N, *adj = (int*) calloc(N*N, sizeof(int));
      for (int i=0; i < G.E; i++) {
        int v1=G.edges[i].V1, v2=G.edges[i].V2;
        adj[N*v1+v2] = 1; adj[N*v2+v1] = 1;
      }

      std::vector<int> s1, s2;
      for (int i=0; i<G.N; i++) {
        for (int j=i+1; j<G.N; j++) {
          if (!adj[N*i + j]) {
            s1.push_back(i);
            s2.push_back(j);
          }
        }
      }
      if (s1.empty()) return 0;

      int rnd = Math::rng_int(0, s1.size());
      G.edges.push_back(Edge(s1[rnd], s2[rnd]));
      G.E++; free(adj); return 1;
    }
    return 0;
  }

  int add_node(Graph& G, bool force) {
    int rnd = Math::rng_int(0, G.N);
    G.edges.push_back(Edge(rnd, G.N));
    G.E++; G.N++; return 1;
  }

  int rm_edge(Graph& G, int* idx, bool reconnect) {

    // make sure there is an edge to remove
    if (G.E <= 1) {
      return 0;
    }

    int e0=G.E, rnd = Math::rng_int(0, G.E);
    G.edges[rnd] = G.edges.back();

    idx[rnd] = idx[e0-1];
    idx[e0-1] = -1;

    G.edges.pop_back(); G.E--;

    // reconnect graph if broken
    if (reconnect) {
      connect(G);
    }

    return 1;

  }

  int rm_node(Graph& G, int* nidx, int* eidx, bool reconnect) {

    // make sure there is an edge to remove
    if (G.N <= 1) {
      return 0;
    }


    int e0 = G.E, n0 = G.N;
    int rnd = Math::rng_int(0, G.N);

    int *t_eidx = (int*) malloc(G.E*sizeof(int));
    memcpy(t_eidx, eidx, e0*sizeof(int));

    for (int i=rnd; i<G.N-1; i++) {
      nidx[i] = nidx[i+1];
    }
    nidx[G.N-1] = -1; G.N--;

    int i, index = 0;
    std::vector<Edge> tmp;

    for (int i=0; i<e0; i++) {
      Edge& e = G.edges[i];
      if (!(e.V1 == rnd || e.V2 == rnd)) {

        if (e.V1 > rnd) e.V1--;
        if (e.V2 > rnd) e.V2--;

        tmp.push_back(e);
        eidx[index] = t_eidx[i];
        index++;
      }
    }
    G.edges = tmp; G.E = tmp.size();


    for (int i=index; i<e0; i++) {
      eidx[i] = -1;
    }

    // reconnect graph if broken
    if (reconnect) {
      connect(G);
    }

    free(t_eidx); return rnd;
  }

  std::ostream& operator<< (std::ostream& os, const Graph& G)
  {
    for (int i=0; i<G.E; i++) {
      const Edge& e = G.edges[i];
      os << "" << e.V1 << ", " << e.V2 << "\n";
    }
    return os;
  }

}
