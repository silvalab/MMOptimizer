#ifndef MODEL_HPP_
#define MODEL_HPP_

#include <vector>
#include <string>
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <stdio.h>
#include <mkl.h>

#include "math_functions.hpp"
#include "graph_functions.hpp"
#include "MarkovChannel.pb.h"


namespace Model {

  extern MarkovChannel::ModelParameter prms;

  struct Model {

    static int count;
    int id;

    Model(int N, double p=0.2);
    Model(Model* m);
    Model();
    ~Model();

    Graph::Graph G;

    double *rs, *rk;
    double *C, *F;
    double *r_vec;

    int n_states() {return G.N;};
    int n_edges() {return G.E;};

    double loss;

  };

  Model* neighbor(Model* m, int num=1);

  double* initial_state(Model& m, double vm, double* s);

  double* transition_matrix(Model& m, double vm, double* Q);

  std::ostream& operator<<(std::ostream& os, const Model& m);

}

#endif
