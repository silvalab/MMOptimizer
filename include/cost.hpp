#ifndef COST_HPP_
#define COST_HPP_

#include <vector>
#include <iostream>
#include <cstdlib>
#include <mkl.h>

#include "Model.hpp"
#include "ChannelProtocol.hpp"
#include "math_functions.hpp"
#include "graph_functions.hpp"
#include "MarkovChannel.pb.h"

double cost(Model::Model& m, std::vector<ChannelProtocol>& protos,
   MarkovChannel::SolverParameter& sparam, std::ostream* os=NULL);

#endif
