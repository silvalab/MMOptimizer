#ifndef SIMULATED_ANNEALING_HPP_
#define SIMULATED_ANNEALING_HPP_


#include "MarkovChannel.pb.h"
#include "Model.hpp"
#include "math_functions.hpp"

#include <vector>
#include <string>


typedef double (*cost_function)(Model::Model*, std::ostream*);

namespace SimulatedAnnealing
{
  void solve(cost_function cost, MarkovChannel::SAParameter& params);
}


#endif
