
#ifndef CHANNEL_PROTOCOL_HPP_
#define CHANNEL_PROTOCOL_HPP_

#include "MarkovChannel.pb.h"
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <iostream>


enum StepType {NONE, PEAK, TAU, TRACE, MIN};
enum DataType {CONDUCTANCE, FLUORIMETRY};


struct Step {
  double dt;
  double vm;

  double stepsize;
  std::vector<double> args;

  StepType stype;
  DataType dtype;
};

struct ChannelProtocol {

  ChannelProtocol(std::string& prototxt);
  std::vector<std::vector<Step> > traces;
  MarkovChannel::ProtocolParameter params;

  std::vector<double> data, vars;
  double v0; int n_traces; bool norm;

};



std::ostream& operator<<(std::ostream& os, const Step& step);

std::ostream& operator<<(std::ostream& os, const std::vector<Step>& steps);

std::ostream& operator<<(std::ostream& os, const ChannelProtocol& proto);

#endif
