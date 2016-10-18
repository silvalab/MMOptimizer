#include "Model.hpp"
#include "cost.hpp"

#include "math_functions.hpp"
#include "graph_functions.hpp"
#include "ChannelProtocol.hpp"
#include "SimulatedAnnealing.hpp"

#include <fstream>
#include <iostream>
#include <vector>
#include <cstdlib>
#include <stdio.h>
#include <limits.h>
#include <math.h>
#include <time.h>

#include <fcntl.h>
#include <google/protobuf/io/coded_stream.h>
#include <google/protobuf/io/zero_copy_stream_impl.h>
#include <google/protobuf/text_format.h>

using namespace std;
using namespace MarkovChannel;

using google::protobuf::io::FileInputStream;
using google::protobuf::io::FileOutputStream;
using google::protobuf::io::ZeroCopyInputStream;
using google::protobuf::io::CodedInputStream;
using google::protobuf::io::ZeroCopyOutputStream;
using google::protobuf::io::CodedOutputStream;
using google::protobuf::Message;


int print_ic (Graph::Graph& G)
{
  double* ic = (double*) calloc(2*G.E*G.N, sizeof(double));
  for(int i=0; i<G.E; i++) {
    Graph::Edge& e = G.edges[i];
    ic[2*G.E*e.V1 + 2*i] = -1;
    ic[2*G.E*e.V2 + 2*i] = 1;
    ic[2*G.E*e.V1 + 2*i+1] = 1;
    ic[2*G.E*e.V2 + 2*i+1] = -1;
  }
  for (int i=0; i<G.N; i++) {
    for (int j=0; j<2*G.E; j++) {
      printf("%8.4f\t", ic[i*2*G.E+j]);
    }
    printf("\n");
  }
  free(ic);
}

void print_mat (double* A, int n, int m) {
  for (int i=0; i<n; i++) {
    for (int j=0; j<m; j++) {
      printf("%8.4f\t", A[i*m+j]);
    }
    printf("\n");
  }
}


MarkovChannel::SolverParameter solver_param;
vector<ChannelProtocol> protos;


double cost_f(Model::Model* m, ostream* os)
{
  return cost(*m, protos, solver_param, os);
}

void load_protocols(char *protolst)
{
  char datapath[PATH_MAX+1];
  realpath(protolst, datapath);
  std::ifstream protofile;
  protofile.open(protolst);

  std::string p1(datapath);
  int i = p1.size() - 1;
  while (p1[i] != '/') i--;
  std::string pth = p1.substr(0, i+1);

  std::string line;
  while ( protofile >> line ) {
    std::cout << pth+line << std::endl;
    std::string proto_file = pth + line;
    protos.push_back(ChannelProtocol(proto_file));
  }
}


int main(int argc, char* argv[])
{

  // google told me to do this
  GOOGLE_PROTOBUF_VERIFY_VERSION;
  Math::init_stream(time(NULL));

  int fd = open(argv[1], O_RDONLY);
  FileInputStream* input = new FileInputStream(fd);

  google::protobuf::TextFormat::Parse(input, &solver_param);
  Model::prms = solver_param.model_param(); delete input;
  load_protocols(argv[2]);

  vector<Model::Model*> models;
  MarkovChannel::SAParameter sa_param = solver_param.sa_param();
  SimulatedAnnealing::solve(cost_f, sa_param);

  // google told me to do this too
  google::protobuf::ShutdownProtobufLibrary();
  exit(0);

}
