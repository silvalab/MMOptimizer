#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "MarkovChannel.pb.h"
#include "ChannelProtocol.hpp"

#include <fcntl.h>
#include <google/protobuf/io/coded_stream.h>
#include <google/protobuf/io/zero_copy_stream_impl.h>
#include <google/protobuf/text_format.h>


using namespace std;

int main(int argc, char* argv[])
{

  GOOGLE_PROTOBUF_VERIFY_VERSION;

  if (argc < 2) {
    cerr << "Usage: " << argv[0] << " CHANNEL_PROTOCOL " << endl;
    return -1;
  }

  string protobuf(argv[1]);
  ChannelProtocol proto(protobuf);

  for (int i = 0; i < proto.n_traces; i++ )
    cout << proto.traces[i] << endl;

  google::protobuf::ShutdownProtobufLibrary();

  return 0;
}
