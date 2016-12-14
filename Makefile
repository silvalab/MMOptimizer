# SRC := src/*.cpp

PROTOC = protoc
CC = gcc
FORT = gfortran

MKLROOT = /opt/intel/mkl

EXTRA_INCLUDE_DIRS = -I/usr/local/include
EXTRA_LINK_LIBS = -L/usr/local/lib -L/opt/intel/lib/intel64

INCLUDE_DIRS = -I include -I proto -I$(MKLROOT)/include $(EXTRA_INCLUDE_DIRS)

CFLAGS = $(INCLUDE_DIRS) -m64 -D USE_MKL -fopenmp

LFLAGS1 = -L$(MKLROOT)/lib/intel64 $(EXTRA_LINK_LIBS)

LIBS = -lgfortran -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lpthread -lm -ldl -lprotobuf -liomp5 -lgsl -lgslcblas -lblas -llapack -lstdc++ -fopenmp

LFLAGS2 = $(LFLAGS1) $(LIBS)


MAIN = MMOptimizer

SRCS = \
	src/ChannelProtocol.cpp \
	src/math_functions.cpp \
	src/graph_functions.cpp \
	src/Model.cpp \
	src/cost.cpp \
	src/SimulatedAnnealing.cpp \
	src/main.cpp \
	src/private/dgpadm.f \

OBJS = \
        src/ChannelProtocol.o \
        src/math_functions.o \
        src/graph_functions.o \
        src/Model.o \
        src/cost.o \
        src/SimulatedAnnealing.o \
        src/main.o \
        src/private/dgpadm.o \


OBJECTS = $($(SRCS:.f=.o) $(SRCS:.cpp=.o) *.o)
PROTOS = proto/MarkovChannel.pb.h

all: $(MAIN)

$(MAIN): $(PROTOS) $(OBJS)
	@echo $(OBJECTS)
	$(CC) -o $@ $(CFLAGS) $(OBJS) proto/MarkovChannel.pb.cc  $(LFLAGS2)

.cpp.o: proto/MarkovChannel.pb.h
	$(CC) $(CFLAGS) -c $< -o $@

.f.o:
	$(FORT) $(INCLUDE_DIRS) -c $< -o $@

proto/MarkovChannel.pb.h: proto/MarkovChannel.proto
	$(PROTOC) -I=proto --cpp_out=proto proto/MarkovChannel.proto


clean:
	$(RM) src/*.o
	$(RM) proto/*.pb.h
	$(RM) proto/*.pb.cc
	$(RM) src/private/*.o


# DO NOT DELETE THIS LINE -- make depend needs it
