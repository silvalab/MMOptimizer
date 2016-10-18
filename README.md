# MarkovChannel
Ion channel Markov model parameter optimization

### Introduction

MarkovChannel is a tool for optimizing Markov models of ion channels, using the matrix exponential to greatly reduce fitting times.  It provides a flexable framework for encoding a range of experimentally collected protject into the model fitting.

### License

MarkovChannel is released under the MIT Licensse (refer to the LICENSE file for details).

### Contents
0. [Requirements: software](#requirements-software)
0. [Compilation](#compilation)
0. [Running the Demos](#running-the-demos)
0. [Understanding the Results](#understanding-the-results)
0. [Training on Other Data](#training-on-other-data)


### Requirements: software

0. 'Protobuf'
    - if you are using ubuntu run `sudo apt-get install libprotobuf-dev protobuf-compiler`
0. 'mkl' Intel Math Kernel Library
    - [mkl](https://software.intel.com/en-us/intel-mkl)
0. 'gsl' Gnu Scientific Library
    - [gsl](https://www.gnu.org/software/gsl/)

### Compilation

0. To compile simply run `make`

### Running the Demos

MarkovChannel comes with two demos; one for a Na<sup>+</sup> channel and one for a K<sup>+</sup> channel.

Before executing either of the examples, you must first set the mkl enviornment variables
* run `source /opt/intel/mkl/bin/mkvars.sh intel64`

To execute the Na<sup>+</sup> optimization
* run `./MarkovChannel solver.prototxt demos/Na+/protocols.lst`

To exectute the K<sup>+</sup> optimization
* run `./MarkovChannel solver.prototxt demos/K+/protocols.lst`

When running these commands, optimization progress will be periodically displayed.  More detailed information and fitted models will be written to the snapshot directory.

### Understanding the Results

Results of the demos are written to the snapshot directory in the "K+\_Demo" and "Na+\_Demo" respectivly.  In each of these subdirectories, you will find iter\_(%d).txt and iter\_(%d).model.  The .txt files contain the model fits at that iteration for each of the protocols.  The .model file provides the model structure and rate parameters that determine the behavior of the model.

Included in MarkovChannel are some MATLAB scripts that can interpret the .model files.


### Training on Other Data

To fit models on other data, you must encode the experimental protocols in the .prototxt format.

For example, the Na+ inactivation protocol is represented as:
```
name: "inac"
source: "inac.dat"
v0: -120.0
normalize: true

step {
  dt: 200.0
  stype: NONE
}

step {
  dt: 2.5
  vm: -20
  stype: PEAK
  stepsizze: 0.05
}
```
with inac.dat
```
-120.0000    1.0000
-110.0000    0.9796
-100.0000    0.9181
 -90.0000    0.7607
 -80.0000    0.4866
 -70.0000    0.2226
 -60.0000    0.0800
 -40.0000    0.0080
 -20.0000    0.0007
        0    0.0001
  20.0000         0
```

The header contains to following fields
* name - the name of the protocol
* source - the data used by the protocol
* v0 - the initial voltage
* normalize - whether to normalize output to [0, 1]

The protocol can then consist of any number of steps.  Each step contains the following fields:
* dt - the duration of the step
* vm - step voltage
* stype - the type of step this can be one of the following forms
	+ NONE - simply perform the step, produce no output
	+ PEAK - record the peak channel conductance
	+ MIN - record the minimum channel conductance
	+ TAU - record rate constants (parameterized by extra_args)
	+ TRACE - record channel conductance after each stepsize
* stepsize - optimal parameter for size of ODE/EXPM stepping
* extra_args - any number of additional step arguments, used for TAU stype

If either 'dt' or 'vm' is missing from the step, then the program treats this value as a variable and searches the .dat file for the missing paramters.  For example, in the inactivation protocol, the first step
```
step {
  dt: 200.0
  stype: NONE
}
```
is missing the 'vm' argument.  The .dat file is then searched and fills in 'vm' with the values of the first column in the .dat file.  So 'vm' becomes {-120, -110, -100, ..., 20}.

In the second step
```
step {
  dt: 2.5
  vm: -20
  stype: PEAK
  stepsizze: 0.05
}
```

stype is set to PEAK.  After the first step, vm is variable, and is set to all the values specified in the first row of the .dat file.  The second step then exposes each of these traces to -20mV for 2.5 ms and records the peak conductance.  The stepsize is 0.05, so a 5 microsecond resolution is used to calculate the peak.  The value of the peak conductance for each of these traces is then compared to the second column of the .dat file to compute the protocol error for the simulated model.

Another example protocol can be seen in rise.prototxt
```
name: "rise_wt"
source: "data/data/rise.txt"
v0: -120.0
normalize: false
weight: 2.0

step {
  dt: 1.6
  stype: TAU
  stepsize: 0.02
  extra_args: 0.1
  extra_args: 0.9
}
```

Here, the extra_args parameters are used.  In this protocol, the amount of time it take each trace to move from .1 to .9 of its peak conductance is measured.

More examples of protocol encodings can be found in the demos folder.


Once all protocols have been encoded in the proper format, then you must create a .lst file listing the location of all the desired protocols relative to the path of the .lst file.

Finally, you must specify a solver.prototxt providing the model/solver parameters
```
solver_mode: SIMULATED_ANNEALING
simulation_mode: ODE

protocol_list: "data/protocol_list.txt"

node_penality: 0.0004
edge_penality: 0.0000
eig_penality: 0.0001

model_param {

  n_prms: 3

  min_states: 3
  max_states: 10
  mu: 0
  std: 2

  mutation {
    add_edge: 0.18
    add_node: 0.05
    rm_edge: 0.05
    rm_node: 0.05

    update_std: 0.2
    update_prob: 0.2

    g_prob: 0.05
    f_prob: 0.05
  }

}

sa_param {

  k_max: 200000
  n_chains: 25
  step: 500
  gamma: 0.99
  t0: 0.0020
  display: 5000
  restart: 0.0001

  snapshot: 1000
  snapshotdir: "snapshots/t2"
```

The default value of the solver.prototxt works for a wide range of protocol settings.  However, some important paramters include:
* node_penality - penalize the model for each node
* edge_penality - penalize the model for each edge
* eig_peanlity - penalize model stiffness
* k_max - number of simulated annealing steps
* gamma - annealing schedule
* snapshotdir - the directory to write the optimized models (make sure this directory actually exists)

Finally, you can begin fitting the model by running `./MarkovChannel solver.prototoxt protocols.lst`


