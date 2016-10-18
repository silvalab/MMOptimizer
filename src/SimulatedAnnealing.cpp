#include "SimulatedAnnealing.hpp"

#include <vector>
#include <string>
#include <math.h>
#include <time.h>

#include <string>
#include <sstream>
#include <fstream>

using namespace std;


namespace SimulatedAnnealing
{
  void solve(cost_function cost, MarkovChannel::SAParameter& params)
  {

    const int n_chains = params.n_chains();
    const int k_max = params.k_max();
    const int step = params.step();
    const int display = params.display();

    const int snapshot = params.snapshot();
    const string snapshotdir = params.snapshotdir();

    const double gamma = params.gamma();
    const double restart = params.restart();

    const double T0 = params.t0();
    double T = T0; // initialize temperature

    vector<Model::Model*> models(n_chains, NULL);
    vector<Model::Model*> argmin(n_chains, NULL);

    vector<double> f_val(n_chains, 1e12);
    vector<double> f_min(n_chains, 1e12);

    clock_t start = clock();

    #pragma omp parallel for
    for ( int i=0; i<n_chains; i++ ) {
      models[i] = new Model::Model(Math::rng_int(3, 10));
      argmin[i] = new Model::Model(models[i]);
      f_val[i] = cost(models[i], NULL); f_min[i] = f_val[i];
    }

    Model::Model* fmin_model = new Model::Model(models[0]);
    double fmin_val = f_val[0];

    for ( int i=0; i<n_chains; i++ ) {
      if ( f_val[i] < fmin_val ) {
        delete fmin_model; fmin_val = f_val[i];
        fmin_model = new Model::Model(models[i]);
      }
    }

    for ( int i=0; i<k_max; i++ ) {

      // update temperature
      if ( i>0 && i % step == 0) {
        T *= gamma;
      }

      #pragma omp parallel for
      for ( int j=0; j<n_chains; j++ ) {

        // generate a random neighbor and compute cost
        Model::Model* neighbor = Model::neighbor(models[j]);
        double fn = cost(neighbor, NULL);

        // use Metroplis-Hasting energy acceptance
        if ( Math::rng_uniform() < exp(-(fn - f_val[j])/T) ) {
          delete models[j]; f_val[j] = fn;
          models[j] = neighbor;
        }
        else {
          delete neighbor;
        }

        if ( f_val[j] < f_min[j] ) {
          delete argmin[j]; f_min[j] = f_val[j];
          argmin[j]= new Model::Model(models[j]);
        }

        // random restarts for stability and fun
        if ( Math::rng_uniform() < restart ) {
          int idx = Math::rng_int(0, n_chains);
          delete models[j]; f_val[j] = f_min[idx];
          models[j] = new Model::Model(argmin[idx]);
        }
      }

      for ( int j=0; j<n_chains; j++ ) {
        if ( f_val[j] < fmin_val ) {
          delete fmin_model; fmin_val = f_val[j];
          fmin_model = new Model::Model(models[j]);
        }
      }

      char buffer[100]; int nc;

      if ( i % params.display() == 0 ) {
        double time_elap = (double)(clock() - start)/CLOCKS_PER_SEC;
        nc = sprintf(buffer, "%8i\t%8.8f\t%8.2f", i, cost(fmin_model, NULL), time_elap);
        std::cout << std::string(buffer, nc) << std::endl;
      }

      if ( i % snapshot == 0 ) {
        ostringstream iss; string snapshot_file;
        iss << snapshotdir << "/iter_" << i << ".model";
        snapshot_file = iss.str();

        ostringstream mss; string modelfit_file;
        mss << snapshotdir << "/iter_" << i << ".txt";
        modelfit_file = mss.str();

        ofstream fss(snapshot_file.c_str());
        fss << *fmin_model << endl;
        fss.close();

        ofstream os(modelfit_file.c_str());
        cost(fmin_model, &os);
      }
    }

    delete fmin_model;
    for ( int i=0; i<n_chains; i++ ) {
      delete models[i];
      delete argmin[i];
    }
  }
}
