#include "cost.hpp"
// extern "C" {
//   #include "intel_ode.h"
// }
#include <math.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>

#include <sstream>

using namespace MarkovChannel;
using namespace std;

extern "C" {
  void dgpadm_(int* ideg, int* m, double* t, double* H, int* ldh,
    double* wsp, int* lwsp, int* ipiv, int* iexp, int* ns, int* iflag);
}


inline double peak(int n, double* x, bool mini=false)
{
  int idx = 0;
  if (mini) {
    idx = cblas_idamin(n, x, 1);
  } else {
    idx = cblas_idamax(n, x, 1);
  }
  double res = x[idx]; return res;
}

inline double tau(int n, double* x, double v1, double v2)
{
  int idx = cblas_idamax(n, x, 1);
  cblas_dscal(n, 1/x[idx], x, 1);
  const double eps=1e-6; double res;

  if (v2 > v1) {
    double i1=-eps, i2=idx+eps;

    for(int t=0; t < min(idx, n-1); t++) {
      if ( i1 <= 0 && x[t+1] > v1 && x[t] <= v1) {
        i1 = t + (v1 - x[t]) / (x[t+1] - x[t]);
      }
      if ( i2 >= idx && x[t+1] > v2 && x[t] <= v2) {
        i2 = t + (v2 - x[t]) / (x[t+1] - x[t]);
      }
    }
    res = i2 - i1;
  }

  else {
    double i1=idx-eps, i2=n+eps;

    for (int t=idx; t<(n-1); t++) {
      if (i1 <= idx && x[t+1] < v1 && x[t] >= v1) {
        i1 = t + (x[t] - v1) / (x[t+1] - x[t]);
      }
      if (i2 >= n && x[t+1] < v2 && x[t] >= v2) {
        i2 = t + (x[t] - v2) / (x[t+1] - x[t]);
      }
    }
    if ( i2 >= n ) {
      res = (i2 - i1) + 10*(x[n-1]-v2);
    }
    else { res = i2 - i1; }
  }

  return res;
}



int func(double t, const double *y, double *f, void *params)
{
  void** prms = (void **) params;
  int N = *(int *) prms[0];
  double *Q = *(double **) prms[1];

  cblas_dgemv(CblasRowMajor, CblasNoTrans,
    N, N, 1.0, Q, N, y, 1, 0.0, f, 1);

  return GSL_SUCCESS;
}

int jac(double t, const double *y, double *dfdy, double *dfdt, void *params)
{
  void** prms = (void **) params;
  int N = *(int *) prms[0];
  double *Q = *(double **) prms[1];

  memcpy(dfdy, Q, N*N*sizeof(double));
  return GSL_SUCCESS;
}

int step_ode(Model::Model& m, vector<Step>& steps,
  double *y, vector<double>& out)
{
  int N = m.n_states();
  double* Q = &y[N];
  void* params[2] = {(void*) &N, (void*) &Q};

  gsl_odeiv2_system sys = {func, jac, N, &params};
  gsl_odeiv2_driver * d =
    gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rkf45,
				  1e-5, 1e-6, 0.0);
  gsl_odeiv2_driver_set_hmin(d, 1e-6);

  for ( int i=0; i<steps.size(); i++ ) {

    double dt=steps[i].dt, vm=steps[i].vm;
    Model::transition_matrix(m, vm, Q);

    if ( steps[i].stype == NONE ) {
      double t=0, ti = dt;
      int status = gsl_odeiv2_driver_apply (d, &t, ti, y);
      if ( status != GSL_SUCCESS ) {
        gsl_odeiv2_driver_free (d); return -1;
      }
    }

    else {

      int n_steps = ceil(dt / steps[i].stepsize);
      double *c_mat = (steps[i].dtype == CONDUCTANCE) ? m.C : m.F;

      double *vals = (double*) malloc (n_steps*sizeof(double));
      vals[0] = cblas_ddot(N, y, 1, c_mat, 1);

      for ( int j=1; j<n_steps; j++ ) {

        double t=0, ti = steps[i].stepsize;
        int status = gsl_odeiv2_driver_apply (d, &t, ti, y);
        if ( status != GSL_SUCCESS ) {
          gsl_odeiv2_driver_free (d);
          free(vals); return -1;
        }
        vals[j] = cblas_ddot(N, y, 1, c_mat, 1);
      }

      if ( steps[i].stype == TAU ) {
        double arg1=steps[i].args[0], arg2=steps[i].args[1];
        double res = tau(n_steps, vals, arg1, arg2);
        out.push_back(res*steps[i].stepsize);
      }

      else if ( steps[i].stype == PEAK ) {
        double res = peak(n_steps, vals);
        out.push_back(res);
      }

      else if ( steps[i].stype == MIN ) {
        double res = peak(n_steps, vals, true);
        out.push_back(res);
      }

      else if ( steps[i].stype == TRACE ) {
        for ( int j=0; j<n_steps; j++ ) {
          out.push_back(vals[j]);
        }
      }
      free(vals);
    }
  }
  gsl_odeiv2_driver_free (d);
  return 1;
}


int step_exp(Model::Model& m, vector<Step>& steps,
  double *y, vector<double>& out)
{
  int N = m.n_states();
  double *H = &y[N], t;

  int iflag, iexp, ns, ideg = 6;

  int lwsp = 4*N*N + ideg + 1;
  double *wsp = (double*) malloc(lwsp*sizeof(double));
  int *ipiv = (int*) malloc(N*sizeof(int));

  for ( int i=0; i<steps.size(); i++ ) {

    double dt=steps[i].dt, vm=steps[i].vm;
    Model::transition_matrix(m, vm, H);

    mkl_dimatcopy('R', 'T', N, N, 1.0, H, N, N);

    double *y0 = (double*) malloc(N*sizeof(double));
    memcpy(y0, y, N*sizeof(double));

    if ( steps[i].stype == NONE ) {
      t = dt;
      dgpadm_(&ideg, &N, &t, H, &N, wsp, &lwsp, ipiv, &iexp, &ns, &iflag);

      if ( iflag < 0 ) {
        free(wsp); free(ipiv); free(y0); return -1;
      }

      cblas_dgemv(CblasColMajor, CblasNoTrans,
        N, N, 1.0, &wsp[iexp-1], N, y0, 1, 0.0, y, 1);
    }

    else {

      int n_steps = ceil(dt / steps[i].stepsize);
      double *c_mat = (steps[i].dtype == CONDUCTANCE) ? m.C : m.F;


      t = steps[i].stepsize;
      dgpadm_(&ideg, &N, &t, H, &N, wsp, &lwsp, ipiv, &iexp, &ns, &iflag);

      if ( iflag < 0 ) {
        free(wsp); free(ipiv); free(y0); return -1;
      }


      double *vals = (double*) malloc (n_steps*sizeof(double));
      vals[0] = cblas_ddot(N, y, 1, c_mat, 1);

      for ( int j=1; j<n_steps; j++ ) {
        cblas_dgemv(CblasColMajor, CblasNoTrans,
          N, N, 1.0, &wsp[iexp-1], N, y0, 1, 0.0, y, 1);

        vals[j] = cblas_ddot(N, y, 1, c_mat, 1);
        memcpy(y0, y, N*sizeof(double));
      }

      if ( steps[i].stype == TAU ) {
        double arg1=steps[i].args[0], arg2=steps[i].args[1];
        double res = tau(n_steps, vals, arg1, arg2);
        out.push_back(res*steps[i].stepsize);
      }

      else if ( steps[i].stype == PEAK ) {
        double res = peak(n_steps, vals);
        out.push_back(res);
      }

      else if ( steps[i].stype == MIN ) {
        double res = peak(n_steps, vals, true);
        out.push_back(res);
      }

      else if ( steps[i].stype == TRACE ) {
        for ( int j=0; j<n_steps; j++ ) {
          out.push_back(vals[j]);
        }
      }
      free(vals);
    }
    free(y0);
  }
  free(wsp);
  free(ipiv);
}


double cost(Model::Model& m, ChannelProtocol& proto,
  SolverParameter sparam, ostream* os = NULL)
{
  std::vector<double> output;
  int ierr, N = m.n_states();

  double *y0 = (double*) malloc(N*sizeof(double));
  double *y  = (double*) malloc((N+N*N)*sizeof(double));
  Model::initial_state(m, proto.v0, y0);


  for ( int i=0; i<proto.n_traces; i++ ) {
    memcpy(y, y0, N*sizeof(double));

    if ( sparam.simulation_mode() == SolverParameter::ODE ) {
      ierr = step_ode(m, proto.traces[i], y, output);
      if ( ierr < 0 ) {
        free(y); free(y0); return 1e6;
      }
    }
    else {
      ierr= step_exp(m, proto.traces[i], y, output);
      if ( ierr < 0 ) {
        free(y); free(y0); return 1e6;
      }
    }
  }

  int idx = 0, n = output.size();
  double dx, err;

  if ( proto.params.normalize() ) {
    double scale = 0;
    for (int i = 0; i < output.size(); i++) {
      scale = max(scale, output[i]);
    }
    for (int i = 0; i < output.size(); i++) {
      output[i] /= scale;
    }
  }

  if ( os ) {
    *os << proto.params.name() << endl; char buffer[32]; int nc;
    for ( int i=0; i<output.size(); i++ ) {
      nc = sprintf(buffer, "%8.4f\t%8.4f\n", output[i], proto.data[i]);
      *os << std::string(buffer, nc);
    }
    *os << endl;
  }


  for (int i = 0; i < output.size(); i++) {
    dx = output[i] - proto.data[i];
    err += (1.0 / n) * dx * dx;
  }

  free(y0); free(y);

  if (isnan(err))
    return 1e6;

  return err;
}


double model_penality(Model::Model& m, SolverParameter& sparam)
{
  double rcond, anorm, penality = 0;
  int N = m.n_states(), E = m.n_edges();

  penality += sparam.node_penality() * N;
  penality += sparam.edge_penality() * E;

  double *Q = (double*) malloc(N*N*sizeof(double));
  double *wr = (double*) malloc(N*sizeof(double));
  double *wi = (double*) malloc(N*sizeof(double));

  int idmax, idmin;

  for (double vm = -120; vm <= 40; vm += 20) {

    Model::transition_matrix(m, vm, Q);

    LAPACKE_dgeev( LAPACK_ROW_MAJOR,
      'N', 'N', N, Q, N, wr, wi, NULL, N, NULL, N );

    idmax = cblas_idamax(N, wr, 1);
    idmin = cblas_idamin(N, wr, 1);

    penality += sparam.eig_penality() * ((wr[idmax]<0) ? -wr[idmax] : wr[idmax]);

  }

  free(Q); free(wr); free(wi);

  return penality;
}


double cost(Model::Model& m, vector<ChannelProtocol>& protos,
  SolverParameter& sparam, ostream* os)
{
  double error = 0;
  for ( int i = 0; i < protos.size(); i++ ) {
    error += protos[i].params.weight() * cost(m, protos[i], sparam, os);
  }
  return  error + model_penality(m, sparam);
}



// void rhs_mat(int *n, double *t, double *y, double *f) {
//   cblas_dgemv(CblasRowMajor, CblasNoTrans,
//     *n, *n, 1.0, Q, *n, y, 1, 0.0, f, 1);
// }
// void jac_mat(int *n, double *t, double *y, double *a) {
//   mkl_domatcopy('R', 'T', *n, *n, 0.0, Q, *n, a, *n);
// }
