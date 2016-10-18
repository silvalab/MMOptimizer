#include "../include/math_functions.hpp"

#define BRNG              VSL_BRNG_MCG31
#define METHOD_GAUSSIAN   VSL_RNG_METHOD_GAUSSIAN_ICDF
#define METHOD_UNIFORM    VSL_RNG_METHOD_UNIFORM_STD

namespace Math {

  VSLStreamStatePtr stream;

  int init_stream(const int seed) {
    return vslNewStream(&stream, BRNG, seed);
  }

  double rng_gaussian(double mu, double sig) {
    double r[1];
    vdRngGaussian(METHOD_GAUSSIAN, stream, 1, r, mu, sig);
    return *r;
  }

  double rng_uniform(double a, double b) {
    double r[1];
    vdRngUniform(METHOD_UNIFORM, stream, 1, r, a, b);
    return *r;
  }

  int rng_int(int a, int b) {
    int r[1];
    viRngUniform(METHOD_UNIFORM, stream, 1, r, a, b);
    return *r;
  }

  int rng_gaussian(int N, double *r, double mu, double sig) {
    return vdRngGaussian(METHOD_GAUSSIAN, stream, N, r, mu, sig);
  }

  int rng_uniform(int N, double *r, double a, double b) {
    return vdRngUniform(METHOD_UNIFORM, stream, N, r, a, b);
  }

  int rng_int(int N, int *r, int a, int b) {
    return viRngUniform(METHOD_UNIFORM, stream, N, r, a, b);
  }

}
