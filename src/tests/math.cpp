#include "math_functions.hpp"
#include <time.h>
#include <stdio.h>

int main(int argc, char* argv[]) {

  Math::init_stream(time(NULL));

  for (int i = 0; i < 5; i++) {
    printf("%.4f\n", Math::rng_gaussian());
  }
  printf("\n");

  for (int i = 0; i < 5; i++) {
    printf("%.4f\n", Math::rng_uniform());
  }

  exit(0);
}
