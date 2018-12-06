#include "c_rambo.h"
#include <stdio.h>
#include <stdlib.h>

#define N_PRT 2
#define MASSES {0, 0}

void rng_init (const unsigned int seed);
double rng_generate ();
void rng_generate_array (const size_t n, double r[]);

void print_array_double (const size_t n, double r[]);

int main() {
  phs_rambo phs;

  double Q[4] = {500, 0, 0, 0};
  double m[N_PRT] = MASSES;
  double r[3 * N_PRT - 4];
  double p[4];

  phs = declare_phs_rambo (N_PRT, Q, m);

  rng_init (1961991);
  rng_generate_array (2, r);

  generate_phs_rambo (phs, 2, r);
  double w = get_weight_phs_rambo (phs);
  printf("Weight: %f\n", w);
  print_array_double (3 * N_PRT - 4, r);

  invert_phs_rambo (phs, r);
  print_array_double (3 * N_PRT - 4, r);

  for (unsigned int i = 0; i < N_PRT; i++) {
    get_p_phs_rambo (phs, i, p);
    print_array_double (4, p);
  }

  write_phs_rambo (phs);
  free_phs_rambo (phs);
}

// Wrapper function for (s)random from stdlib.h.
// We want to use double random numbers ∈ (0, 1).

void rng_init (const unsigned int seed) {
  srandom (seed);
}

double rng_generate () {
  return (double) random () / (double) RAND_MAX;
}

void rng_generate_array (const size_t n, double r[]) {
  for(size_t i = 0; i < n; i++) {
    r[i] = rng_generate ();
  }
}

void print_array_double (const size_t n, double r[]) {
  for (size_t i = 0; i < n; i++) {
    printf(" %f", r[i]);
  }
  printf("\n");
}
