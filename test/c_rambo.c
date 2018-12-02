#include "c_rambo.h"
#include <stdio.h>

int main() {
  void * phs;
  const int n = 2;
  double Q[4] = {500, 0, 0, 0};
  double m[2] = {0, 0};
  phs = declare_phs_rambo (n, Q, m);
  write_phs_rambo (phs);

  double r[2] = {0.1, 0.9};
  generate_phs_rambo (phs, 2, r);
  double w = get_weight_phs_rambo (phs);

  printf("Weight: %f\n", w);

  write_phs_rambo (phs);

  free_phs_rambo (phs);
}
