#include <math.h>
#include <unordered_set>
#include <vector>
#include <array>
#include <assert.h>
#include <stdio.h>
#include <cfloat>
#include <random>
#include <chrono>
#include "params.hpp"
#include "vec.hpp"

const int num_iter = 100000000;
default_random_engine rng;
uniform_real_distribution<flt_t> u;
flt_t r;

flt init_a, init_b;

#define CALC (atan2(a, b))

flt_t run_scalar()
{
  auto t0 = chrono::high_resolution_clock::now();

  for (int sample = 0; sample < nsample; sample++) {
    flt_t a = init_a[sample], b = init_b[sample];
    for (int iter = 0; iter < num_iter; iter++) {
      a = CALC;
    }
    r += a;
  }

  auto t1 = chrono::high_resolution_clock::now();
  flt_t scalar_s = 1e-6 * chrono::duration_cast<chrono::microseconds>(t1 - t0).count();

  return scalar_s;
}

flt_t run_flt()
{
  flt a = init_a, b = init_b;

  auto t0 = chrono::high_resolution_clock::now();

  for (int iter = 0; iter < num_iter; iter++)
    a = CALC;

  for (int sample = 0; sample < nsample; sample++)
    r += a[sample];

  auto t1 = chrono::high_resolution_clock::now();
  flt_t flt_s = 1e-6 * chrono::duration_cast<chrono::microseconds>(t1 - t0).count();

  return flt_s;
}

flt_t run_fastor()
{
  Tensor<flt_t, nsample> a = init_a.v, b = init_b.v;

  auto t0 = chrono::high_resolution_clock::now();

  for (int iter = 0; iter < num_iter; iter++)
    a = CALC;

  for (int sample = 0; sample < nsample; sample++)
    r += a[sample];

  auto t1 = chrono::high_resolution_clock::now();
  flt_t fastor_s = 1e-6 * chrono::duration_cast<chrono::microseconds>(t1 - t0).count();

  return fastor_s;
}

int main(int argc, char **argv)
{
  flt_conds_vec.push_back(true);

  for (int sample = 0; sample < nsample; sample++) {
    init_a[sample] = u(rng);
    init_b[sample] = u(rng);
  }

  flt_t scalar_s = run_scalar();
  flt_t flt_s = run_flt();
  flt_t fastor_s = run_fastor();

  printf("scalar: %.4fs, 1.0\n", scalar_s);
  printf("flt: %.4fs, %.2f\n", flt_s, scalar_s / flt_s);
  printf("fastor: %.4fs, %.2f\n", fastor_s, scalar_s / fastor_s);

  return r;
}
