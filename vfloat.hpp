#pragma once

#include "Fastor/Fastor.h"
#include <math.h>

using namespace Fastor;
using namespace std;

typedef FLT_T flt_t;

typedef Tensor<flt_t, nsample> flt_v;

vector<Tensor<bool, nsample>> flt_conds_vec;
bool have_true_cond = true;

#define CONCAT(a, b) a##b
#define EXPAND_CONCAT(a, b) CONCAT(a, b)
#define UNIQUE_NAME(base) EXPAND_CONCAT(base, __LINE__)

#ifndef SKIP_FALSE_BRANCHES
#define SKIP_FALSE_BRANCHES 1
#endif

#if SKIP_FALSE_BRANCHES == 1
#define vec_if(COND) {\
 if (any_of(COND)) {\
    auto flt_conds = flt_conds_vec.back();\
    flt_conds *= (COND)(flt_conds_vec.back());\
    have_true_cond = any_of(flt_conds);\
    if (have_true_cond)\
      flt_conds_vec.emplace_back(flt_conds);\
  } else {\
    have_true_cond = false;\
  }\
}\
if (have_true_cond)


#define vec_fi {}\
  if (have_true_cond)\
    flt_conds_vec.pop_back();\
  have_true_cond = true;


#else
#define vec_if(COND) { auto flt_conds = flt_conds_vec.back(); flt_conds *= COND; flt_conds_vec.emplace_back(flt_conds); }
#define vec_fi {} flt_conds_vec.pop_back();
#endif

#define flt_sqrt(X) evaluate(sqrt(X))

class vfloat {
public:
  flt_v v;
  vfloat() {};
  vfloat(const flt_v& v_) : v(v_) {};
  vfloat(const flt_t x) : v(x) {};
  flt_t operator[](const size_t idx) const { return v[idx]; }
  flt_t& operator[](const size_t idx) { return v[idx]; }
  vfloat operator+(const vfloat& other) const { return evaluate(v + other.v); }
  vfloat operator-(const vfloat& other) const { return evaluate(v - other.v); }
  vfloat operator*(const vfloat& other) const { return evaluate(v * other.v); }
  vfloat operator/(const vfloat& other) const { return evaluate(v / other.v); }

  vfloat operator+(const flt_t& other) const { return evaluate(v + other); }
  vfloat operator-(const flt_t& other) const { return evaluate(v - other); }
  vfloat operator*(const flt_t& other) const { return evaluate(v * other); }
  vfloat operator/(const flt_t& other) const { return evaluate(v / other); }

  Tensor<bool, nsample> operator<(const vfloat& other) const { return v < other.v; }
  Tensor<bool, nsample> operator<=(const vfloat& other) const { return v <= other.v; }
  Tensor<bool, nsample> operator>(const vfloat& other) const { return v > other.v; }
  Tensor<bool, nsample> operator>=(const vfloat& other) const { return v >= other.v; }

  Tensor<bool, nsample> operator<(const flt_t& other) const { return v < other; }
  Tensor<bool, nsample> operator<=(const flt_t& other) const { return v <= other; }
  Tensor<bool, nsample> operator>(const flt_t& other) const { return v > other; }
  Tensor<bool, nsample> operator>=(const flt_t& other) const { return v >= other; }

  vfloat operator+=(const flt_t& other) {
    v += flt_v(other)(flt_conds_vec.back());
    return *this;
  };

  vfloat operator+=(const int& other) {
    v += flt_v(other)(flt_conds_vec.back());
    return *this;
  };

  vfloat operator+=(vfloat&& other) {
    v += other.v(flt_conds_vec.back());
    return *this;
  };

  vfloat operator=(const vfloat& other) {
    flt_v& other_ = (flt_v&)other;
    v(flt_conds_vec.back()) = other_(flt_conds_vec.back());
    return *this;
  };
};

class vfloat2 {
public:
  vfloat v[2];
  vfloat2() {};
  vfloat2(const vfloat& v0, const vfloat& v1) : v{v0, v1} {};
  vfloat2(const flt_t x) : v{x, x} {};
  vfloat operator[](const size_t idx) const { return v[idx]; }
  vfloat& operator[](const size_t idx) { return v[idx]; }
  vfloat2 operator+(const vfloat2& other) const { return { v[0] + other.v[0], v[1] + other.v[1] }; }
  vfloat2 operator-(const vfloat2& other) const { return { v[0] - other.v[0], v[1] - other.v[1] }; }
  vfloat2 operator*(const vfloat2& other) const { return { v[0] * other.v[0], v[1] * other.v[1] }; }
  vfloat2 operator/(const vfloat2& other) const { return { v[0] / other.v[0], v[1] / other.v[1] }; }

  vfloat2 operator+(const vfloat& other) const { return { v[0] + other.v, v[1] + other.v }; }
  vfloat2 operator-(const vfloat& other) const { return { v[0] - other.v, v[1] - other.v }; }
  vfloat2 operator*(const vfloat& other) const { return { v[0] * other.v, v[1] * other.v }; }
  vfloat2 operator/(const vfloat& other) const { return { v[0] / other.v, v[1] / other.v }; }

  vfloat2 operator+(const flt_t& other) const { return { v[0] + other, v[1] + other }; }
  vfloat2 operator-(const flt_t& other) const { return { v[0] - other, v[1] - other }; }
  vfloat2 operator*(const flt_t& other) const { return { v[0] * other, v[1] * other }; }
  vfloat2 operator/(const flt_t& other) const { return { v[0] / other, v[1] / other }; }

  vfloat2& operator=(const vfloat2& other) {
    // fastor views seem to have trouble with the const above, but it's needed here
    flt_v& v0 = (flt_v&)other.v[0];
    flt_v& v1 = (flt_v&)other.v[1];
    v[0].v(flt_conds_vec.back()) = v0(flt_conds_vec.back());
    v[1].v(flt_conds_vec.back()) = v1(flt_conds_vec.back());
    return *this;
  }

  template<typename OTHER_T> vfloat2 operator+=(OTHER_T&& other) {
    v[0].v += other.v[0].v(flt_conds_vec.back());
    v[1].v += other.v[1].v(flt_conds_vec.back());
    return *this;
  };
  template<typename OTHER_T> vfloat2 operator-=(const OTHER_T& other) { *this = *this - other; return *this; };
  template<typename OTHER_T> vfloat2 operator*=(const OTHER_T& other) { *this = *this * other; return *this; };
  template<typename OTHER_T> vfloat2 operator/=(const OTHER_T& other) { *this = *this / other; return *this; };
};

class gvfloat {
public:
  vfloat v[nent];
  gvfloat() {};
  gvfloat(const flt_t x) { for (int e = 0; e < nent; e++) v[e] = x; };
  gvfloat(const gvfloat& other) { *this = other; };
  vfloat& operator[](const size_t idx) { return v[idx]; }
  gvfloat operator+(const vfloat& other) const { gvfloat r; for (int e = 0; e < nent; e++) r.v[e] = v[e] + other.v[e]; return r; }
  gvfloat operator-(const vfloat& other) const { gvfloat r; for (int e = 0; e < nent; e++) r.v[e] = v[e] - other.v[e]; return r; }
  gvfloat operator*(const vfloat& other) const { gvfloat r; for (int e = 0; e < nent; e++) r.v[e] = v[e] * other.v[e]; return r; }
  gvfloat operator/(const vfloat& other) const { gvfloat r; for (int e = 0; e < nent; e++) r.v[e] = v[e] / other.v[e]; return r; }
};

class gvfloat2 {
public:
  vfloat2 v[nent];
  gvfloat2() {};
  gvfloat2(const flt_t x) { for (int e = 0; e < nent; e++) v[e] = x; }
  vfloat2& operator[](const size_t idx) { return v[idx]; }
  gvfloat2 operator+(const gvfloat2& other) const { gvfloat2 r; for (int e = 0; e < nent; e++) r.v[e] = v[e] + other.v[e]; return r; }
  gvfloat2 operator-(const gvfloat2& other) const { gvfloat2 r; for (int e = 0; e < nent; e++) r.v[e] = v[e] - other.v[e]; return r; }
  gvfloat2 operator*(const gvfloat2& other) const { gvfloat2 r; for (int e = 0; e < nent; e++) r.v[e] = v[e] * other.v[e]; return r; }
  gvfloat2 operator/(const gvfloat2& other) const { gvfloat2 r; for (int e = 0; e < nent; e++) r.v[e] = v[e] / other.v[e]; return r; }
};

vfloat norm(vfloat2&& v) { return flt_sqrt(v[0].v * v[0].v + v[1].v * v[1].v); }
vfloat norm(vfloat2& v) { return flt_sqrt(v[0].v * v[0].v + v[1].v * v[1].v); }

vfloat dot(vfloat2& v0, vfloat2& v1) { return vfloat(v0[0].v * v1[0].v + v0[1].v * v1[1].v); }

vfloat operator+(const flt_t& lhs, const vfloat& rhs) { return evaluate(lhs + rhs.v); };
vfloat operator-(const flt_t& lhs, const vfloat& rhs) { return evaluate(lhs - rhs.v); };
vfloat operator*(const flt_t& lhs, const vfloat& rhs) { return evaluate(lhs * rhs.v); };
vfloat operator/(const flt_t& lhs, const vfloat& rhs) { return evaluate(lhs / rhs.v); };

vfloat2 operator+(const flt_t& lhs, const vfloat2& rhs) { return { lhs + rhs[0], lhs + rhs[1] }; }
vfloat2 operator-(const flt_t& lhs, const vfloat2& rhs) { return { lhs - rhs[0], lhs - rhs[1] }; }
vfloat2 operator*(const flt_t& lhs, const vfloat2& rhs) { return { lhs * rhs[0], lhs * rhs[1] }; }
vfloat2 operator/(const flt_t& lhs, const vfloat2& rhs) { return { lhs / rhs[0], lhs / rhs[1] }; }

flt_t sum(vfloat&& v) { flt_t r = 0.0; for (int rep = 0; rep < nsample; rep++) r += v[rep]; return r; }
flt_t sum(vfloat& v) { flt_t r = 0.0; for (int rep = 0; rep < nsample; rep++) r += v[rep]; return r; }

flt_t cond_density() {
  int c = 0;
  auto &conds = flt_conds_vec.back();
  for (int rep = 0; rep < nsample; rep++)
    c += conds[rep];
  return (double)c / nsample;
}

vfloat exp(const vfloat& v) { return vfloat(exp(v.v)); }
vfloat abs(const vfloat& v) { return vfloat(abs(v.v)); }
vfloat atan2(const vfloat& x, const vfloat& y) { return vfloat(atan2(x.v, y.v)); }
