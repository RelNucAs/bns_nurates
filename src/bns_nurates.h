//
// Created by maitraya on 6/2/23.
//


#ifndef BNS_NURATES_SRC_BNS_NURATES_H_
#define BNS_NURATES_SRC_BNS_NURATES_H_

// ================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file bns_nurates.h
//  \brief essential data structures for the library

// -------------------------------------------------
// Quadrature specific data structures

// Quadrature enum
// Holds the type of quadrature
enum Quadrature { kGauleg, kGaulag };
typedef enum Quadrature Quadrature;

// MyQuadrature struct
// Contains quadrature information, supports 1d/2d

/* Members:
 *
 * type:  type of quadrature
 * dim:   dimension of quadrature (1d/2d)
 * n:     number of points in scheme
 * x1:    lower end point in x
 * x2:    upper end point in x
 * y1:    lower endpoint in y (optional)
 * y2:    upper endpoint in y (optional)
 * x:     points for quadrature scheme (x)
 * y:     points for quadrature scheme (y, optional)
 * w:     weights for quadrature scheme
 */
struct MyQuadrature {
  enum Quadrature type;
  int n;
  int dim;
  double alpha;
  double x1;
  double x2;
  double y1;
  double y2;
  double *x;
  double *y;
  double *w;
};
typedef struct MyQuadrature MyQuadrature;

struct MyFunction {
  int dim;
  double (*function)(double var, void *params);
  void *params;
};
typedef struct MyFunction MyFunction;

struct MyKernel {
  double absorption_e;
  double production_e;
  double absorption_x;
  double production_x;
};
typedef struct MyKernel MyKernel;

struct MyParams {
  double omega;
  double omega_prime;
  double temp;
  double nb;
  double m_N;
  int lmax;
  double filt;
  double eta;
  double cos_theta;
};
typedef struct MyParams MyParams;

// @FIXME: decide with Maitraya what to use for this
// Temporary struct for storing output of emissivity/absorptivity functions
struct MyOpacity {
  double em;
  double ab;
};
typedef struct MyOpacity MyOpacity;


#endif //BNS_NURATES_SRC_BNS_NURATES_H_
