#ifndef BNS_NURATES_SRC_INTEGRATION_INTEGRATION_H_
#define BNS_NURATES_SRC_INTEGRATION_INTEGRATION_H_

//=================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file  integration.h
//  \brief header file for all integration routines

#include "bns_nurates.h"

// routines for generating quadratures
void GaussLegendre(MyQuadrature *quad);
void GaussLegendreMultiD(MyQuadrature *quad);
void GaussLaguerre(MyQuadrature *quad);

// routines for saving and loading quadratures
void SaveQuadrature(char *filedir, MyQuadrature *quad);
void LoadQuadrature(char *filedir, MyQuadrature *quad);

// routines for integrating functions
double DoIntegration(const int n, const double *wtarray, const double *fnarray);
double GaussLegendreIntegrateZeroInf(MyQuadrature *quad, MyFunction *func, double t);
MyKernelQuantity GaussLegendreIntegrateZeroInfSpecial(MyQuadrature *quad, MyFunctionSpecial *func, double t);
double GaussLaguerreIntegrateZeroInf(MyQuadrature *quad, MyFunction *func);
MyQuadratureIntegrand GaussLegendreIntegrateFixedSplit2D(MyQuadrature *quad, MyFunctionMultiD *func, double t);
MyQuadratureIntegrand GaussLegendreIntegrateFixedSplit1D(MyQuadrature *quad, MyFunctionMultiD *func, double t);
MyQuadratureIntegrand GaussLegendreIntegrate2D(MyQuadrature *quad, MyFunctionMultiD *func, double *tx, double *ty);
MyQuadratureIntegrand GaussLegendreIntegrate1D(MyQuadrature *quad, MyFunctionMultiD *func, double *t);
MyQuadratureIntegrand GaussLegendreIntegrate2DMatrix(MyQuadrature *quad, M1Matrix *mat, double t);

#endif //BNS_NURATES_SRC_INTEGRATION_INTEGRATION_H_
