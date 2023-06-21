#ifndef BNS_NURATES_SRC_INTEGRATION_INTEGRATION_H_
#define BNS_NURATES_SRC_INTEGRATION_INTEGRATION_H_

//=================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file  integration.h
//  \brief header file for all integration routines

#include "../bns_nurates.h"

// routines for generating quadratures
void GaussLegendre(MyQuadrature *quad);
void GaussLaguerre(MyQuadrature *quad);

// routines for saving and loading quadratures
void SaveQuadrature(char *filedir, MyQuadrature *quad);
void LoadQuadrature(char *filedir, MyQuadrature *quad);

// routines for integrating functions
double GaussLegendreIntegrateZeroInf(MyQuadrature *quad, MyFunction *func, double t);
MyOpacity GaussLegendreIntegrateZeroInfSpecial(MyQuadrature *quad, MyFunction *func, double t);
double GaussLaguerreIntegrateZeroInf(MyQuadrature *quad, MyFunction *func);

#endif //BNS_NURATES_SRC_INTEGRATION_INTEGRATION_H_
