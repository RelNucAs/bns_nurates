//
// Created by maitraya on 6/2/23.
//

#ifndef BNS_NURATES_SRC_INTEGRATION_INTEGRATION_H_
#define BNS_NURATES_SRC_INTEGRATION_INTEGRATION_H_
#include "../bns_nurates.h"

void GaussLegendre(MyQuadrature *quad);
void GaussLaguerre(MyQuadrature quad, const double alpha);

void SaveQuadrature(const int n, const int dim, enum Quadrature type, const double x1, const double x2, const double y1, const double y2, const double alpha);
void LoadQuadrature(MyQuadrature quad, const int n, const int dim, enum Quadrature type, const double x1, const double x2, const double y1, const double y2, const double alpha);
#endif //BNS_NURATES_SRC_INTEGRATION_INTEGRATION_H_
