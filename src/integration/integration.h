//
// Created by maitraya on 6/2/23.
//

#ifndef BNS_NURATES_SRC_INTEGRATION_INTEGRATION_H_
#define BNS_NURATES_SRC_INTEGRATION_INTEGRATION_H_
#include "../bns_nurates.h"

void GaussLegendre(MyQuadrature *quad);
void GaussLaguerre(MyQuadrature *quad);

void SaveQuadrature(char* filedir, MyQuadrature *quad);
void LoadQuadrature(char* filedir, MyQuadrature *quad);
#endif //BNS_NURATES_SRC_INTEGRATION_INTEGRATION_H_
