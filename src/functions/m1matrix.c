// ================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file  m1matrix.c
//  \brief Routines for operating on M1Matrix structures

#include "functions.h"
#include "bns_nurates.h"

void InitializeM1MatrixSingleFlavor(M1Matrix *mat, const int n, const int idx) {
  mat->m1_mat_ab[idx] = (double **) malloc(sizeof(double *) * 2 * n);
  mat->m1_mat_em[idx] = (double **) malloc(sizeof(double *) * 2 * n);

  for (int i = 0; i < 2 * n; i++) {
    mat->m1_mat_ab[idx][i] = (double *) malloc(sizeof(double) * 2 * n);
    mat->m1_mat_em[idx][i] = (double *) malloc(sizeof(double) * 2 * n);

    for (int j = 0; j < 2 * n; j++) {
      mat->m1_mat_ab[idx][i][j] = 0.;
      mat->m1_mat_em[idx][i][j] = 0.;
    }
  }

  return;
}

void FreeM1MatrixSingleFlavor(M1Matrix *mat, const int n, const int idx) {
  for (int i = 0; i < 2 * n; i++) {
    free(mat->m1_mat_ab[idx][i]);
    free(mat->m1_mat_em[idx][i]);
  }

  free(mat->m1_mat_ab[idx]);
  free(mat->m1_mat_em[idx]);
   
  return;
}

void InitializeM1Matrix(M1Matrix *mat, const int n) {
  for (int idx = 0; idx < total_num_species; idx++) {
    InitializeM1MatrixSingleFlavor(mat, n, idx);
  }

  return;
}

void FreeM1Matrix(M1Matrix *mat, const int n) {
  for (int idx = 0; idx < total_num_species; idx++) {
    FreeM1MatrixSingleFlavor(mat, n, idx);
  }

  return;
}