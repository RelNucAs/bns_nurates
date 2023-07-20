//
// Created by maitraya on 6/2/23.
//

#ifndef BNS_NURATES_TESTS_TESTS_H_
#define BNS_NURATES_TESTS_TESTS_H_


// tests with the quadrature generation and input/output
void PrintGaussLegendreQuadrature(const int n, const double x1, const double x2);
void TestGaussLegendreQuadrature();
void TestGaussLegendreQuadratureMultiD();
void PrintGaussLaguerreQuadrature(const int n, const double alpha);
void TestGaussLaguerreQuadrature();
void TestQuadratureInputOutput(char *filedir);

// tests with integration and comparison with GSL
void TestQuadratureWithGSL();
void TestIntegrationMultiD();

// tests for the pair process
void TestPairT();
void TestPairF();
void TestPairG();
void TestPairPsi();
void TestPairPhi();
void TestPairKernels();
void TestPairOpacities();

// tests for the bremsstrahlung process
void TestBremKernelS(char *filedir);
void TestBremKernelG(char *filedir);

#endif //BNS_NURATES_TESTS_TESTS_H_


