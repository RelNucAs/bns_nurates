#include <math.h>

#include "functions.h"

/*
Evaluate modified Bessel functions of order v = 1
Description:                
              The differential equation 

                       2
                   2  d w       dw      2   2
                  x . --- + x . --- - (x + v ).w = 0
                        2       dx
                      dx
                      
              has two solutions called modified Bessel functions
              Iv(x) and Kv(x).
              The routines bessi1 and bessk1 return the I and K for 
              v = 1.
*/

//static
double bessi1(const double x)
/*------------------------------------------------------------*/
/* PURPOSE: Evaluate modified Bessel function In(x) and n=1.  */
/*------------------------------------------------------------*/
{
   double ax,ans;
   double y;


   if ((ax=fabs(x)) < 3.75) {
      y=x/3.75,y=y*y;
      ans=ax*(0.5+y*(0.87890594+y*(0.51498869+y*(0.15084934
         +y*(0.2658733e-1+y*(0.301532e-2+y*0.32411e-3))))));
   } else {
      y=3.75/ax;
      ans=0.2282967e-1+y*(-0.2895312e-1+y*(0.1787654e-1
         -y*0.420059e-2));
      ans=0.39894228+y*(-0.3988024e-1+y*(-0.362018e-2
         +y*(0.163801e-2+y*(-0.1031555e-1+y*ans))));
      ans *= (exp(ax)/sqrt(ax));
   }
   return x < 0.0 ? -ans : ans;
}


//static
double bessk1(const double x)
/*------------------------------------------------------------*/
/* PURPOSE: Evaluate modified Bessel function Kn(x) and n=1.  */
/*------------------------------------------------------------*/
{
   double y,ans;

   if (x <= 2.0) {
      y=x*x/4.0;
      ans=(log(x/2.0)*bessi1(x))+(1.0/x)*(1.0+y*(0.15443144
         +y*(-0.67278579+y*(-0.18156897+y*(-0.1919402e-1
         +y*(-0.110404e-2+y*(-0.4686e-4)))))));
   } else {
      y=2.0/x;
      ans=(exp(-x)/sqrt(x))*(1.25331414+y*(0.23498619
         +y*(-0.3655620e-1+y*(0.1504268e-1+y*(-0.780353e-2
         +y*(0.325614e-2+y*(-0.68245e-3)))))));
   }
   return ans;
}

/* Numerical recipes code
 
const double k1pi[]={0.5,5.598072040178741e-2,1.818666382168295e-3,
2.397509908859959e-5,1.239567816344855e-7};
const double k1qi[]={9.870202601341150e-1,1.292092053534579e-2,
5.881933053917096e-5};
const double k1p[]={-3.079657578292062e-1,-8.109417631822442e-2,
-3.477550948593604e-3,-5.385594871975406e-5,-3.110372465429008e-7};
const double k1q[]={9.861813171751389e-1,1.375094061153160e-2,
6.774221332947002e-5};
const double k1pp[]={1.253314137315502,1.457171340220454e1,
6.063161173098803e1,1.147386690867892e2,1.040442011439181e2,
4.356596656837691e1,7.265230396353690,3.144418558991021e-1};
const double k1qq[]={1.0,1.125154514806458e1,4.427488496597630e1,
7.616113213117645e1,5.863377227890893e1,1.850303673841586e1,
1.857244676566022,2.538540887654872e-2};

double k1(const double x) {
    // Returns the modiﬁed Bessel function K1(x) for positive real x.
    if (x <= 1.0) {  // Use two rational approximations.
        const double z=x*x;
        const double term = poly(k1pi,4,z)*log(x)/poly(k1qi,2,1.-z);
        return x*(poly(k1p,4,z)/poly(k1q,2,1.-z)+term)+1./x;
    } else {         // Rational approximation with e^{-x}/sqrt(x) factored out.
        const double z=1.0/x;
        return exp(-x)*poly(k1pp,7,z)/(poly(k1qq,7,z)*sqrt(x));
    }
}

inline double poly(const double *cof, const int n, const double x) {
    // Common code: Evaluate a polynomial.
    double ans = cof[n];
    for (int i=n-1;i>=0;i--) ans = ans*x+cof[i];
    return ans;
}
*/