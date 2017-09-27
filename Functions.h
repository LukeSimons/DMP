#include <string>

#include <cmath>
#include <math.h>
#include <stdio.h>

// Empirical fit to secondary electron emission equation as in Stangeby
double sec(double Te, char material);

void WarnOnce(bool &T, std::string Message);

double maxwellian(double E, double T);

// Ion backscattering
double ionback(double E, char isotope, char material, int flag);

double backscatter(double Te, double Ti, double mi, double Vion, char material, double &RE, double &RN);


double round_to_digits(double value, int digits);

/* Lambert W function. 
   Was ~/C/LambertW.c written K M Briggs Keith dot Briggs at bt dot com 97 May 21.  
   Revised KMB 97 Nov 20; 98 Feb 11, Nov 24, Dec 28; 99 Jan 13; 00 Feb 23; 01 Apr 09

   Computes Lambert W function, principal branch.
   See LambertW1.c for -1 branch.

   Returned value W(z) satisfies W(z)*exp(W(z))=z
   test data...
      W(1)= 0.5671432904097838730
      W(2)= 0.8526055020137254914
      W(20)=2.2050032780240599705
   To solve (a+b*R)*exp(-c*R)-d=0 for R, use
   R=-(b*W(-exp(-a*c/b)/b*d*c)+a*c)/b/c

   Test: 
     gcc -DTESTW LambertW.c -o LambertW -lm && LambertW
   Library:
     gcc -O3 -c LambertW.c 
*/

double LambertW(const double z);


