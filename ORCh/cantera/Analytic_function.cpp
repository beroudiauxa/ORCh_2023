#include "Analytic_function.h"

int Mechanism = 2;

void Reduced(doublereal *w)
{
  std::cout << "Reduced function for 1D applications" << std::endl;
}


void Reduced(doublereal *w, doublereal Temp, doublereal *mass, doublereal density, const doublereal *mw)
{
   std::cout << "Reduced function for 0D applications" << std::endl;
}

