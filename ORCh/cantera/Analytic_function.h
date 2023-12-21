#ifndef Analytic_function_H
#define Analytic_function_H

#include <iostream>
#include <sstream>

//Cantera properties
#include <Cantera.h>
//#include <IdealGasMix.h>
//#include <equilibrium.h>
//#include <transport.h>
//#include <zerodim.h>
//#include <user.h>

using namespace Cantera;
//using namespace Cantera_CXX;
//using namespace User;
  
   void Reduced(doublereal *w);
   void Reduced(doublereal *w, doublereal Temp, doublereal *mass, doublereal density, const doublereal *mw);

#endif



