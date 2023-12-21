// GL => OK

#ifndef CT_USER_H
#define CT_USER_H

#include <Cantera.h>
#include <onedim.h>
#include <IdealGasMix.h>
#include <equilibrium.h>
#include <transport.h>
#include <thermo.h>
#include <zerodim.h>
#include <iostream>
#include <iomanip>

using namespace Cantera;
using namespace std;
using namespace Cantera_CXX;

extern doublereal m_phi;



namespace User
{

  void test(int t);
  double flame_distance(double T1, double T2, double P1, double P2, double* Y1, double* Y2, int nsp);
  void findStoechiometryFromTmax(IdealGasMix* Fuel, IdealGasMix* Oxidizer, double* zst);
  void findStoechiometryFromCompo(IdealGasMix* Fuel, IdealGasMix* Oxidizer, double* zst);
  void getYfYo(IdealGasMix* mixture, double& Yf, double& Yo);
  void remesh_and_converge(Sim1D* flame, double criterion=1.0e-3, int loglevel=0);
  void goto_state(Sim1D* flame, double alpha=1.0e-6, int loglevel=0);
  void createFlameFromScratch(Sim1D* flame, int nint0=5, double lz=0.003, int loglevel=0);

}

#endif  // CT_USER_H
