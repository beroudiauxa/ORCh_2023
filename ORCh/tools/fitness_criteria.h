#ifndef FITNESS_CRITERIA_H
#define FITNESS_CRITERIA_H

#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>


//Read and Write mechanisms
#include "../read_write/read.h"
#include "../read_write/write.h"
#include "../read_write/write_QSS.h"
#include "../read_write/write_QSS_FORTRAN.h"

//Cantera properties
#include <Cantera.h>
#include <IdealGasMix.h>
#include <equilibrium.h>
#include <transport.h>
#include <zerodim.h>
#include <user.h>
using namespace Cantera;
using namespace Cantera_CXX;
using namespace User;

#include "../cantera/flamemodel.h"


using namespace std;

//---------------------
       


      void fit_function_0D(string mech_ref, string mech, double& fitness, int nbInlets, vector<string> listTargets, vector<string> trajectory_ref, int nbToKeep, string step, int rank);


      void fit_function_1D(string mech_ref, string mech, double& fitness, int nbFlames, vector<string> listTargets, vector<string> trajectory_ref, int nbToKeep, string step, int rank);



#endif

