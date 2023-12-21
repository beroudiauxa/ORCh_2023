#ifndef GNUPLOT_H
#define GNUPLOT_H

#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>

//Tools


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
#include "fitness_criteria.h"

using namespace std;

//---------------------
   void Script_gnuplot ( string step, 
                         string initial_mech, 
                         vector<string> speciesToPlot, 
                         string configuration, 
                         string mech_desc, 
                         int nbSpeciesToKeep, 
                         bool plot_U, 
                         bool plot_T, 
                         vector<string> trajectory_ref, 
                         string mech_ref, 
                         int nbInlets,  
                         vector<string> listTarget,
                         string  outputSchemeName,
                         int rank);

#endif

