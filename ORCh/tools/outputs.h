#ifndef outputs_H
#define outputs_H

//General properties
#include <sstream>
#include <fstream>


#include "../read_write/read.h"

//Cantera properties
#include <Cantera.h>
#include <user.h>
using namespace Cantera;
using namespace Cantera_CXX;
using namespace User;


   void OutputDeterministicTrajectories (int nInlets, 
            int nbIterations, 
            vector<Species_ORCh*> listSpecies, 
            string outputName, 
            vector<vector<vector<double> > > Ym_Trajectories, 
            vector<vector<double> > T_Trajectories, 
            vector<double> time_store);


   void output_datas (string outputName,
                      vector<vector<double> > Y, vector<double> T ,
                      vector<double> U, 
                      vector<double> position, 
                      double position_max_wdot, 
                      IdealGasMix* mixture);

   void plotFitness(string path);

#endif
