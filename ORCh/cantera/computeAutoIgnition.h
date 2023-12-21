#ifndef computeAutoIgnition_H
#define computeAutoIgnition_H

//General properties
#include <complex>
#include <time.h>
#include <ctime>
#include <fstream>

//ORCh
#include "flamemodel.h"

//Cantera properties
#include <Cantera.h>
#include <user.h>

using namespace Cantera;
using namespace Cantera_CXX;
using namespace User;


   void computeAutoIgnition(string mech, string outputName, string mech_desc, vector<AutoIgnition*> listIgnitions);

#endif
