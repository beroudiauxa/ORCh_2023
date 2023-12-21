#ifndef computePremixedFlames_H
#define computePremixedFlames_H

//General properties
#include <sstream>
#include <time.h>


//ORCh
#include "../drgep/drgep.h"
#include "flamemodel.h"
#include "../tools/outputs.h"

//Cantera properties
#include <Cantera.h>
#include <user.h>
#include <onedim.h>

using namespace Cantera;
using namespace Cantera_CXX;
using namespace User;




   void output_datas (string outputName, vector<vector<double> > Y, vector<double> T , vector<double> U, vector<double> position, double position_max_wdot, IdealGasMix* mixture);

   void computePremixedFlames(string mech, string outputName, string mech_desc, vector<PremixedFlames*> listFlames, vector<bool> Targets, string step, vector<bool> &SpeciesIntoReactants, 
                              vector<vector<double> > &R_AD_Premixed, vector<vector<double> > &max_j_on_Target, vector<vector<double> > &max_jf_on_Target, vector<vector<double> > &max_jr_on_Target, 
                              vector<vector<double> > &QSS_Criteria, bool Proceed_DRGEP_analysis);

#endif
