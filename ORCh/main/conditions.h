#ifndef conditions_H
#define conditions_H

#include "../cantera/flamemodel.h"
#include "../read_write/QSSscenario.h"
#include "../optimisation/OptimScenario.h"

   void conditions(
                  int &debuglevel, 
                  vector<string> &speciesToPlot,
		  string &mech_ref, 
                  string &mech, 
                  string &mech_desc, 
                  string &configuration, 
                  vector<MultipleInlet*> &listInlets, 
                  vector<PremixedFlames*> &listFlames, 
                  vector<string> &listTargets, 
                  string &step, 
                  bool &new_mixing,
		  bool &plot_T,
		  bool &plot_U, 
                  vector<QSSscenario*> &listQSSscenarios, 
                  OptimScenario* &listOptimScenarios,
		  vector<string> &trajectory_ref);
#endif
