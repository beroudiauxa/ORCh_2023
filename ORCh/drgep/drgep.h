#ifndef DRGEP_H
#define DRGEP_H

//General properties
#include <iostream>
using namespace std;

//Cantera properties
#include <Cantera.h>
#include <user.h>
using namespace Cantera;
using namespace Cantera_CXX;
using namespace User;

//---------------------
class drgep 
{
   public:
   //constructeur
   drgep();

   virtual void drgep_0D_species(IdealGasMix *mixture, vector<bool> Targets, vector<vector<double> > &R_AD_Trajectories, int n, double time) const;

   virtual void drgep_0D_reactions(IdealGasMix *mixture, vector<vector<double> > &rj_for_k) const;

   virtual void drgep_1D_species(IdealGasMix *mixture, StFlow* flow, int ino, vector<bool> Targets, vector<vector<double> > &R_AD_Trajectories, double position_max_wdot) const;

   virtual void drgep_1D_reactions(IdealGasMix *mixture, StFlow* flow, int ino, vector<vector<double> > &rj_for_k, vector<vector<double> > &rjf_for_k, vector<vector<double> > &rjr_for_k) const;

   virtual void drgep_species(IdealGasMix *mixture, vector<bool> Targets, vector<vector<double> > &R_AD_Trajectories, vector<double> fwdRates, vector<double> revRates, string outputName) const;

   virtual void drgep_reactions(IdealGasMix *mixture, vector<vector<double> > &rj_for_k, vector<vector<double> > &rjf_for_k, vector<vector<double> > &rjr_for_k, vector<double> fwdRates, vector<double> revRates) const;

   //destructeur 
   virtual ~drgep();

   private:

};

#endif

