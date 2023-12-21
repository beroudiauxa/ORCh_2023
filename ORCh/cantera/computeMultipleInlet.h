#ifndef computeMultipleInlet_H
#define computeMultipleInlet_H

//General properties
#include <mpi.h>
#include <iostream>
#include <complex>
#include <time.h>
#include <ctime>
#include <fstream>
#include "mpi.h"



//ORCh
#include "flamemodel.h"
#include "../drgep/drgep.h"
#include "particle.h"
#include "../read_write/read.h"


//Cantera properties
#include <Cantera.h>
#include <user.h>

using namespace Cantera;
using namespace Cantera_CXX;
using namespace User;

// Huu-Tri@20200724 : Add Tensorflow libraries "cppflow" - Use Tensorflow C-API to load ANN model and predict in C++
// /home/fkissel/workdir/orch/ORCh/cppflow
// CppFlow : https://github.com/serizba/cppflow
#include "../cppflow/include/Model.h"
#include "../cppflow/include/Tensor.h"


#define PI 3.14159265359

using namespace std;

//--------------------
class computeMultipleInlet
{
   public:
   //constructeur
   computeMultipleInlet();

   virtual void getMultipleInlet(string mech, string mech_desc, vector<MultipleInlet*> listInlets, vector<bool> Targets,
                          bool new_mixing, string step, vector<vector<vector<double> > > &R_AD_Trajectories, vector<vector<double> > &max_j_on_Target,
                          vector<vector<vector<double> > > &Ym_Trajectories_store, vector<vector<vector<double> > > &Production_Trajectories_ref,
                          vector<vector<vector<double> > > &Consumption_Trajectories_ref, vector<vector<double> > &T_Trajectories_store, vector<double> &time_store, vector<bool> &SpeciesIntoReactants);

   virtual void Next_Time_Step_with_drgep(vector<bool> Targets, double P, vector<double> &Ym, double &Hm, double &Tm, double delta_t,
                          vector<vector<double> > &R_AD_Trajectories, vector<vector<double> > &max_j_on_Target, string step);

   virtual void Next_Time_Step(double P, vector<double> &Ym, double &Hm, double &Tm, double delta_t);

   virtual void Next_Time_Step(double P, vector<double> &Ym, double &Hm, double &Tm, double delta_t,
                    vector<vector<vector<double> > > &Production_Trajectories_ref, vector<vector<vector<double> > > &Consumption_Trajectories_ref, int nInlet, int nLine);

   virtual void getMixedGasesComposition(vector<MultipleInlet*> listInlets, string step);

   virtual void Reacting(vector<Particle*> &listParticles, int nsp, double dt, double Pressure);
   virtual void ReactingParallel(vector<Particle*> &listParticles, int nsp, double dt, double Pressure);




   //destructeur
   virtual ~computeMultipleInlet();

   private:
	int Ifi_rank, Ila_rank, nb_var_loc, nproc;
	int *RecvCounts, *Disp;

	IdealGasMix *mixture;
};


#endif

