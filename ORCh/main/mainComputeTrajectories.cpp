#include "mainComputeTrajectories.h"


void ComputeTrajectories(int debuglevel, vector<string> speciesToPlot, vector<string> listTargets,  string configuration, string mech, string mech_desc, vector<MultipleInlet*> listInlets, vector<PremixedFlames*> listFlames, vector<bool> Targets, bool new_mixing, bool plot_T, bool plot_U, vector<string> trajectory_ref,string mech_ref, int rank)
{




   IdealGasMix *mixture  = new IdealGasMix(mech,mech_desc);

   int nbLines = 0;
   int nbInlets = 0;
   if (listInlets.size() > 1)
   {
      nbInlets = listInlets.size()-1; //minus 1 because last inlet is always for burned gases
      nbLines =  dynamic_cast <Characteristics_MultipleInlet*> (listInlets[nbInlets])->m_NbIterations;
   }


   int nsp_ref = mixture->nSpecies();
   int nreac_ref = mixture->nReactions();

   vector<double> time_store (nbLines, 0.0);
   vector<vector<vector<double> > > Ym_Trajectories_ref (nbInlets, vector<vector<double> > (nbLines, vector<double> (nsp_ref, 0.0)));
   vector<vector<double> > T_Trajectories_ref (nbInlets, vector<double> (nbLines, 0.0));
   vector<vector<vector<double> > > Production_Trajectories_ref (nbInlets, vector<vector<double> > (nbLines, vector<double> (nsp_ref, 0.0)));
   vector<vector<vector<double> > > Consumption_Trajectories_ref (nbInlets, vector<vector<double> > (nbLines, vector<double> (nsp_ref, 0.0)));

   vector<vector<vector<double> > > R_AD_Trajectories (nbInlets, vector<vector<double> > (nsp_ref, vector<double>(nsp_ref, 0.0)));
   vector<vector<double> > max_j_on_Target (nsp_ref, vector<double> (nreac_ref,0.0));
   vector<vector<double> > max_jf_on_Target (nsp_ref, vector<double> (nreac_ref,0.0));
   vector<vector<double> > max_jr_on_Target (nsp_ref, vector<double> (nreac_ref,0.0));


   vector<bool> SpeciesIntoReactants (nsp_ref, false);


   vector<vector<double> > R_AD_Premixed (nsp_ref, vector<double>(nsp_ref, 0.0));
   vector<vector<double> > QSS_Criteria (nsp_ref, vector<double> (2, 0.0));


   vector<Species_ORCh*> listSpecies_ref;

   Read *r = new Read();
   r->Read_species(mech, listSpecies_ref);


   if (configuration == "MultipleInlet")
   {
      computeMultipleInlet *c = new computeMultipleInlet();
      c->getMultipleInlet(mech, mech_desc, listInlets, Targets, new_mixing, "ComputeTrajectories", R_AD_Trajectories, max_j_on_Target, Ym_Trajectories_ref, Production_Trajectories_ref, Consumption_Trajectories_ref, T_Trajectories_ref, time_store, SpeciesIntoReactants);

      OutputDeterministicTrajectories (listInlets.size(), nbLines, listSpecies_ref , "./outputs/Stochastic/Trajectory_", Ym_Trajectories_ref, T_Trajectories_ref, time_store);
   }

   if (configuration == "PremixedFlames")
   {
      computePremixedFlames(mech, "./outputs/Premixed/Premixed_", mech_desc, listFlames, Targets, "ComputeTrajectories", SpeciesIntoReactants, R_AD_Premixed, max_j_on_Target, max_jf_on_Target, max_jr_on_Target, QSS_Criteria, false);
   }


   Script_gnuplot("ComputeTrajectories", mech, speciesToPlot, configuration, mech_desc, nsp_ref, plot_U, plot_T, trajectory_ref, mech_ref, nbInlets, listTargets, mech, rank ) ;
  
   string path="outputs/";
   if (configuration == "MultipleInlet")
      path.append("Stochastic/");
   if (configuration == "PremixedFlames")
      path.append("Premixed/");
   string graph = "cd ";
   graph.append(path).append("; ").append("gnuplot makeGnu.gnu");
   int ret = system(graph.c_str());
   if (ret == -1) cout << "ERROR system";



/*
   if (configuration == "AutoIgnition")
   {
      computeAutoIgnition(mech, "./outputs/AutoIgnition_", mech_desc, listIgnitions);
   }
*/

}
















