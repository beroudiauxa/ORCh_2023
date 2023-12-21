#include "mainDRGEPReactions.h"
#include "../tools/tools.h"

//----------------------------------------
//   <mainDRGEPReactions> 
//----------------------------------------
//   DESCRIPTION:
//         Run the DRGEP analysis by
//            1- Computing reference trajectories
//            2- Computing the reactions DRGEP coefficients (R_AD)
//            3- 
//   IN: 
//         configuration: Either "Premixed" or "MultipleInlet" case
//   OUT:
//
//----------------------------------------
void drgepReactions(int debuglevel, vector<string> speciesToPlot, vector<string> listTargets, string configuration, string initial_mech, string mech_desc, vector<MultipleInlet*> listInlets, vector<PremixedFlames*> listFlames, vector<bool> Targets, bool new_mixing, bool plot_T, bool plot_U, vector<string> trajectory_ref, string mech_ref, int rank)
{



   IdealGasMix *mixture  = new IdealGasMix(initial_mech,mech_desc);

   int nsp_init = mixture->nSpecies();
   int nreac_init = mixture->nReactions();
   int nbFlame = listFlames.size()-1;


   stringstream s_nbFlame;
   s_nbFlame << nbFlame;


   vector<Species_ORCh*> listSpecies_init;
   vector<Reaction_ORCh*> listReactions_init;

   Read *r = new Read();
   r->Read_species(initial_mech, listSpecies_init);
   r->Read_reactions(initial_mech, listReactions_init);


   int nbIterations = 0;
   int nbInlets = 0;

   if (listInlets.size() > 1)
   {
      nbInlets = listInlets.size()-1; //minus 1 because last inlet is always for burned gases
      nbIterations =  dynamic_cast <Characteristics_MultipleInlet*> (listInlets[nbInlets])->m_NbIterations;
   }


   vector<vector<vector<double> > > R_AD_Trajectories (nbInlets, vector<vector<double> > (nsp_init, vector<double>(nsp_init, 0.0)));
   vector<vector<double> > R_AD_Premixed (nsp_init, vector<double>(nsp_init, 0.0));
   vector<vector<double> > max_j_on_Target (nsp_init, vector<double> (nreac_init,0.0));
   vector<vector<double> > max_jf_on_Target (nsp_init, vector<double> (nreac_init,0.0));
   vector<vector<double> > max_jr_on_Target (nsp_init, vector<double> (nreac_init,0.0));
   vector<vector<vector<double> > > Production_Trajectories_ref (nbInlets, vector<vector<double> > (nbIterations, vector<double> (nsp_init, 0.0)));
   vector<vector<vector<double> > > Consumption_Trajectories_ref (nbInlets, vector<vector<double> > (nbIterations, vector<double> (nsp_init, 0.0)));
   vector<vector<double> > T_Trajectories_ref (nbInlets, vector<double> (nbIterations, 0.0));
   vector<vector<vector<double> > > Ym_Trajectories_ref (nbInlets, vector<vector<double> > (nbIterations, vector<double> (nsp_init, 0.0)));
   vector<double> time_store (nbIterations, 0.0);
   vector<bool> SpeciesIntoReactants (nsp_init, false);
   vector<vector<double> > QSS_Criteria (nsp_init, vector<double> (2, 0.0));

   if (configuration == "MultipleInlet")
   {

      string outputName = "./outputs/Stochastic/Ref_DRGEP_Reactions";
      stringstream s_nbReactionsInit;
      s_nbReactionsInit << nreac_init;
      outputName.append(s_nbReactionsInit.str()).append("_");

      computeMultipleInlet *c = new computeMultipleInlet();
      c->getMultipleInlet(initial_mech,
                          mech_desc,
                          listInlets,
                          Targets,
                          new_mixing,
                          "DRGEP_Reactions",
                          R_AD_Trajectories,
                          max_j_on_Target,
                          Ym_Trajectories_ref,
                          Production_Trajectories_ref,
                          Consumption_Trajectories_ref,
                          T_Trajectories_ref,
                          time_store,
                          SpeciesIntoReactants);

      OutputDeterministicTrajectories (listInlets.size(), nbIterations, listSpecies_init, outputName, Ym_Trajectories_ref, T_Trajectories_ref, time_store);

   } //end if (configuration == "MultipleInlet")

   if (configuration == "PremixedFlames")
   {
      string outputName = "./outputs/Premixed/Ref_DRGEP_Reactions";
      stringstream s_nbReactionsInit;
      s_nbReactionsInit << nreac_init;
      outputName.append(s_nbReactionsInit.str());

      computePremixedFlames(initial_mech,
                            outputName,
                            mech_desc,
                            listFlames,
                            Targets,
                            "DRGEP_Reactions",
                            SpeciesIntoReactants,
                            R_AD_Premixed,
                            max_j_on_Target,
                            max_jf_on_Target,
                            max_jr_on_Target,
                            QSS_Criteria,
                            true);

   } //end if (configuration == "PremixedFlames")

   double SortReactions [nreac_init][2];
   double SortForwardReactions [nreac_init][2];
   double SortReverseReactions [nreac_init][2];

   for (int j=0; j<nreac_init; j++)
   {
      SortReactions[j][1] = j;
      SortForwardReactions[j][1] = j;
      SortReverseReactions[j][1] = j;

      SortReactions[j][0] = 0.0;
      SortForwardReactions[j][0] = 0.0;
      SortReverseReactions[j][0] = 0.0;

      for (int ka=0; ka<nsp_init; ka++)
      {
         if (SortReactions[j][0] < max_j_on_Target[ka][j])
             SortReactions[j][0] = max_j_on_Target[ka][j];

         if (SortForwardReactions[j][0] < max_jf_on_Target[ka][j])
             SortForwardReactions[j][0] = max_jf_on_Target[ka][j];

         if (SortReverseReactions[j][0] < max_jr_on_Target[ka][j])
             SortReverseReactions[j][0] = max_jr_on_Target[ka][j];
      }
	
      /* Huu-Tri Nguyen - Add important reaction - To modify for different scheme - 2020.03.06*/
/*      if(j== 5-1)	// Number in .xml file - 1
      {
	SortReactions[j][0] = 1.0;
 	SortForwardReactions[j][0] = 1.0;
      	SortReverseReactions[j][0] = 1.0;
      }
*/


     // End add Huu-Tri Nguyen - 2020.03.06

   }

      

   qsort (SortReactions, sizeof(SortReactions)/sizeof(SortReactions[0]), sizeof(SortReactions[0]), compare_numbers);
   qsort (SortForwardReactions, sizeof(SortForwardReactions)/sizeof(SortForwardReactions[0]), sizeof(SortForwardReactions[0]), compare_numbers);
   qsort (SortReverseReactions, sizeof(SortReverseReactions)/sizeof(SortReverseReactions[0]), sizeof(SortReverseReactions[0]), compare_numbers);

   ofstream log("DRGEP_Reactions.log");

   log << "mech = " << initial_mech << endl;
   log << "mech_ref = " << mech_ref << endl;
   log << "trajectory_ref = " << trajectory_ref[0] << endl;
   log << "DRGEP coefficients ---------->" << endl;



  /* 
      for (int j=0; j<nreac_init; j++)
      {
   cout << endl  << "Forward reaction " << SortForwardReactions[j][1]+1 << "  " << SortForwardReactions[j][0] << endl;
   log << "Forward reaction " << SortForwardReactions[j][1]+1 << "  " << SortForwardReactions[j][0] << endl;
      }

      for (int j=0; j<nreac_init; j++)
      {
         cout << "Reverse reaction " << SortReverseReactions[j][1]+1 << "  " << SortReverseReactions[j][0] << endl;
         log << "Reverse reaction " << SortReverseReactions[j][1]+1 << "  " << SortReverseReactions[j][0] << endl;
      }
*/

   if (rank ==0)
   {
      cout << endl  << endl << endl << "-------DRGEP coefficients-------" << endl;
      cout << "--------------------------------" << endl;
   }


   for (int j=0; j<nreac_init; j++)
   {
      if (rank ==0)
         cout <<  "Reaction " << SortReactions[j][1]+1 << "  " << SortReactions[j][0] << endl;
      log << "Reaction " << SortReactions[j][1]+1 << "  " << SortReactions[j][0] << endl;
   }

   log.close();
   if (rank ==0)
      cout << "--------------------------------" << endl << endl;



   for (int nbReactionsToKeep=nreac_init-1; nbReactionsToKeep>10; nbReactionsToKeep--)
   {
      vector<bool> Species_to_add (nsp_init, true);
      vector<bool> Reactions_to_add (nreac_init, false);

      for (int j=0; j<nbReactionsToKeep; j++)
         Reactions_to_add[int(SortReactions[nreac_init-1-j][1])] = true;

      string outputSchemeName = "./outputs/mechanisms/drgepReactions";
      stringstream s_nbReactionsToKeep;
      s_nbReactionsToKeep << nbReactionsToKeep;
      outputSchemeName.append(s_nbReactionsToKeep.str()).append(".xml");

      Write *w = new Write();
      w->Write_xml_file(initial_mech, 
            mech_desc, 
            outputSchemeName, 
            Species_to_add, 
            Reactions_to_add, 
            false, 
            vector<double> (nreac_init, 0.0), 
            vector<double> (nreac_init, 0.0), 
            vector<double> (nreac_init, 0.0), 
            vector<Species_ORCh*> (), 
            vector<Reaction_ORCh*> ());
/*
      // Ecriture des schemas V2 correspondants
      string outputSchemeNameV2 = "./outputs/mechanisms/V2drgepReactions";
      outputSchemeNameV2.append(s_nbReactionsToKeep.str()).append(".xml");

      vector<bool> Reactions_to_addV2 (nreac_init, false);

      for (int k=0; k<nbReactionsToKeep+1; k++)
         Reactions_to_addV2[int(SortReactions[nreac_init-k-1][1])] = true;

      Reactions_to_addV2[int(SortReactions[nreac_init-nbReactionsToKeep][1])] = false;

      Write *w2 = new Write();
      w2->Write_xml_file(initial_mech, 
            mech_desc, 
            outputSchemeNameV2, 
            Species_to_add, 
            Reactions_to_addV2, 
            false, 
            vector<double> (nreac_init, 0.0), 
            vector<double> (nreac_init, 0.0), 
            vector<double> (nreac_init, 0.0), 
            vector<Species_ORCh*> (), 
            vector<Reaction_ORCh*> ());
*/


   }

   double T_final_ref(0);
   double T_final(0);
   double Te = 3000; //50K limit between the final ref temperature and the final reduced temperature
   double fitness_final = 0.0;
   int last_test = 0; //determine how many final tests will occur 



   for (int nbReactionsToKeep=nreac_init-1; nbReactionsToKeep>10; nbReactionsToKeep--)
   {

      string outputName = "./outputs/Stochastic/Reduced_DRGEP_Reactions";
      string outputSchemeName = "./outputs/mechanisms/drgepReactions";
      stringstream s_nbReactionsToKeep;
      s_nbReactionsToKeep << nbReactionsToKeep;
      outputSchemeName.append(s_nbReactionsToKeep.str()).append(".xml");

      string saveScheme;
/*
      if (fitness_final > Te)
      {
         nbReactionsToKeep++;

         string outputSchemeNameV2 = "./outputs/mechanisms/V2drgepReactions";
         s_nbReactionsToKeep.seekp(ios::beg);
         s_nbReactionsToKeep << nbReactionsToKeep;
         outputSchemeNameV2.append(s_nbReactionsToKeep.str()).append(".xml");


         saveScheme = outputSchemeName;
         outputSchemeName = outputSchemeNameV2;
         //count
         last_test++;
         cout << "Passage dans la boucle de reprise avec test du schema " << outputSchemeName << endl;
      }
*/

      vector<Species_ORCh*> listSpecies_drgep;
      Read *r = new Read();
      r->Read_species(outputSchemeName, listSpecies_drgep);
      int nsp_drgep = listSpecies_drgep.size();

      if (configuration == "MultipleInlet")
      {

         outputName.append(s_nbReactionsToKeep.str()).append("_");

         vector<vector<vector<double> > > Ym_Trajectories (nbInlets, vector<vector<double> > (nbIterations, vector<double> (nsp_drgep, 0.0)));
         vector<vector<double> > T_Trajectories (nbInlets, vector<double> (nbIterations, 0.0));

         computeMultipleInlet *c = new computeMultipleInlet();
         c->getMultipleInlet(outputSchemeName, 
                     mech_desc, 
                     listInlets, 
                     Targets, 
                     false /*no new mixing*/, 
                     "ComputeTrajectories", 
                     R_AD_Trajectories, 
                     max_j_on_Target, 
                     Ym_Trajectories, 
                     Production_Trajectories_ref, 
                     Consumption_Trajectories_ref, 
                     T_Trajectories, 
                     time_store, 
                     SpeciesIntoReactants);

         T_final = T_Trajectories[0][nbIterations-1];
         fitness_final = abs(T_final - T_final_ref);

         OutputDeterministicTrajectories (listInlets.size(), nbIterations, listSpecies_drgep, outputName, Ym_Trajectories, T_Trajectories, time_store);

      } //end if configuration == "MultipleInlet"


      if (configuration == "PremixedFlames")
      {
         string outputName = "./outputs/Premixed/Reduced_DRGEP_Reactions";
         outputName.append(s_nbReactionsToKeep.str()).append("_").append(s_nbFlame.str());

         computePremixedFlames(outputSchemeName,
                     outputName,
                     mech_desc,
                     listFlames,
                     Targets,
                     "ComputeTrajectories",
                     SpeciesIntoReactants,
                     R_AD_Premixed,
                     max_j_on_Target,
                     max_jf_on_Target,
                     max_jr_on_Target,
                     QSS_Criteria,
                     false);
      }

   Script_gnuplot( "DRGEP_Reactions" , initial_mech, speciesToPlot, configuration, mech_desc, nbReactionsToKeep, plot_U, plot_T, trajectory_ref, mech_ref, nbInlets, listTargets, outputSchemeName, rank ) ;



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
   if (last_test == 1 && fitness_final < Te)
   {
      cout << "DRGEP Reactions step terminated. For next step, use the scheme \"V2drgepReactions" << nbReactionsToKeep+1  <<".xml\" located in output/mechanisms." << endl;
      break;
   }
   else if (last_test == 1 && fitness_final >= Te)
   {
      cout << "DRGEP Reactions step terminated. For next step, use the scheme \"drgepSpecies" << nbReactionsToKeep+1  <<".xml\" located in output/mechanisms." << endl;
      break;
   }
*/

   } // end for nbToKeep






}







