#include "mainDRGEPSpecies.h"

//----------------------------------------
//   <mainDRGEPSpecies> 
//----------------------------------------
//   DESCRIPTION:
//         Run the DRGEP analysis by
//            1- Computing reference trajectories
//            2- Computing the species DRGEP coefficients (R_AD)
//            3- 
//   IN: 
//         configuration: Either "Premixed" or "MultipleInlet" case
//   OUT:
//
//----------------------------------------
void drgepSpecies(int debuglevel, vector<string> speciesToPlot, vector<string> listTargets,  string configuration, string initial_mech, string mech_desc, vector<MultipleInlet*> listInlets, vector<PremixedFlames*> listFlames, vector<bool> Targets, bool new_mixing, bool plot_T, bool plot_U, vector<string> trajectory_ref,string mech_ref)
{

   cout << endl;
                        cout << "----------------------------------------------------------------------------------------------------------------------------------" << endl << endl;
                        cout << "------------------------------------------------ STARTING DRGEP SPECIES STEP -----------------------------------------------------" << endl << endl;
                        cout << "----------------------------------------------------------------------------------------------------------------------------------" << endl << endl;



   print_to_screen_with_newline("MECHANISM:------------------------------------------------------------------------------------------------------------------------", 0, debuglevel);
   print_to_screen("            Reading initial mechanism \"", 0, debuglevel);
   print_to_screen(initial_mech, 0, debuglevel);
   print_to_screen("\" with description \"", 0, debuglevel);
   print_to_screen(mech_desc, 0, debuglevel);
   print_to_screen("\" ----------> ", 0, debuglevel);
 


   IdealGasMix *mixture  = new IdealGasMix(initial_mech,mech_desc);
   print_to_screen_with_newline("OK", 0, debuglevel);

   int nsp_init = mixture->nSpecies();
   int nreac_init = mixture->nReactions();
   int nbFlame = listFlames.size()-1;   

   stringstream s_nbFlame;
   s_nbFlame << nbFlame;

   print_to_screen("               Number of species: ", 0, debuglevel);
   print_to_screen_with_newline(nsp_init, 0, debuglevel);
   print_to_screen("               Number of reactions: ", 0, debuglevel);
   print_to_screen_with_newline(nreac_init, 0, debuglevel);
   print_to_screen_with_newline("", 0, debuglevel);

   vector<Species*> listSpecies_init;
   vector<Reaction*> listReactions_init;

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

   double sort_R_AD [nsp_init][2]; //sort_R_AD[k][0] -> interaction coefficient //sort_R_AD[k][1] -> species number
   int e = 1; //epsilon
   double T_final_ref(0);

   ofstream log("DRGEP_Species.log");

   if (configuration == "MultipleInlet")
   {
      for (int k=0; k<nsp_init; k++)
      {
         sort_R_AD[k][0] = 0.0;
         sort_R_AD[k][1] = k;
      }

      string outputName = "./outputs/Stochastic/Ref_DRGEP_Species";
      stringstream s_nbSpeciesInit;
      s_nbSpeciesInit << nsp_init;
      outputName.append(s_nbSpeciesInit.str()).append("_");

      computeMultipleInlet *c = new computeMultipleInlet();
      c->getMultipleInlet(initial_mech, 
                          mech_desc, 
                          listInlets, 
                          Targets, 
                          new_mixing, 
                          "DRGEP_Species", 
                          R_AD_Trajectories, 
                          max_j_on_Target, 
                          Ym_Trajectories_ref, 
                          Production_Trajectories_ref, 
                          Consumption_Trajectories_ref, 
                          T_Trajectories_ref, 
                          time_store, 
                          SpeciesIntoReactants);

      OutputDeterministicTrajectories (listInlets.size(), nbIterations, listSpecies_init, outputName, Ym_Trajectories_ref, T_Trajectories_ref, time_store);

         T_final_ref = T_Trajectories_ref[0][nbIterations-1];

      for (unsigned int n=0; n<nbInlets; n++)
      {
         for (int ka=0; ka<nsp_init; ka++)
         {
            if (Targets[ka])
            {
               for (int kb=0; kb<nsp_init; kb++)
               {
                  if (sort_R_AD[kb][0] < R_AD_Trajectories[n][ka][kb])
                     sort_R_AD[kb][0] = R_AD_Trajectories[n][ka][kb];
               }
            }
         }
      }



   } //end if (configuration == "MultipleInlet")

   if (configuration == "PremixedFlames")
   {
      for (int k=0; k<nsp_init; k++)
      {
         sort_R_AD[k][0] = 0.0;
         sort_R_AD[k][1] = k;
      }

      string outputName = "./outputs/Premixed/Ref_DRGEP_Species";
      stringstream s_nbSpeciesInit;
      s_nbSpeciesInit << nsp_init;
      outputName.append(s_nbSpeciesInit.str());

      computePremixedFlames(initial_mech, 
                            outputName, 
                            mech_desc, 
                            listFlames, 
                            Targets, 
                            "DRGEP_Species", 
                            SpeciesIntoReactants, 
                            R_AD_Premixed, 
                            max_j_on_Target, 
                            max_jf_on_Target, 
                            max_jr_on_Target, 
                            QSS_Criteria, 
                            true);


      for (int ka=0; ka<nsp_init; ka++) //BE CAREFUL, FOR NOW DRGEP ANALYSIS ON PREMIXED FLAMES IS CODED FOR ONLY ONE FLAME, JUST NEED TO LOOP ON THE NFLAMES
      {
         if (Targets[ka])
         {
            for (int kb=0; kb<nsp_init; kb++)
            {
               if (sort_R_AD[kb][0] < R_AD_Premixed[ka][kb])
                  sort_R_AD[kb][0] = R_AD_Premixed[ka][kb];
            }
         }
      }
   } //end if (configuration == "PremixedFlames")

   //If the species is within the targets or within the initial reactants its R_AD coefficient is set to one to ensure it is kept while reducing the mechanism  
   for (int k=0; k<nsp_init; k++)
   {
      if (SpeciesIntoReactants[k])
         sort_R_AD[k][0] = 1;

      if (Targets[k])
         sort_R_AD[k][0] = 1;
   }



  
   qsort (sort_R_AD, sizeof(sort_R_AD)/sizeof(sort_R_AD[0]), sizeof(sort_R_AD[0]), compare_numbers);



   log << "mech = " << initial_mech << endl;
   log << "mech_ref = " << mech_ref << endl;
   log << "trajectory_ref = " << trajectory_ref[0] << endl;
   log << "DRGEP coefficients ---------->" << endl; 

   for (int k=0; k<nsp_init; k++)
   {
      cout << sort_R_AD[k][0] << "  " << listSpecies_init[sort_R_AD[k][1]]->m_Name << endl;
      log << sort_R_AD[k][0] << "  " << listSpecies_init[sort_R_AD[k][1]]->m_Name << endl;
   }



   log.close();
   cout << endl;

   for (int nbSpeciesToKeep=nsp_init-1; nbSpeciesToKeep>9; nbSpeciesToKeep--)
   {
      vector<bool> Species_to_add (nsp_init, false);
      vector<bool> Reactions_to_add (nreac_init, true);

      for (int k=0; k<nbSpeciesToKeep; k++)
         Species_to_add[int(sort_R_AD[nsp_init-1-k][1])] = true;

      string outputSchemeName = "./outputs/mechanisms/drgepSpecies";
      stringstream s_nbSpeciesToKeep;
      s_nbSpeciesToKeep << nbSpeciesToKeep;
      outputSchemeName.append(s_nbSpeciesToKeep.str()).append(".xml");

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
           vector<Species*> (),
           vector<Reaction*> ());

      //Ecriture des schema V2 correspondants
      string outputSchemeNameV2 = "./outputs/mechanisms/V2drgepSpecies";
      outputSchemeNameV2.append(s_nbSpeciesToKeep.str()).append(".xml");

      vector<bool> Species_to_addV2 (nsp_init, false);

      for (int k=0; k<nbSpeciesToKeep+1; k++)
         Species_to_addV2[int(sort_R_AD[nsp_init-k-1][1])] = true;

      Species_to_addV2[int(sort_R_AD[nsp_init-nbSpeciesToKeep][1])] = false;

      Write *w2 = new Write();
      w2->Write_xml_file(initial_mech,
           mech_desc,
           outputSchemeNameV2,
           Species_to_addV2,
           Reactions_to_add,
           false,
           vector<double> (nreac_init, 0.0),
           vector<double> (nreac_init, 0.0),
           vector<double> (nreac_init, 0.0),
           vector<Species*> (),
           vector<Reaction*> ());

   }


   double T_final(0);
   double Te = 500000; //50K limit between the final ref temperature and the final reduced temperature
   double fitness_final = 0.0;
   int last_test = 0; //determine how many final tests will occur 

   for (int nbSpeciesToKeep=nsp_init-1; nbSpeciesToKeep>9; nbSpeciesToKeep--)
   {

      string outputName = "./outputs/Stochastic/Reduced_DRGEP_Species";
      string outputSchemeName = "./outputs/mechanisms/drgepSpecies";
      stringstream s_nbSpeciesToKeep;
      s_nbSpeciesToKeep << nbSpeciesToKeep;
      outputSchemeName.append(s_nbSpeciesToKeep.str()).append(".xml");

      string saveScheme;

      if (fitness_final > Te)
      {
         nbSpeciesToKeep++;        
         s_nbSpeciesToKeep.seekp(ios::beg); 
         s_nbSpeciesToKeep << nbSpeciesToKeep;
         string outputSchemeNameV2 = "./outputs/mechanisms/V2drgepSpecies";
         outputSchemeNameV2.append(s_nbSpeciesToKeep.str()).append(".xml");
        

         saveScheme = outputSchemeName;
         outputSchemeName = outputSchemeNameV2;
         //count
         last_test++;
         cout << "Passage dans la boucle de reprise avec test du schema " << outputSchemeName << endl;
      }

      vector<Species*> listSpecies_drgep;
      Read *r = new Read();
      r->Read_species(outputSchemeName, listSpecies_drgep);
      int nsp_drgep = listSpecies_drgep.size();

      if (configuration == "MultipleInlet")
      {


         outputName.append(s_nbSpeciesToKeep.str()).append("_");

         vector<vector<vector<double> > > Ym_Trajectories (nbInlets, vector<vector<double> > (nbIterations, vector<double> (nsp_drgep, 0.0)));
         vector<vector<double> > T_Trajectories (nbInlets, vector<double> (nbIterations, 0.0));

         
         computeMultipleInlet *c = new computeMultipleInlet();
         c->getMultipleInlet(outputSchemeName, 
                             mech_desc, 
                             listInlets, 
                             Targets, 
                             false, 
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

         string outputName = "./outputs/Premixed/Reduced_DRGEP_Species";
         outputName.append(s_nbSpeciesToKeep.str()).append("_").append(s_nbFlame.str());

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
      
  



        Script_gnuplot("DRGEP_Species", initial_mech, speciesToPlot, configuration, mech_desc, nbSpeciesToKeep, plot_U, plot_T, trajectory_ref, mech_ref, nbInlets, listTargets,outputSchemeName ) ;




      string path="outputs/";
      if (configuration == "MultipleInlet")
         path.append("Stochastic/");
      if (configuration == "PremixedFlames")
         path.append("Premixed/");
      string graph = "cd ";
      graph.append(path).append("; ").append("gnuplot makeGnu.gnu");
      int ret = system(graph.c_str());
      if (ret == -1) cout << "ERROR system";

      if (configuration == "MultipleInlet")
      {
	
         if (last_test == 1 && fitness_final < Te)
         {
            cout << "DRGEP Species step terminated. For next step, use the scheme \"V2drgepSpecies" << nbSpeciesToKeep+1  <<".xml\" located in output/mechanisms." << endl; 
            break;
         }
         else if (last_test == 1 && fitness_final >= Te)
         {
            cout << "DRGEP Species step terminated. For next step, use the scheme \"drgepSpecies" << nbSpeciesToKeep+1  <<".xml\" located in output/mechanisms." << endl; 
            break;
         }
      }
  }// end for nbToKeep









}





