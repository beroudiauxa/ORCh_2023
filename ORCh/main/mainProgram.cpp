#include "mainProgram.h"
#include "mainDRGEPSpecies.h"
#include "mainDRGEPReactions.h"
#include "mainComputeTrajectories.h"
#include "mainLumping.h"














int comparison_main (const void *x, const void *y)
{
   double xx = *(double*)x;
   double yy = *(double*)y;
   if (xx < yy) return -1;
   if (xx > yy) return 1;
   return 0;
}



int main(int argc, char *argv[])
{
   MPI_Init(&argc, &argv);

   int debuglevel;
   vector<string> speciesToPlot;
   string mech_ref; string mech; string mech_desc;
   vector<MultipleInlet*> listInlets;
   vector<PremixedFlames*> listFlames;
   vector<AutoIgnition*> listIgnitions;
   vector<QSSscenario*> listQSSscenarios;
   vector<string> listTargets;
   OptimScenario* listOptimScenarios;

   string step = "";
   string configuration = "";
   bool new_mixing;
   bool plot_T;
   bool plot_U;
   vector<string> trajectory_ref; 

   conditions(debuglevel, speciesToPlot, mech_ref, mech, mech_desc, configuration, listInlets, listFlames, listTargets, step, new_mixing, plot_T, plot_U, listQSSscenarios, listOptimScenarios, trajectory_ref);

   int nInlets = listInlets.size();
   

   IdealGasMix *mixture  = new IdealGasMix(mech,mech_desc);
   int nsp_ref = mixture->nSpecies();
   int nreac_ref = mixture->nReactions();

   vector<Species_ORCh*> listSpecies_ref;
   vector<Reaction_ORCh*> listReactions;

   Read *r = new Read();
   r->Read_species(mech, listSpecies_ref);
   r->Read_reactions(mech, listReactions);

   vector<bool> Targets (nsp_ref, false);
   for (int k=0; k<nsp_ref; k++)
   {
      for (unsigned int kt=0; kt<listTargets.size(); kt++)
      {
         if (listSpecies_ref[k]->m_Name == listTargets[kt])
            Targets[k] = true;
      }
   }


   int nbLines = 0;
   int nbInlets = 0;

   if (listInlets.size() > 1)
   {
      nbInlets = listInlets.size()-1; //minus 1 because last inlet is always for burned gases
      nbLines =  dynamic_cast <Characteristics_MultipleInlet*> (listInlets[nbInlets])->m_NbIterations;
   }


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


   //Treat parallel stuff
   int rank;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);

   if (rank == 0)
   {
      cout << endl;
      cout << "                  _______                   _____                    _____                    _____          " << endl;
      cout << "                 /::\\    \\                 /\\    \\                  /\\    \\                  /\\    \\         " << endl;
      cout << "                /::::\\    \\               /::\\    \\                /::\\    \\                /::\\____\\        " << endl;
      cout << "               /::::::\\    \\             /::::\\    \\              /::::\\    \\              /:::/    /        " << endl;
      cout << "              /::::::::\\    \\           /::::::\\    \\            /::::::\\    \\            /:::/    /         " << endl;
      cout << "             /:::/~~\\:::\\    \\         /:::/\\:::\\    \\          /:::/\\:::\\    \\          /:::/    /          " << endl;
      cout << "            /:::/    \\:::\\    \\       /:::/__\\:::\\    \\        /:::/  \\:::\\    \\        /:::/____/           " << endl;
      cout << "           /:::/    / \\:::\\    \\     /::::\\   \\:::\\    \\      /:::/    \\:::\\    \\      /::::\\    \\           " << endl;
      cout << "          /:::/____/   \\:::\\____\\   /::::::\\   \\:::\\    \\    /:::/    / \\:::\\    \\    /::::::\\    \\    " << endl;
      cout << "         |:::|    |     |:::|    | /:::/\\:::\\   \\:::\\____\\  /:::/    /   \\:::\\    \\  /:::/\\:::\\    \\  " << endl;
      cout << "         |:::|____|     |:::|    |/:::/  \\:::\\   \\:::|    |/:::/____/     \\:::\\____\\/:::/  \\:::\\    \\" << endl;
      cout << "          \\:::\\    \\   /:::/    / \\::/   |::::\\  /:::|____|\\:::\\    \\      \\::/    /\\::/    \\:::\\    \\" << endl;
      cout << "           \\:::\\    \\ /:::/    /   \\/____|:::::\\/:::/    /  \\:::\\    \\      \\/____/  \\/____/ \\:::\\    \\ " << endl;
      cout << "            \\:::\\    /:::/    /          |:::::::::/    /    \\:::\\    \\                       \\:::\\____\\ " << endl;
      cout << "             \\:::\\__/:::/    /           |::|\\::::/    /      \\:::\\    \\                       \\:::|    |     " << endl;
      cout << "              \\::::::::/    /            |::| \\::/____/        \\:::\\    \\                      /:::|____|    " << endl;
      cout << "               \\::::::/    /             |::|  ~|               \\:::\\    \\                    /:::/    /     " << endl;
      cout << "                \\::::/    /              |::|   |                \\:::\\    \\                  /:::/    /      " << endl;
      cout << "                 \\::/____/               \\::|   |                 \\:::\\____\\                /:::/    /       " << endl;
      cout << "                  ~~                      \\:|   |                  \\::/    /                \\::/    /        " << endl;
      cout << "                                           \\|___|                   \\/____/                  \\/____/   by jaouen" << endl;
      cout << endl;
      cout << endl;
      cout << endl;
   }
                                                                                                    









//---------------------------------------------------------------------------//
//                                                                           //
//----------Proceed to the DRGEP reduction on the number of species----------//
//                                                                           //
//---------------------------------------------------------------------------//


if (step == "DRGEP_Species")
{
   if (rank == 0)
   {
      cout << endl;
      cout << "----------------------------------------------------------------------------------------------------------------" << endl << endl;
      cout << "--------------------------------------- STARTING DRGEP SPECIES STEP --------------------------------------------" << endl << endl;
      cout << "----------------------------------------------------------------------------------------------------------------" << endl << endl;

      print_to_screen_with_newline("MECHANISM:------------------------------------------------------------------------------------------------------", 0, debuglevel);
      print_to_screen("            Reading initial mechanism \"", 0, debuglevel);
      print_to_screen(mech, 0, debuglevel);
      print_to_screen("\" with description \"", 0, debuglevel);
      print_to_screen(mech_desc, 0, debuglevel);
      print_to_screen("\" ----------> ", 0, debuglevel);
      print_to_screen_with_newline("OK", 0, debuglevel);

      int nsp_init = mixture->nSpecies();
      int nreac_init = mixture->nReactions();

      print_to_screen("               Number of species: ", 0, debuglevel);
      print_to_screen_with_newline(nsp_init, 0, debuglevel);
      print_to_screen("               Number of reactions: ", 0, debuglevel);
      print_to_screen_with_newline(nreac_init, 0, debuglevel);
      print_to_screen_with_newline("", 0, debuglevel);

   }

   drgepSpecies(debuglevel,speciesToPlot, listTargets, configuration, mech, mech_desc, listInlets, listFlames, Targets, new_mixing, plot_T, plot_U, trajectory_ref, mech_ref, rank);
}
 
if (step == "DRGEP_Reactions")
{

   if (rank == 0)
   {
      cout << endl;
      cout << "----------------------------------------------------------------------------------------------------------------" << endl << endl;
      cout << "-------------------------------------- STARTING DRGEP REACTIONS STEP -------------------------------------------" << endl << endl;
      cout << "----------------------------------------------------------------------------------------------------------------" << endl << endl;

      print_to_screen_with_newline("MECHANISM:------------------------------------------------------------------------------------------------------", 0, debuglevel);
      print_to_screen("            Reading initial mechanism \"", 0, debuglevel);
      print_to_screen(mech, 0, debuglevel);
      print_to_screen("\" with description \"", 0, debuglevel);
      print_to_screen(mech_desc, 0, debuglevel);
      print_to_screen("\" ----------> ", 0, debuglevel);
      print_to_screen_with_newline("OK", 0, debuglevel);
   
      int nsp_init = mixture->nSpecies();
      int nreac_init = mixture->nReactions();
    
      print_to_screen("               Number of species: ", 0, debuglevel);
      print_to_screen_with_newline(nsp_init, 0, debuglevel);
      print_to_screen("               Number of reactions: ", 0, debuglevel);
      print_to_screen_with_newline(nreac_init, 0, debuglevel);
      print_to_screen_with_newline("", 0, debuglevel);
   }

   drgepReactions(debuglevel, speciesToPlot, listTargets, configuration, mech, mech_desc, listInlets, listFlames, Targets, new_mixing, plot_T, plot_U, trajectory_ref, mech_ref, rank);
}


if (step == "ComputeTrajectories")
 {
    if (rank == 0)
    {
       cout << endl;
       cout << "----------------------------------------------------------------------------------------------------------------" << endl << endl;
       cout << "----------------------------------- STARTING THE TRAJECTORIES COMPUTING ----------------------------------------" << endl << endl;
       cout << "----------------------------------------------------------------------------------------------------------------" << endl << endl;

         print_to_screen("            Reading mechanism \"", 0, debuglevel);
         print_to_screen(mech, 0, debuglevel);
         print_to_screen("\" with description \"", 0, debuglevel);
         print_to_screen(mech_desc, 0, debuglevel);
         print_to_screen("\" ----------> ", 0, debuglevel);
         print_to_screen_with_newline("OK", 0, debuglevel);
      
         int nsp_init = mixture->nSpecies();
         int nreac_init = mixture->nReactions();
       
         print_to_screen("               Number of species: ", 0, debuglevel);
         print_to_screen_with_newline(nsp_init, 0, debuglevel);
         print_to_screen("               Number of reactions: ", 0, debuglevel);
         print_to_screen_with_newline(nreac_init, 0, debuglevel);
         print_to_screen_with_newline("", 0, debuglevel);
    }

  ComputeTrajectories(debuglevel,speciesToPlot, listTargets, configuration, mech, mech_desc, listInlets, listFlames, Targets, new_mixing, plot_T, plot_U, trajectory_ref, mech_ref, rank);
 }

if (step == "Lumping")
   Lumping(debuglevel, configuration, mech, trajectory_ref[0], mech_desc);



//---------------------------------------------------------------------------//
//                                                                           //
//----------Proceed to the QSS criteria computation--------------------------//
//                                                                           //
//---------------------------------------------------------------------------//


if (step == "computeQSSCriteria")
{

   if (rank == 0)
   {
      cout << endl;
      cout << "----------------------------------------------------------------------------------------------------------------" << endl << endl;
      cout << "--------------------------------------- COMPUTING QSS CRITERIA STEP --------------------------------------------" << endl << endl;
      cout << "----------------------------------------------------------------------------------------------------------------" << endl << endl;

      print_to_screen_with_newline("MECHANISM:------------------------------------------------------------------------------------------------------", 0, debuglevel);
      print_to_screen("            Reading initial mechanism \"", 0, debuglevel);
      print_to_screen(mech, 0, debuglevel);
      print_to_screen("\" with description \"", 0, debuglevel);
      print_to_screen(mech_desc, 0, debuglevel);
      print_to_screen("\" ----------> ", 0, debuglevel);
      print_to_screen_with_newline("OK", 0, debuglevel);
   
      int nsp_init = mixture->nSpecies();
      int nreac_init = mixture->nReactions();
    
      print_to_screen("               Number of species: ", 0, debuglevel);
      print_to_screen_with_newline(nsp_init, 0, debuglevel);
      print_to_screen("               Number of reactions: ", 0, debuglevel);
      print_to_screen_with_newline(nreac_init, 0, debuglevel);
      print_to_screen_with_newline("", 0, debuglevel);

   }

   if (configuration == "MultipleInlet")
   {
      ofstream log("computeQSSCriteria.log");
      log << endl;


      computeMultipleInlet *c = new computeMultipleInlet();
      c->getMultipleInlet(mech, mech_desc, listInlets, Targets, new_mixing, step, R_AD_Trajectories, max_j_on_Target, Ym_Trajectories_ref, Production_Trajectories_ref, Consumption_Trajectories_ref, T_Trajectories_ref, time_store, SpeciesIntoReactants);


      string refPath = "./analytic_schemes/Ref_QSS";
      stringstream s_nsp_ref;
      s_nsp_ref << nsp_ref;
      refPath.append(s_nsp_ref.str()).append("_");

      OutputDeterministicTrajectories (nInlets, nbLines, listSpecies_ref , refPath , Ym_Trajectories_ref, T_Trajectories_ref, time_store);

   for (int k=0; k<nsp_ref; k++)
      QSS_Criteria[k][1] = k;

   for (int n=0; n<nInlets-1; n++)
   {
      vector<double> inte_Net (nsp_ref, 0.0);
      vector<double> inte_Production (nsp_ref, 0.0);
      vector<double> inte_Consumption (nsp_ref, 0.0);

      for (int i=0; i<nbLines; i++)
      {
         for (int k=0; k<nsp_ref; k++)
         {
            inte_Net[k] += abs(Production_Trajectories_ref[n][i][k]-Consumption_Trajectories_ref[n][i][k])*(time_store[2]-time_store[1]);
            inte_Production[k] += Production_Trajectories_ref[n][i][k]*(time_store[2]-time_store[1]);
            inte_Consumption[k] += Consumption_Trajectories_ref[n][i][k]*(time_store[2]-time_store[1]);
         }
      }

      for (int k=0; k<nsp_ref; k++)
      {

         if (QSS_Criteria[k][0] < inte_Net[k]/max(inte_Production[k], inte_Consumption[k]))
            QSS_Criteria[k][0] = inte_Net[k]/max(inte_Production[k], inte_Consumption[k]);
      }
   }

   for (int k=0; k<nsp_ref; k++)
   {
      if (listSpecies_ref[k]->m_Name == "N2")
         QSS_Criteria[k][0] = 1;

      if (Targets[k])
      {
         QSS_Criteria[k][0] = 1;
      }
   }


      log << "mech = " << mech ;
      log << endl  << endl << "--------QSS coefficients--------" << endl;
      log << "--------------------------------" << endl;
      
      if (rank == 0)
      { 
         cout << endl  << endl << endl << "--------QSS coefficients--------" << endl;
         cout << "--------------------------------" << endl;
      }


   for (int k=0; k<nsp_ref; k++)
   {
      if (rank == 0)
         cout << "Species " << listSpecies_ref[int(QSS_Criteria[k][1])]->m_Name << "  " << QSS_Criteria[k][0] << endl;
      log << "Species " << listSpecies_ref[int(QSS_Criteria[k][1])]->m_Name << "  " << QSS_Criteria[k][0] << endl;

   }

   if (rank == 0)
      cout << "--------------------------------" << endl << endl;

   log.close();

   Write_QSS *w_qss = new Write_QSS();
   w_qss->Check_Non_Linearity(mech, QSS_Criteria,rank);
   }


   if (configuration == "PremixedFlames")
   {


      string refPath = "./analytic_schemes/Ref_QSS";
      stringstream s_nsp_ref;
      s_nsp_ref << nsp_ref;
      refPath.append(s_nsp_ref.str());

      computePremixedFlames(mech, refPath, mech_desc, listFlames, Targets, step, SpeciesIntoReactants, R_AD_Premixed, max_j_on_Target, max_jf_on_Target, max_jr_on_Target, QSS_Criteria, false);

      Write_QSS *w_qss = new Write_QSS();
      w_qss->Check_Non_Linearity(mech, QSS_Criteria, rank);

   }

}






//---------------------------------------------------------------------------//
//                                                                           //
//----------Proceed to the QSS file creation---------------------------------//
//                                                                           //
//---------------------------------------------------------------------------//


if (step == "getQSSfile")
{

   if (rank == 0)
   {
      cout << endl;
      cout << "----------------------------------------------------------------------------------------------------------------" << endl << endl;
      cout << "---------------------------------------- STARTING GET QSS FILE STEP --------------------------------------------" << endl << endl;
      cout << "----------------------------------------------------------------------------------------------------------------" << endl << endl;

      print_to_screen_with_newline("MECHANISM:------------------------------------------------------------------------------------------------------", 0, debuglevel);
      print_to_screen("            Reading initial mechanism \"", 0, debuglevel);
      print_to_screen(mech, 0, debuglevel);
      print_to_screen("\" with description \"", 0, debuglevel);
      print_to_screen("\" ----------> ", 0, debuglevel);
      print_to_screen_with_newline("OK", 0, debuglevel);
   
      int nsp_init = mixture->nSpecies();
      int nreac_init = mixture->nReactions();
    
      print_to_screen("               Number of species: ", 0, debuglevel);
      print_to_screen_with_newline(nsp_init, 0, debuglevel);
      print_to_screen("               Number of reactions: ", 0, debuglevel);
      print_to_screen_with_newline(nreac_init, 0, debuglevel);
      print_to_screen_with_newline("", 0, debuglevel);

   }

   string Dim;
   if (configuration == "MultipleInlet" || configuration == "AutoIgnition")
   {
      Dim = "0D";
      
   }

   if (configuration == "PremixedFlames")
   {
      Dim = "1D";
   }


   for (unsigned int s=0; s<listQSSscenarios.size(); s++)
   {
      vector<bool> QSS_Species (nsp_ref, false);
     

      for (unsigned int k=0; k<listQSSscenarios[s]->m_QSSspecies.size(); k++)
      {
         for (int kName=0; kName<nsp_ref; kName++)
         {
            if (listSpecies_ref[kName]->m_Name == listQSSscenarios[s]->m_QSSspecies[k])
               QSS_Species[kName] = true;
         }
      }

      vector<bool> Species_to_add (nsp_ref, true);
      for (int k=0; k<nsp_ref; k++)
      {
         if (QSS_Species[k])
            Species_to_add[k] = false;
      }

      stringstream s_nbQSS;
      s_nbQSS << s;

      string path = "./analytic_schemes/RefQSSAnalysis";
      path.append(s_nbQSS.str());


      string copy_inputs = "cp *.ini ";
      copy_inputs.append(path).append("/");

      Optim *o = new Optim();
      if (configuration == "MultipleInlet")
         o->CreateAnalyticDirectory(configuration, path, listInlets.size(), "QSS", trajectory_ref);
      if (configuration == "PremixedFlames")
         o->CreateAnalyticDirectory(configuration, path, listFlames.size(), "QSS", trajectory_ref);
        
      string chem_file = path;
      chem_file.append("/mech_QSS");
      chem_file.append(".h");

      

      Write *w = new Write();
      Write_QSS *w_qss = new Write_QSS();

      int ret0 = system(copy_inputs.c_str());
      if (ret0 == -1) cout << "ERROR system";

      w_qss->Write_QSS_file(Dim, mech, chem_file.c_str(), QSS_Species, false, vector<double> (nreac_ref, 0.0), vector<double> (nreac_ref, 0.0), vector<double> (nreac_ref, 0.0));

      string scheme_file = path;	
      scheme_file.append("/scheme");
      scheme_file.append(".xml");

      w->Write_xml_for_Analytic_Applications(mech, mech_desc, scheme_file.c_str(), Species_to_add);

      string execute = "cd ";
      execute.append(path).append("; sh launch_cantera.sh");
      int ret = system(execute.c_str());
      if (ret == -1) cout << "ERROR with launch calculation";


   }




}




//---------------------------------------------------------------------------//
//                                                                           //
//----------Proceed to the optimisation of the chemical constants------------//
//                                                                           //
//---------------------------------------------------------------------------//


if (step == "Optimisation")
{

   if (rank == 0)
   {
      cout << endl;
      cout << "----------------------------------------------------------------------------------------------------------------" << endl << endl;
      cout << "---------------------------------------- STARTING OPTIMISATION STEP --------------------------------------------" << endl << endl;
      cout << "----------------------------------------------------------------------------------------------------------------" << endl << endl;
   }


   vector<bool> QSS_Species (nsp_ref, false);

   for (unsigned int k=0; k<listQSSscenarios[0]->m_QSSspecies.size(); k++)
   {
      for (int kName=0; kName<nsp_ref; kName++)
      {
         if (listSpecies_ref[kName]->m_Name == listQSSscenarios[0]->m_QSSspecies[k])
            QSS_Species[kName] = true;
      }
   }

   string copy_inputs0 = "cp *.ini ./analytic_schemes";
   int ret01 = system(copy_inputs0.c_str());
   if (ret01 == -1) cout << "ERROR system";


//copy the trajectory ref into the analytic_scheme directory
if (!trajectory_ref[0].empty() )
{
   for (int i=0; i<trajectory_ref.size();i++)
   {
      stringstream s_trajectory_ref;
      s_trajectory_ref << trajectory_ref[i];
      string copy_ref = "cp ./outputs/";
      if (configuration == "MultipleInlet")         
         copy_ref.append("Stochastic/").append(s_trajectory_ref.str()).append("_* ./analytic_schemes ");
      if (configuration == "PremixedFlames")
         copy_ref.append("Premixed/").append(s_trajectory_ref.str()).append(".dat ./analytic_schemes ");
      int ret = system(copy_ref.c_str());
      if (ret == -1) cout << "ERROR system";
   }
}

   Optim *o = new Optim();
   if (configuration == "MultipleInlet")
      o->Optimise(mech_ref, mech, mech_desc, QSS_Species, listOptimScenarios, listInlets.size(), listTargets, configuration, trajectory_ref);
   if (configuration == "PremixedFlames")
      o->Optimise(mech_ref, mech, mech_desc, QSS_Species, listOptimScenarios, listFlames.size(), listTargets, configuration, trajectory_ref);





}



//---------------------------------------------------------------------------//
//                                                                           //
//----------Proceed to the QSS file writing in FORTRAN-----------------------//
//                                                                           //
//---------------------------------------------------------------------------//


if (step == "getQSSfileFORTRAN")
{

   cout << "Converting the xml mech file : "<< mech <<" into FORTRAN, in the analytic_schemes/Ref directory" << endl;

   vector<bool> QSS_Species (nsp_ref, false);

   for (unsigned int k=0; k<listQSSscenarios[0]->m_QSSspecies.size(); k++)
   {
      for (int kName=0; kName<nsp_ref; kName++)
      {
         if (listSpecies_ref[kName]->m_Name == listQSSscenarios[0]->m_QSSspecies[k])
            QSS_Species[kName] = true;
      }
   }

   Write_QSS_FORTRAN *w_qss_fortran = new Write_QSS_FORTRAN();
   w_qss_fortran->Write_QSS_file_in_FORTRAN(mech, "./analytic_schemes/Ref/analytical.f90", QSS_Species);
}




   MPI_Finalize();

   return 0;

}



