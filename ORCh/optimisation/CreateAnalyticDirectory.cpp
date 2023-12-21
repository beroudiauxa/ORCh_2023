#include "optimisation.h"
#include <stdlib.h>
#include <fstream>


void Optim::CreateAnalyticDirectory(string configuration, string fileName, int nbInlets_nbFlames, string step, vector<string> trajectory_ref)
{

//-----------Remove existing Ref0 Ref1 ... directories----------//

   string eraseFolder = "rm -rf ";
   eraseFolder.append(fileName);

   int ret = system(eraseFolder.c_str());
   if (ret == -1) cout << "ERROR system";


//-----------Create new Ref0 Ref1 ... directories----------//

   string createFolder = "mkdir ";
   createFolder.append(fileName);

   ret = system(createFolder.c_str());
   if (ret == -1) cout << "ERROR system";


//-----------Create launch_cantera.sh file----------//

   string create_launch_cantera = fileName;
   create_launch_cantera.append("/launch_cantera.sh");
   ofstream launch_cantera (create_launch_cantera.c_str());
   launch_cantera << "make clean" << endl;
   launch_cantera << "make -j" << endl;
   launch_cantera << " echo \"Analytic trajectories calculating \""   << endl;

   if (configuration == "MultipleInlet")
   {
      if (step == "Optimisation")  
         launch_cantera << "./computeAnalyticFlame" << endl;

      if (step == "QSS")
         launch_cantera << "mpirun -np 8 computeAnalyticFlame" << endl;
   }

   if (configuration == "PremixedFlames")
//      launch_cantera << "timeout 60s ./computeAnalyticFlame" << endl;
      launch_cantera << "./computeAnalyticFlame" << endl;
   if (configuration == "AutoIgnition")
      launch_cantera << "timeout 60s ./computeAnalyticFlame" << endl;

   launch_cantera.close();



   //-----------Copy the reference trajectories and the input file----------//

      
      string copy_ref;
      if (configuration == "PremixedFlames")
      {  
         for (int i=0; i<nbInlets_nbFlames; i++) 
         {   
            stringstream s_traj_ref;
            s_traj_ref << trajectory_ref[i];
            copy_ref = "cp outputs/Premixed/";
            copy_ref.append(s_traj_ref.str()).append(".dat ");
         
            copy_ref.append(fileName);
            ret = system(copy_ref.c_str());
            if (ret == -1) cout << "ERROR system";
         }
      }

     if (configuration == "MultipleInlet")
     {
         for (int i=0; i<nbInlets_nbFlames-1; i++) 
         {   
            stringstream s_traj_ref;
            copy_ref = "cp outputs/Stochastic/";
            stringstream s_inlet;
            s_inlet << i;
            s_traj_ref << trajectory_ref[0];
            copy_ref.append(s_traj_ref.str()).append("_").append(s_inlet.str()).append(".dat ");
            copy_ref.append(fileName);
            int ret_ = system(copy_ref.c_str());
            if (ret_ == -1) cout << "ERROR system";

         }
     }         


      if (step == "Optimisation")
      {
         string copy_input = "cp analytic_schemes/*.ini ";
         copy_input.append(fileName);
         ret = system(copy_input.c_str());
         if (ret == -1) cout << "ERROR system";
      }
//----------Create the file for 0D computations----------//
//-------------------------------------------------------//
   if (configuration == "MultipleInlet")
   {
   //-----------Copy Selection_read.dat file----------//
   
      string copy_Selection_read = "cp Selection.dat ";
      copy_Selection_read.append(fileName);
      ret = system(copy_Selection_read.c_str());
      if (ret == -1) cout << "ERROR system";

   

   //-----------Create computeAnalyticFlame.cpp file----------//
        
      string create_computeAnalyticFlame = fileName;
      create_computeAnalyticFlame.append("/computeAnalyticFlame.cpp");
      ofstream computeAnalyticFlame (create_computeAnalyticFlame.c_str());
   
      computeAnalyticFlame << "#include <mpi.h>" << endl;
      computeAnalyticFlame << "#include \"../../../../ORCh/tools/outputs.h\"" << endl;
      computeAnalyticFlame << "#include \"../../../../ORCh/tools/gnuplot.h\"" << endl;
      computeAnalyticFlame << "#include \"../../../../ORCh/tools/fitness_criteria.h\"" << endl;
      computeAnalyticFlame << "#include \"../../../../ORCh/read_write/read.h\"" << endl; 
      computeAnalyticFlame << "#include \"../../../../ORCh/read_write/reaction.h\"" << endl; 
      computeAnalyticFlame << "#include \"../../../../ORCh/read_write/species.h\"" << endl; 
      computeAnalyticFlame << "" << endl;
      computeAnalyticFlame << "#include \"mech_QSS.h\"" << endl;
      computeAnalyticFlame << "" << endl;
      computeAnalyticFlame << "#include \"../../../../ORCh/main/conditions.h\"" << endl;
      computeAnalyticFlame << "#include \"../../../../ORCh/cantera/computeMultipleInlet.h\"" << endl;
      computeAnalyticFlame << "" << endl;
      computeAnalyticFlame << "#include <Cantera.h>" << endl;
      computeAnalyticFlame << "#include <IdealGasMix.h>    // defines class IdealGasMix" << endl;
      computeAnalyticFlame << "#include <equilibrium.h>    // chemical equilibrium" << endl;
      computeAnalyticFlame << "#include <transport.h>      // transport properties" << endl;
      computeAnalyticFlame << "#include <zerodim.h>" << endl;
      computeAnalyticFlame << "#include <user.h>" << endl;
      computeAnalyticFlame << "using namespace Cantera;" << endl;
      computeAnalyticFlame << "using namespace Cantera_CXX;" << endl;
      computeAnalyticFlame << "using namespace User;" << endl;
      computeAnalyticFlame << "" << endl;
      computeAnalyticFlame << "IdealGasMix* mixture;" << endl;
      computeAnalyticFlame << "string mech;" << endl;
      computeAnalyticFlame << "" << endl;
      computeAnalyticFlame << "" << endl;
  if (step == "QSS")
      computeAnalyticFlame << "int main(int argc, char *argv[])" << endl;
  else
      computeAnalyticFlame << "int main()" << endl;
      computeAnalyticFlame << "{" << endl;
      computeAnalyticFlame << "" << endl;
  
  if (step == "QSS")
     computeAnalyticFlame << "   MPI_Init(&argc, &argv);" << endl;
     computeAnalyticFlame << "   int rank=0;" << endl;
  if (step == "QSS")
     computeAnalyticFlame << "   MPI_Comm_rank(MPI_COMM_WORLD, &rank);" << endl;
      computeAnalyticFlame << "" << endl;
     // computeAnalyticFlame << "   cout << \"\" << endl;" << endl;
      computeAnalyticFlame << "" << endl;
      computeAnalyticFlame << "   string mech_ref; string mech; string mech_desc;" << endl;
      computeAnalyticFlame << "   vector<MultipleInlet*> listInlets;" << endl;
      computeAnalyticFlame << "   vector<PremixedFlames*> listFlames;" << endl;
      computeAnalyticFlame << "   vector<AutoIgnition*> listIgnitions;" << endl;
      computeAnalyticFlame << "   vector<string> listTargets;" << endl;
      computeAnalyticFlame << "   bool new_mixing;" << endl;
      computeAnalyticFlame << "   string step = \"\";" << endl;
      computeAnalyticFlame << "   string configuration = \"\";" << endl;
      computeAnalyticFlame << "   vector<QSSscenario*> listQSSscenarios;" << endl;
      computeAnalyticFlame << "   OptimScenario* listOptimScenarios;" << endl;
      computeAnalyticFlame << "   int debuglevel;" << endl;
      computeAnalyticFlame << "   vector<string> speciesToPlot;" << endl;
      computeAnalyticFlame << "   bool plot_T;" << endl;
      computeAnalyticFlame << "   bool plot_U;" << endl;
      computeAnalyticFlame << "   vector<string> trajectory_ref;" << endl;
      computeAnalyticFlame << "" << endl;
      computeAnalyticFlame << "   conditions(debuglevel, speciesToPlot, mech_ref, mech, mech_desc, configuration, listInlets, listFlames, listIgnitions, listTargets, step, new_mixing, plot_T, plot_U, listQSSscenarios, listOptimScenarios, trajectory_ref);" << endl;
      computeAnalyticFlame << "" << endl;
      computeAnalyticFlame << "   string temp = \"../../\";" << endl;
      computeAnalyticFlame << "   stringstream s_mech_ref;" << endl;
      computeAnalyticFlame << "   s_mech_ref << mech_ref;" << endl;
      computeAnalyticFlame << "   temp.append(s_mech_ref.str());" << endl;
      computeAnalyticFlame << "   mech_ref = temp;" << endl;
      computeAnalyticFlame << "" << endl;
      computeAnalyticFlame << "" << endl;
      computeAnalyticFlame << "" << endl;
      computeAnalyticFlame << "  string initial_mech = \"../../\"; " << endl;
      computeAnalyticFlame << "  stringstream s_mech;" << endl;
      computeAnalyticFlame << "  s_mech << mech;" << endl;
      computeAnalyticFlame << "  initial_mech.append(s_mech.str());" << endl;
      computeAnalyticFlame << "" << endl;
      computeAnalyticFlame << "   mech = \"scheme.xml\";" << endl;
      computeAnalyticFlame << "" << endl;
      computeAnalyticFlame << "   mixture  = new IdealGasMix(mech, mech_desc);" << endl;
      computeAnalyticFlame << "   vector<Species_ORCh*> listSpecies;" << endl;
      computeAnalyticFlame << "   Read *r = new Read();" << endl;
      computeAnalyticFlame << "   r->Read_species(mech, listSpecies);" << endl;
      computeAnalyticFlame << "" << endl;
      computeAnalyticFlame << "   vector<bool> Targets (mixture->nSpecies(), false);" << endl;
      computeAnalyticFlame << "   vector<bool> SpeciesIntoReactants (mixture->nSpecies(), false);" << endl;
      computeAnalyticFlame << "   int nbLines =  dynamic_cast <Characteristics_MultipleInlet*> (listInlets[listInlets.size()-1])->m_NbIterations;" << endl;
      computeAnalyticFlame << "   vector<vector<vector<double> > > R_AD_Trajectories (listInlets.size()-1, vector<vector<double> > (mixture->nSpecies(), vector<double>(mixture->nSpecies(), 0.0)));" << endl;
      computeAnalyticFlame << "   vector<vector<double> > max_j_on_Target (mixture->nSpecies(), vector<double> (mixture->nReactions(),0.0));" << endl;
      computeAnalyticFlame << "   vector<vector<vector<double> > > Ym_Trajectories_ref (listInlets.size()-1, vector<vector<double> > (nbLines, vector<double> (mixture->nSpecies(), 0.0)));" << endl;
      computeAnalyticFlame << "   vector<vector<vector<double> > > Production_Trajectories_ref (listInlets.size()-1, vector<vector<double> > (nbLines, vector<double> (mixture->nSpecies(), 0.0)));" << endl;
      computeAnalyticFlame << "   vector<vector<vector<double> > > Consumption_Trajectories_ref (listInlets.size()-1, vector<vector<double> > (nbLines, vector<double> (mixture->nSpecies(), 0.0)));" << endl;
      computeAnalyticFlame << "   vector<vector<double> > T_Trajectories_ref (listInlets.size()-1, vector<double> (nbLines, 0.0));" << endl;
      computeAnalyticFlame << "   vector<double> time_store (nbLines, 0.0);" << endl;
      computeAnalyticFlame << "" << endl;
      computeAnalyticFlame << "   computeMultipleInlet *c = new computeMultipleInlet();" << endl;
      computeAnalyticFlame << "   c->getMultipleInlet(mech, mech_desc, listInlets, Targets, new_mixing, step, R_AD_Trajectories, max_j_on_Target, Ym_Trajectories_ref," << endl;
      computeAnalyticFlame << "                        Production_Trajectories_ref, Consumption_Trajectories_ref, T_Trajectories_ref, time_store, SpeciesIntoReactants);" << endl;
      computeAnalyticFlame << "" << endl;
      computeAnalyticFlame << "" << endl;
      computeAnalyticFlame << "   string redPath = \"Reduced_QSS\";" << endl;
      computeAnalyticFlame << "   stringstream s_nsp;" << endl;
      computeAnalyticFlame << "   s_nsp << listSpecies.size();" << endl;
      computeAnalyticFlame << "   redPath.append(s_nsp.str()).append(\"_\");" << endl;
      computeAnalyticFlame << "" << endl;
      computeAnalyticFlame << "   OutputDeterministicTrajectories (listInlets.size(), nbLines, listSpecies, redPath, Ym_Trajectories_ref, T_Trajectories_ref, time_store);" << endl;
      computeAnalyticFlame << "   Script_gnuplot (\"QSS\", initial_mech, speciesToPlot, configuration, mech_desc, listSpecies.size() , plot_U, plot_T, trajectory_ref, mech_ref, listInlets.size()-1, listTargets, mech, rank); " << endl;
      computeAnalyticFlame << "" << endl;
      computeAnalyticFlame << "   string graph = \"gnuplot makeGnu.gnu\";" << endl;
      computeAnalyticFlame << "   int ret = system(graph.c_str());" << endl;
      computeAnalyticFlame << "   if (ret == -1) cout << \"ERROR system\";" << endl;
      computeAnalyticFlame << "" << endl;

   if (step == "QSS")
      computeAnalyticFlame << "   MPI_Finalize();" << endl;

      computeAnalyticFlame << "}" << endl;
      computeAnalyticFlame.close();
   }


//----------Create the file for 1D computations----------//
//-------------------------------------------------------//
   if (configuration == "PremixedFlames")
   {
      string create_computeAnalyticFlame = fileName;
      create_computeAnalyticFlame.append("/computeAnalyticFlame.cpp");
      ofstream computeAnalyticFlame (create_computeAnalyticFlame.c_str());

      computeAnalyticFlame << "#include \"../../../../ORCh/tools/outputs.h\"" << endl;
      computeAnalyticFlame << "#include \"../../../../ORCh/tools/gnuplot.h\"" << endl;
      computeAnalyticFlame << "#include \"../../../../ORCh/tools/fitness_criteria.h\"" << endl;
      computeAnalyticFlame << "#include \"../../../../ORCh/read_write/read.h\"" << endl; 
      computeAnalyticFlame << "#include \"../../../../ORCh/read_write/reaction.h\"" << endl; 
      computeAnalyticFlame << "#include \"../../../../ORCh/read_write/species.h\"" << endl; 
      computeAnalyticFlame << "#include \"mech_QSS.h\"" << endl;
      computeAnalyticFlame << "#include \"../../../../ORCh/main/conditions.h\"" << endl;
      computeAnalyticFlame << endl;
      computeAnalyticFlame << "#include <Cantera.h> "   << endl;
      computeAnalyticFlame << "#include <user.h>"   << endl;
      computeAnalyticFlame << "using namespace Cantera;"   << endl;
      computeAnalyticFlame << "using namespace Cantera_CXX;"   << endl;
      computeAnalyticFlame << "using namespace User;"   << endl;
      computeAnalyticFlame << endl;
      computeAnalyticFlame << "IdealGasMix* mixture;" << endl;
      computeAnalyticFlame << "StFlow* flow;" << endl;
      computeAnalyticFlame << endl;
      computeAnalyticFlame << "   "   << endl;
      computeAnalyticFlame << "int main()"   << endl;
      computeAnalyticFlame << "{"   << endl;
      computeAnalyticFlame << "   string mech_ref; string mech; string mech_desc; string configuration;" << endl;
      computeAnalyticFlame << "   vector<MultipleInlet*> listInlets;" << endl;
      computeAnalyticFlame << "   vector<PremixedFlames*> listFlames;" << endl;
      computeAnalyticFlame << "   vector<AutoIgnition*> listIgnitions;" << endl;
      computeAnalyticFlame << "   vector<string> listTargets;" << endl;
      computeAnalyticFlame << "   bool new_mixing;" << endl;
      computeAnalyticFlame << "   string step = \"\";" << endl;
      computeAnalyticFlame << "   vector<QSSscenario*> listQSSscenarios;" << endl;
      computeAnalyticFlame << "   OptimScenario* listOptimScenarios;" << endl;
      computeAnalyticFlame << "   int debuglevel;" << endl;
      computeAnalyticFlame << "   vector<string> speciesToPlot;" << endl;
      computeAnalyticFlame << "   bool plot_T;" << endl;
      computeAnalyticFlame << "   bool plot_U;" << endl;
      computeAnalyticFlame << "   vector<string> trajectory_ref;" << endl;
      computeAnalyticFlame << endl;
      computeAnalyticFlame << "   conditions(debuglevel, speciesToPlot, mech_ref, mech, mech_desc, configuration, listInlets, listFlames, listIgnitions, listTargets, step, new_mixing, plot_T, plot_U, listQSSscenarios, listOptimScenarios, trajectory_ref);" << endl;
      computeAnalyticFlame << endl; 
      computeAnalyticFlame << "  string initial_mech = \"../../\"; " << endl;
      computeAnalyticFlame << "  stringstream s_mech;" << endl;
      computeAnalyticFlame << "  s_mech << mech;" << endl;
      computeAnalyticFlame << "  initial_mech.append(s_mech.str());" << endl;
      computeAnalyticFlame << "" << endl;
      computeAnalyticFlame << "" << endl;
      computeAnalyticFlame << "   string temp = \"../../\";" << endl;
      computeAnalyticFlame << "   stringstream s_mech_ref;" << endl;
      computeAnalyticFlame << "   s_mech_ref << mech_ref;" << endl;
      computeAnalyticFlame << "   temp.append(s_mech_ref.str());" << endl;
      computeAnalyticFlame << "   mech_ref = temp;" << endl;
      computeAnalyticFlame << "" << endl;
      computeAnalyticFlame << "   mech = \"scheme.xml\";" << endl;
      computeAnalyticFlame << endl;
      computeAnalyticFlame << "   mixture  = new IdealGasMix(mech,mech_desc);" << endl;
      computeAnalyticFlame << "   flow = new StFlow(mixture);" << endl;
      computeAnalyticFlame << "   flow->setFreeFlow();" << endl;
      computeAnalyticFlame << endl;
      computeAnalyticFlame << "   Inlet1D* inlet = new Inlet1D;" << endl;
      computeAnalyticFlame << "   Outlet1D* outlet = new Outlet1D;" << endl;
      computeAnalyticFlame << " " << endl;
      computeAnalyticFlame << "   vector<Domain1D*> domains;" << endl;
      computeAnalyticFlame << "   domains.push_back(inlet);" << endl;
      computeAnalyticFlame << "   domains.push_back(flow);" << endl;
      computeAnalyticFlame << "   domains.push_back(outlet);" << endl;
      computeAnalyticFlame << "   int flowdomain=1; // 0 is inlet, 1 is flow, 2 is outlet" << endl;
      computeAnalyticFlame << "   Sim1D* flame = new Sim1D(domains);" << endl;
      computeAnalyticFlame << endl;
      computeAnalyticFlame << "   //----------Number of species----------//" << endl;
      computeAnalyticFlame << "   int nsp = mixture->nSpecies();" << endl;
      computeAnalyticFlame << "   int nreac = mixture->nReactions();" << endl;
      computeAnalyticFlame << "" << endl;
      computeAnalyticFlame << "   for (int n=0;n<listFlames.size();n++)" << endl;
      computeAnalyticFlame << "   {" << endl;
      computeAnalyticFlame << "      //----------Set properties----------//" << endl;
      computeAnalyticFlame << "      flow->solveEnergyEqn();" << endl;
      computeAnalyticFlame << "      Transport* trmix = newTransportMgr(\"Mix\", mixture);" << endl;
      computeAnalyticFlame << "      flow->setTransport(*trmix);" << endl;
      computeAnalyticFlame << "      flow->setKinetics(*mixture);" << endl;
      computeAnalyticFlame << "      string initial_flame = \"../../\";" << endl;
      computeAnalyticFlame << "      initial_flame.append(listFlames[n]->m_Initial_Flame);" << endl;
      computeAnalyticFlame << "      cout << endl <<  \"initial_flame \" << initial_flame << endl;" << endl;
      computeAnalyticFlame << "      flame->restore(initial_flame, \"st_flame\");" << endl;
      computeAnalyticFlame << endl;
      computeAnalyticFlame << endl;
      computeAnalyticFlame << "      double *Ym = new double[nsp];" << endl;
      computeAnalyticFlame << "      double Hm = 0.0;" << endl;
      computeAnalyticFlame << endl;
      computeAnalyticFlame << "      //----------Composition mixture----------//" << endl;
      computeAnalyticFlame << "      // Initialize Fuel and Oxidizer" << endl;
      computeAnalyticFlame << "      IdealGasMix *fuel     = new IdealGasMix(mech,mech_desc);" << endl;
      computeAnalyticFlame << "      IdealGasMix *oxidizer = new IdealGasMix(mech,mech_desc);" << endl;
      computeAnalyticFlame << "      IdealGasMix *stoichio = new IdealGasMix(mech,mech_desc);" << endl;
      computeAnalyticFlame << endl;
      computeAnalyticFlame << "      // fuel state" << endl;
      computeAnalyticFlame << "       if (listFlames[0]->m_Yf_str != \"\")" << endl;
      computeAnalyticFlame << "      fuel->setState_TPY(listFlames[n]->m_T_fuel, listFlames[n]->m_Pressure, listFlames[n]->m_Yf_str);" << endl;
      computeAnalyticFlame << "       if (listFlames[0]->m_Xf_str != \"\")" << endl;
      computeAnalyticFlame << "      fuel->setState_TPX(listFlames[n]->m_T_fuel, listFlames[n]->m_Pressure, listFlames[n]->m_Xf_str);" << endl;
      computeAnalyticFlame << "      double Hf  = fuel->enthalpy_mass();" << endl;
      computeAnalyticFlame << "      double *Yf = new double[nsp];" << endl;
      computeAnalyticFlame << "      double *Xf = new double[nsp];" << endl;
      computeAnalyticFlame << "      fuel->getMassFractions(Yf);" << endl;
      computeAnalyticFlame << "      fuel->getMoleFractions(Xf);" << endl;
      computeAnalyticFlame << endl;
      computeAnalyticFlame << "      // oxidizer state" << endl;
      computeAnalyticFlame << "      if (listFlames[0]->m_Yo_str != \"\")" << endl;
      computeAnalyticFlame << "      oxidizer->setState_TPY(listFlames[n]->m_T_oxidizer, listFlames[n]->m_Pressure, listFlames[n]->m_Yo_str);" << endl;
      computeAnalyticFlame << "      if (listFlames[0]->m_Xo_str != \"\")" << endl;
      computeAnalyticFlame << "      oxidizer->setState_TPX(listFlames[n]->m_T_oxidizer, listFlames[n]->m_Pressure, listFlames[n]->m_Xo_str);" << endl;
      computeAnalyticFlame << "      double Ho  = oxidizer->enthalpy_mass();" << endl;
      computeAnalyticFlame << "      double *Yo = new double[nsp];" << endl;
      computeAnalyticFlame << "      double *Xo = new double[nsp];" << endl;
      computeAnalyticFlame << "      oxidizer->getMassFractions(Yo);" << endl;
      computeAnalyticFlame << "      oxidizer->getMoleFractions(Xo);" << endl;
      computeAnalyticFlame << endl;
      computeAnalyticFlame << "      // get stoichiometry" << endl;
      computeAnalyticFlame << "      double zst;" << endl;
      computeAnalyticFlame << "      findStoechiometryFromCompo(fuel,oxidizer,&zst);" << endl;
      computeAnalyticFlame << endl;
      computeAnalyticFlame << "      // get fuel and oxidizer mass fraction in fuel" << endl;
      computeAnalyticFlame << "      double Yf_F,Yo_F;" << endl;
      computeAnalyticFlame << "      getYfYo(fuel,Yf_F,Yo_F);" << endl;
      computeAnalyticFlame << endl;
      computeAnalyticFlame << "      // get fuel and oxidizer mass fraction in oxidizer" << endl;
      computeAnalyticFlame << "      double Yf_O,Yo_O;" << endl;
      computeAnalyticFlame << "      getYfYo(oxidizer,Yf_O,Yo_O);" << endl;
      computeAnalyticFlame << endl;
      computeAnalyticFlame << "      // Stoichiometry state" << endl;
      computeAnalyticFlame << "      double *Yst = new double[nsp];" << endl;
      computeAnalyticFlame << "      for (int ns=0;ns<nsp;ns++)" << endl;
      computeAnalyticFlame << "         Yst[ns] = zst*Yf[ns]+(1.0-zst)*Yo[ns];" << endl;
      computeAnalyticFlame << "      stoichio->setMassFractions(Yst);" << endl;
      computeAnalyticFlame << "      double Hst = zst*Hf+ (1.0-zst)*Ho;" << endl;
      computeAnalyticFlame << "      stoichio->setState_HP(Hst,listFlames[n]->m_Pressure);" << endl;
      computeAnalyticFlame << endl;
      computeAnalyticFlame << "     // get fuel and oxidizer mass fraction in stoichiometric mixture" << endl;
      computeAnalyticFlame << "     double Yf_st,Yo_st;" << endl;
      computeAnalyticFlame << "     getYfYo(stoichio,Yf_st,Yo_st);" << endl;
      computeAnalyticFlame << endl;
      computeAnalyticFlame << endl;
      computeAnalyticFlame << "     // compute useful quantities" << endl;
      computeAnalyticFlame << "     double s = Yo_st/Yf_st;" << endl;
      computeAnalyticFlame << "     double Z1_O = s*Yf_O - Yo_O;" << endl;
      computeAnalyticFlame << "     double Z1_F = s*Yf_F - Yo_F;" << endl;
      computeAnalyticFlame << endl;
      computeAnalyticFlame << "     double z = (Yo_O*(listFlames[n]->m_Phi-1.0) - Z1_O)/((Z1_F-Z1_O)-(listFlames[n]->m_Phi-1.0)*(Yo_F-Yo_O));" << endl;
      computeAnalyticFlame << endl;
      computeAnalyticFlame << "     // Mixture state" << endl;
      computeAnalyticFlame << "     for (int ns=0;ns<nsp;ns++)" << endl;
      computeAnalyticFlame << "        Ym[ns] = z*Yf[ns]+(1.0-z)*Yo[ns];" << endl;
      computeAnalyticFlame << "     Hm = z*Hf+ (1.0-z)*Ho;" << endl;
      computeAnalyticFlame << endl;
      computeAnalyticFlame << "     mixture->setMassFractions(Ym);" << endl;
      computeAnalyticFlame << "     mixture->setState_HP(Hm, listFlames[n]->m_Pressure);" << endl;
      computeAnalyticFlame << endl;
      computeAnalyticFlame << "      try" << endl;
      computeAnalyticFlame << "      {" << endl;
      computeAnalyticFlame << endl;
      computeAnalyticFlame << "           cout<<\"               Adjusting flame ...\"<<endl;" << endl;
      computeAnalyticFlame << "         {" << endl;
      computeAnalyticFlame << "           double ratio=2;" << endl;
      computeAnalyticFlame << "           double slope=0.05;" << endl;
      computeAnalyticFlame << "           double curve=0.8;" << endl;
      computeAnalyticFlame << "           int loglevel = 1;" << endl;
      computeAnalyticFlame << "           double coeff_start=1.0e-3;" << endl;
      computeAnalyticFlame << "           flame->setRefineCriteria(flowdomain,ratio,slope,curve,-1);" << endl;
      computeAnalyticFlame << "           flame->set_max_isteps(30);" << endl;
      computeAnalyticFlame << "           flame->setJacAge(30);" << endl;
      computeAnalyticFlame << "           flame->newton().set_jacMax(30);" << endl;
      computeAnalyticFlame << "           flame->newton().setOptions(30);" << endl;
      computeAnalyticFlame << "           double mesh_crit=1.0e-3;" << endl;
      computeAnalyticFlame << "           goto_state(flame,coeff_start,loglevel);" << endl;
      computeAnalyticFlame << "           remesh_and_converge(flame,mesh_crit,loglevel);" << endl;
      computeAnalyticFlame << "         }" << endl;
      computeAnalyticFlame << endl;
      computeAnalyticFlame << "         cout<<\"               Converging flame ...\"<<endl;" << endl;
      computeAnalyticFlame << "         {" << endl;
      computeAnalyticFlame << "           double ratio=2;" << endl;
      computeAnalyticFlame << "           double slope=0.05;" << endl;
      computeAnalyticFlame << "           double curve=0.8;" << endl;
      computeAnalyticFlame << "           int loglevel = 1;" << endl;
      computeAnalyticFlame << "           double mesh_crit=1.0e-3;" << endl;
      computeAnalyticFlame << "           flame->setRefineCriteria(flowdomain,ratio,slope,curve,-1);" << endl;
      computeAnalyticFlame << "           flame->set_max_isteps(30);" << endl;
      computeAnalyticFlame << "           flame->setJacAge(30);" << endl;
      computeAnalyticFlame << "           flame->newton().set_jacMax(20);" << endl;
      computeAnalyticFlame << "           flame->newton().setOptions(20);" << endl;
      computeAnalyticFlame << "           remesh_and_converge(flame,mesh_crit,loglevel);" << endl;
      computeAnalyticFlame << "         }" << endl;
      computeAnalyticFlame << "      }" << endl;
      computeAnalyticFlame << "      catch(...)" << endl;
      computeAnalyticFlame << "      {" << endl;
      computeAnalyticFlame << "      cout << \"An exception occurred. Flame not converged\" << endl; "<< endl;
      computeAnalyticFlame << "         throw -1;" << endl;
      computeAnalyticFlame << "      }" << endl;
      computeAnalyticFlame << endl;
      computeAnalyticFlame << "cout << \"               Success !\" << endl << endl;" << endl;
      computeAnalyticFlame << endl;
      computeAnalyticFlame << "      int npoints = flow->nPoints();" << endl;
      computeAnalyticFlame << endl;
      computeAnalyticFlame << "      double position_max_wdot = 0;" << endl;
      computeAnalyticFlame << "      double max_wdot = 0;" << endl;
      computeAnalyticFlame << "      vector<double> max_wdot_k (nsp, 0.0);" << endl;
      computeAnalyticFlame << endl;
      computeAnalyticFlame << endl;
//      computeAnalyticFlame << "      string data = \"wdot.dat\";" << endl;
//      computeAnalyticFlame << "      ofstream dat(data.c_str());" << endl;
//      computeAnalyticFlame << "      dat << \"(32:wdot_CO) (33:wdot_NC10H22) (34:wdot_NO)\" << endl;"<< endl;
      computeAnalyticFlame << endl;
      computeAnalyticFlame << endl;
      computeAnalyticFlame << "      for (int ino=0; ino<npoints; ino++)" << endl;
      computeAnalyticFlame << "      {" << endl;
      computeAnalyticFlame << "         for (int k=0; k<nsp; k++)" << endl;
      computeAnalyticFlame << "         {" << endl;
      computeAnalyticFlame << "            if (fabs(max_wdot_k[k]) < fabs(flow->masswdot(k,ino)))" << endl;
      computeAnalyticFlame << "            {" << endl;
      computeAnalyticFlame << "               max_wdot_k[k] = flow->masswdot(k,ino);" << endl;
      computeAnalyticFlame << "            }" << endl;
      computeAnalyticFlame << endl;
      computeAnalyticFlame << "            if (max_wdot < fabs(flow->masswdot(k,ino)))" << endl;
      computeAnalyticFlame << "            {" << endl;
      computeAnalyticFlame << "               max_wdot = flow->masswdot(k,ino);" << endl;
      computeAnalyticFlame << "               position_max_wdot = flow->grid(ino);" << endl;
      computeAnalyticFlame << "            }" << endl;
//      computeAnalyticFlame << "            if (mixture->speciesName(k) == \"CO\")" << endl;
//      computeAnalyticFlame << "               dat << max_wdot_k[k] << \"  \";" << endl;
//      computeAnalyticFlame << "            if (mixture->speciesName(k) == \"NC10H22\")" << endl;
//      computeAnalyticFlame << "               dat << max_wdot_k[k] << \"  \";" << endl;
//      computeAnalyticFlame << "            if (mixture->speciesName(k) == \"NO\")" << endl;
//      computeAnalyticFlame << "               dat << max_wdot_k[k] << \"  \";" << endl;
      computeAnalyticFlame << "         }" << endl;
//      computeAnalyticFlame << "         dat << endl;" << endl;
      computeAnalyticFlame << "      }" << endl;
      computeAnalyticFlame << endl;
      computeAnalyticFlame << "      int ino_start = 0;" << endl;
      computeAnalyticFlame << "      int ino_stop = 0;" << endl;
      computeAnalyticFlame << endl;
      computeAnalyticFlame << "      for (int ino=0; ino<npoints; ino++)" << endl;
      computeAnalyticFlame << "      {" << endl;
      computeAnalyticFlame << "         for (int k=0; k<nsp; k++)" << endl;
      computeAnalyticFlame << "         {" << endl;
      computeAnalyticFlame << "            if (ino_start == 0 && fabs(flow->masswdot(k,ino))>0.05*fabs(max_wdot_k[k]))" << endl;
      computeAnalyticFlame << "            {" << endl;
      computeAnalyticFlame << "               ino_start = ino;" << endl;
      computeAnalyticFlame << "            }" << endl;
      computeAnalyticFlame << "         }" << endl;
      computeAnalyticFlame << "      }" << endl;
      computeAnalyticFlame << endl;
      computeAnalyticFlame << "      for (int ino=npoints-1; ino>0; ino--)" << endl;
      computeAnalyticFlame << "      {" << endl;
      computeAnalyticFlame << "         for (int k=0; k<nsp; k++)" << endl;
      computeAnalyticFlame << "         {" << endl;
      computeAnalyticFlame << "            if (ino_stop == 0 && fabs(flow->masswdot(k,ino)>0.05*fabs(max_wdot_k[k])))" << endl;
      computeAnalyticFlame << "            {" << endl;
      computeAnalyticFlame << "               ino_stop = ino;" << endl;
      computeAnalyticFlame << "            }" << endl;
      computeAnalyticFlame << "         }" << endl;
      computeAnalyticFlame << "      }" << endl;
      computeAnalyticFlame << endl;
      computeAnalyticFlame << "      vector<vector<double> > position (listFlames.size(), vector<double> (npoints, 0.0));" << endl;
      computeAnalyticFlame << "      vector<vector<vector<double> > > Y (listFlames.size(), vector<vector<double> > (npoints, vector<double> (nsp, 0.0)));" << endl;
      computeAnalyticFlame << "      vector<vector<double> > T (listFlames.size(), vector<double> (npoints, 0.0));" << endl;
      computeAnalyticFlame << "      vector<vector<double> > U (listFlames.size(), vector<double> (npoints, 0.0));" << endl;
      computeAnalyticFlame << endl;
      computeAnalyticFlame << "      for (int ino=0; ino<npoints; ino++)" << endl;
      computeAnalyticFlame << "      {" << endl;
      computeAnalyticFlame << "         position[n][ino] = flow->grid(ino);" << endl;
      computeAnalyticFlame << "         T[n][ino] = flame->value(flowdomain, flow->componentIndex(\"T\"), ino);" << endl;
      computeAnalyticFlame << "         U[n][ino] = flame->value(flowdomain, flow->componentIndex(\"u\"), ino);" << endl;
      computeAnalyticFlame << "         for (int k=0; k<nsp; k++)" << endl;
      computeAnalyticFlame << "         {" << endl;
      computeAnalyticFlame << "            string sp_name = mixture->speciesName(k);" << endl;
      computeAnalyticFlame << "            int sp_index = flow->componentIndex(sp_name);" << endl;
      computeAnalyticFlame << "            Y[n][ino][k] = flame->value(flowdomain, sp_index, ino);" << endl;
      computeAnalyticFlame << "         }" << endl;
      computeAnalyticFlame << "      }" << endl;
      computeAnalyticFlame << endl;
      computeAnalyticFlame << " " << endl;
      computeAnalyticFlame << " " << endl;
      computeAnalyticFlame << "      string redPath = \"Reduced_QSS\";" << endl;
      computeAnalyticFlame << "      stringstream s_nsp;" << endl;
      computeAnalyticFlame << "      s_nsp << nsp;" << endl;
      computeAnalyticFlame << "      stringstream s_n;" << endl;
      computeAnalyticFlame << "      s_n << n;" << endl;
      computeAnalyticFlame << "      redPath.append(s_nsp.str()).append(\"_\").append(s_n.str());" << endl;
      computeAnalyticFlame << " " << endl;
      computeAnalyticFlame << "      output_datas (redPath, Y[n], T[n], U[n], position[n], position_max_wdot, mixture);" << endl;
      computeAnalyticFlame << "   } //end loop nbFlames" << endl;
      computeAnalyticFlame << " " << endl;
      computeAnalyticFlame << "   int rank = 0; " << endl;
      computeAnalyticFlame << "   Script_gnuplot (\"QSS\", initial_mech, speciesToPlot, configuration, mech_desc, nsp , plot_U, plot_T, trajectory_ref, mech_ref, listFlames.size()-1, listTargets, mech, rank); " << endl;
      computeAnalyticFlame << "   string graph = \"gnuplot makeGnu.gnu\";" << endl;
      computeAnalyticFlame << "   int ret = system(graph.c_str());" << endl;
      computeAnalyticFlame << "   if (ret == -1) cout << \"ERROR system\";" << endl;
      computeAnalyticFlame << endl;
      computeAnalyticFlame << "}"   << endl;








   }


   if (configuration == "AutoIgnition")
   {
      string create_computeAnalyticFlame = fileName;
      create_computeAnalyticFlame.append("/computeAnalyticFlame.cpp");
      ofstream computeAnalyticFlame (create_computeAnalyticFlame.c_str());

      computeAnalyticFlame << "#include \"mech_QSS.h\"" << endl;
      computeAnalyticFlame << "" << endl;
      computeAnalyticFlame << "#include \"../../../../ORCh/main/conditions.h\"" << endl;
      computeAnalyticFlame << "" << endl;
      computeAnalyticFlame << "#include <Cantera.h>" << endl;
      computeAnalyticFlame << "#include <IdealGasMix.h>    // defines class IdealGasMix" << endl;
      computeAnalyticFlame << "#include <equilibrium.h>    // chemical equilibrium" << endl;
      computeAnalyticFlame << "#include <transport.h>      // transport properties" << endl;
      computeAnalyticFlame << "#include <zerodim.h>" << endl;
      computeAnalyticFlame << "#include <user.h>" << endl;
      computeAnalyticFlame << "#include <sstream>" << endl;
      computeAnalyticFlame << "" << endl;
      computeAnalyticFlame << "" << endl;
      computeAnalyticFlame << "using namespace Cantera;" << endl;
      computeAnalyticFlame << "using namespace Cantera_CXX;" << endl;
      computeAnalyticFlame << "using namespace User;" << endl;
      computeAnalyticFlame << "" << endl;
      computeAnalyticFlame << "IdealGasMix* mixture;" << endl;
      computeAnalyticFlame << "" << endl;
      computeAnalyticFlame << "int main()" << endl;
      computeAnalyticFlame << "{" << endl;
      computeAnalyticFlame << "   string mech_ref; string mech; string mech_desc; string configuration;" << endl;
      computeAnalyticFlame << "   vector<MultipleInlet*> listInlets;" << endl;
      computeAnalyticFlame << "   vector<PremixedFlames*> listFlames;" << endl;
      computeAnalyticFlame << "   vector<AutoIgnition*> listIgnitions;" << endl;
      computeAnalyticFlame << "   vector<string> listTargets;" << endl;
      computeAnalyticFlame << "   bool new_mixing;" << endl;
      computeAnalyticFlame << "   string step = \"\";" << endl;
      computeAnalyticFlame << "   vector<QSSscenario*> listQSSscenarios;" << endl;
      computeAnalyticFlame << "   OptimScenario* listOptimScenarios;" << endl;
      computeAnalyticFlame << "   int debuglevel;" << endl;
      computeAnalyticFlame << "   vector<string> speciesToPlot;" << endl;
      computeAnalyticFlame << "   bool plot_T;" << endl;
      computeAnalyticFlame << "   bool plot_U;" << endl;
      computeAnalyticFlame << "   string trajectory_ref;" << endl;
      computeAnalyticFlame << "" << endl;
      computeAnalyticFlame << "   conditions(debuglevel, speciesToPlot, mech_ref, mech, mech_desc, configuration, listInlets, listFlames, listIgnitions, listTargets, step, new_mixing, plot_T, plot_U, listQSSscenarios, listOptimScenarios, trajectory_ref);" << endl;
      computeAnalyticFlame << "" << endl;
      computeAnalyticFlame << "   mech = \"scheme.xml\";" << endl;
      computeAnalyticFlame << "" << endl;
      computeAnalyticFlame << "   mixture  = new IdealGasMix(mech, mech_desc);" << endl;
      computeAnalyticFlame << "   int nsp = mixture->nSpecies();"   << endl;
      computeAnalyticFlame << endl;
      computeAnalyticFlame << "   double *Ym = new double[nsp];" << endl;
      computeAnalyticFlame << "   for (int k=0; k<nsp; k++)" << endl;
      computeAnalyticFlame << "      Ym[k] = 0;" << endl;
      computeAnalyticFlame << "   double Hm = 0.0;"   << endl;
      computeAnalyticFlame << "   "   << endl;
      computeAnalyticFlame << "   "   << endl;
      computeAnalyticFlame << endl;


//temporary
     computeAnalyticFlame << "   ofstream i_delay (\"i_delay.dat\");" << endl; 
  
     computeAnalyticFlame << "   double t_phi_init = 0.001; "   << endl;
     computeAnalyticFlame << "   for (int nbPhi=1; nbPhi<105; nbPhi++)" << endl; 
     computeAnalyticFlame << "   {" << endl;
     computeAnalyticFlame << "     double t_phi = t_phi_init*pow(1.1,nbPhi);" << endl;
     computeAnalyticFlame << "     double Y_O2_fuel = 0.194221;" << endl;
     computeAnalyticFlame << "     double Y_N2_fuel = 0.589443;" << endl;
     computeAnalyticFlame << "     double Y_H2O_fuel = 0.00211403;" << endl;
     computeAnalyticFlame << "     double Y_CH4_fuel = 0.214222;" << endl;
     computeAnalyticFlame << "     double T_fuel = 320;" << endl;
     computeAnalyticFlame << endl;
     computeAnalyticFlame << "     double Y_O2_coflow = 0.142234;" << endl;
     computeAnalyticFlame << "     double Y_N2_coflow = 0.757491;" << endl;
     computeAnalyticFlame << "     double Y_H2O_coflow = 0.100097;" << endl;
     computeAnalyticFlame << "     double Y_CH4_coflow = 0.00017827;" << endl;
     computeAnalyticFlame << "     double T_coflow = 1350;" << endl;
     computeAnalyticFlame << endl;
     computeAnalyticFlame << "     double Z_1 = (4*Y_CH4_fuel-t_phi*Y_O2_fuel)/(t_phi*Y_O2_coflow-4*Y_CH4_coflow)/(1+(4*Y_CH4_fuel-t_phi*Y_O2_fuel)/(t_phi*Y_O2_coflow-4*Y_CH4_coflow));" << endl;
     computeAnalyticFlame << "     double Z_2 = 1-Z_1;" << endl;
     computeAnalyticFlame << endl;
     computeAnalyticFlame << "     double Y_O2 = Y_O2_fuel*Z_2+Y_O2_coflow*Z_1;" << endl;
     computeAnalyticFlame << "     double Y_H2O = Y_H2O_fuel*Z_2+Y_H2O_coflow*Z_1;" << endl;
     computeAnalyticFlame << "     double Y_CH4 = Y_CH4_fuel*Z_2+Y_CH4_coflow*Z_1;" << endl;
     computeAnalyticFlame << "     double Y_N2 = Y_N2_fuel*Z_2+Y_N2_coflow*Z_1;" << endl;
     computeAnalyticFlame << endl;
     computeAnalyticFlame << "" << endl;
     computeAnalyticFlame << "      for (int k=0; k<nsp; k++)" << endl;
     computeAnalyticFlame << "      {" << endl;
     computeAnalyticFlame << "         if (mixture->speciesName(k) == \"CH4\")" << endl;
     computeAnalyticFlame << "            Ym[k] = Y_CH4;" << endl;
     computeAnalyticFlame << endl;
     computeAnalyticFlame << "         if (mixture->speciesName(k) == \"O2\")" << endl;
     computeAnalyticFlame << "            Ym[k] = Y_O2;" << endl;
     computeAnalyticFlame << endl;
     computeAnalyticFlame << "         if (mixture->speciesName(k) == \"H2O\")" << endl;
     computeAnalyticFlame << "            Ym[k] = Y_H2O;" << endl;
     computeAnalyticFlame << endl;
     computeAnalyticFlame << "         if (mixture->speciesName(k) == \"N2\")" << endl;
     computeAnalyticFlame << "            Ym[k] = Y_N2;" << endl;
     computeAnalyticFlame << "      }" << endl;
     computeAnalyticFlame << "      cout << \"Z1 \" << Z_1 << endl;" << endl;
     computeAnalyticFlame << "      cout << \"Z2 \" << Z_2 << endl;" << endl;
     computeAnalyticFlame << endl;
     computeAnalyticFlame << "      double Tm = T_fuel*Z_2+T_coflow*Z_1;" << endl;

    computeAnalyticFlame << "   for (int n=0; n<listIgnitions.size(); n++) " << endl;
    computeAnalyticFlame << "   { " << endl;
    computeAnalyticFlame << "      int nbLines = listIgnitions[n]->m_max_t/listIgnitions[n]->m_delta_t; " << endl;
    computeAnalyticFlame << "      cout << \"nbLines \" << nbLines << endl; " << endl;
    computeAnalyticFlame << endl;
    computeAnalyticFlame << "      double *Y = new double[nsp]; " << endl;
    computeAnalyticFlame << "      double H = 0.0; " << endl;
    computeAnalyticFlame << "      double T = 0.0; " << endl;
    computeAnalyticFlame << endl;
    computeAnalyticFlame << "      //---Particles store--- " << endl;
    computeAnalyticFlame << "      vector<vector<double> > Y_Store(nbLines, vector<double>(nsp)); " << endl;
    computeAnalyticFlame << "      vector<double> H_Store(nbLines, 0.0); " << endl;
    computeAnalyticFlame << "      vector<double> T_Store(nbLines, 0.0); " << endl;
    computeAnalyticFlame << "      vector<double> t_Store(nbLines, 0.0); " << endl;
    computeAnalyticFlame << endl;
    computeAnalyticFlame << "      IdealGasMix *IgnitionMixture = new IdealGasMix(mech, mech_desc); " << endl;
    computeAnalyticFlame << "      int nsp = IgnitionMixture->nSpecies(); " << endl;
    computeAnalyticFlame << "      double t = listIgnitions[n]->m_delta_t; " << endl;
    computeAnalyticFlame << "      IgnitionMixture->setMassFractions(Ym); " << endl;
    computeAnalyticFlame << endl;
    computeAnalyticFlame << "      IgnitionMixture->setState_TP(Tm, listIgnitions[n]->m_Pressure); " << endl;
    computeAnalyticFlame << endl;
    computeAnalyticFlame << "      t = listIgnitions[n]->m_delta_t; " << endl;
    computeAnalyticFlame << endl;
    computeAnalyticFlame << "      string path = \"Ignite\"; " << endl;
    computeAnalyticFlame << "      stringstream s_Phi; " << endl;
    computeAnalyticFlame << "      s_Phi << t_phi; " << endl;
    computeAnalyticFlame << "      path.append(s_Phi.str()).append(\"_\").append(\"13S_3QSS.dat\");" << endl;
    computeAnalyticFlame << "      ofstream Ignite (path.c_str()); " << endl;
    computeAnalyticFlame << "" << endl;
    computeAnalyticFlame << "      for (int i=0; i<nbLines-1; i++) " << endl;
    computeAnalyticFlame << "      { " << endl;
    computeAnalyticFlame << "         ConstPressureReactor reac; " << endl;
    computeAnalyticFlame << "         reac.insert(*IgnitionMixture); " << endl;
    computeAnalyticFlame << "         ReactorNet sim; " << endl;
    computeAnalyticFlame << "         sim.addReactor(&reac,false); " << endl;
    computeAnalyticFlame << "         sim.advance(listIgnitions[n]->m_delta_t); " << endl;
    computeAnalyticFlame << "" << endl;
    computeAnalyticFlame << "         H  = IgnitionMixture->enthalpy_mass(); " << endl;
    computeAnalyticFlame << "         T = IgnitionMixture->temperature(); " << endl;
    computeAnalyticFlame << "         IgnitionMixture->getMassFractions(Y); " << endl;
    computeAnalyticFlame << "" << endl;
    computeAnalyticFlame << "         Ignite << t << \"  \" << T << endl; " << endl;
    computeAnalyticFlame << "" << endl;
    computeAnalyticFlame << "         t = t + listIgnitions[n]->m_delta_t; " << endl;
    computeAnalyticFlame << "         t_Store[i+1] = t; " << endl;
    computeAnalyticFlame << "         T_Store[i+1] = T; " << endl;
    computeAnalyticFlame << "      }" << endl;
    computeAnalyticFlame << "      bool done = true; " << endl;
    computeAnalyticFlame << "      for (int i=0; i<nbLines; i++) " << endl;
    computeAnalyticFlame << "      { " << endl;
    computeAnalyticFlame << "         if (T_Store[i] > T_Store[1]+0.5*(T_Store[nbLines-1]-T_Store[1]) && done) " << endl;
    computeAnalyticFlame << "         { " << endl;
    computeAnalyticFlame << "            cout << \"Ignition delay \" << t_Store[i] << \"  phi \" << t_phi << \" T_middle \" << T_Store[i] << \"  Tmin\" << T_Store[0] << \"  Tmax\" << T_Store[nbLines-1] << endl; " << endl;
    computeAnalyticFlame << "            i_delay << t_phi << \"  \" << t_Store[i] << endl; " << endl;
    computeAnalyticFlame << "            done = false; " << endl;
    computeAnalyticFlame << "         } " << endl;
    computeAnalyticFlame << "      } " << endl;
    computeAnalyticFlame << "   }" << endl;

      computeAnalyticFlame << endl;
      computeAnalyticFlame << endl;
      computeAnalyticFlame << endl;
      computeAnalyticFlame << endl;
      computeAnalyticFlame << endl;
      computeAnalyticFlame << "   }" << endl;
      computeAnalyticFlame << "}" << endl;

   }





   if (configuration == "MultipleInlet" || configuration == "AutoIgnition")
   {

      string create_Makefile = fileName;
      create_Makefile.append("/Makefile");
      ofstream Makefile (create_Makefile.c_str());

      Makefile << "" << endl;
      Makefile << "ORCh= ../../../../ORCh" << endl;
      Makefile << "" << endl;
      Makefile << "all: computeAnalyticFlame" << endl;
      Makefile << "" << endl;
      Makefile << "computeAnalyticFlame: computeAnalyticFlame.cpp" << endl;
      Makefile << "\t@echo \"Compiling files ...\"" << endl;
      Makefile << "\t@mpicxx -w -std=c++11 -L$(MPI_LIB) -I$(BOOSTPATH) -I$(GTCOMB_CT_HOME)/include -L$(GTCOMB_CT_HOME)/build/lib -L$(GTCOMB_CT_HOSTTYPE)  -I$(HDF5_INC) ../../conditions.cpp $(ORCh)/read_write/read.cpp $(ORCh)/read_write/reaction.cpp $(ORCh)/tools/outputs.cpp $(ORCh)/tools/gnuplot.cpp  $(ORCh)/tools/fitness_criteria.cpp  $(ORCh)/read_write/species.cpp $(ORCh)/cantera/computeMultipleInlet.cpp $(ORCh)/cantera/computePremixedFlames.cpp $(ORCh)/cantera/computeAutoIgnition.cpp $(ORCh)/cantera/particle.cpp $(ORCh)/cantera/flamemodel.cpp $(ORCh)/drgep/drgep.cpp $(ORCh)/read_write/QSSscenario.cpp $(ORCh)/optimisation/OptimScenario.cpp $(ORCh)/emst/emst_subs.o $(ORCh)/emst/emst.o -o computeAnalyticFlame computeAnalyticFlame.cpp -O3 -Wall -lcantera  -luser -lm -lhdf5 -lhdf5_hl -lgfortran" << endl;
      Makefile << "" << endl;
      Makefile << "clean:" << endl;
      Makefile << "\t@rm -f computeAnalyticFlame" << endl;
      Makefile.close();
   }
   
   if (configuration == "PremixedFlames")
   {
      string create_Makefile = fileName;
      create_Makefile.append("/Makefile");
      ofstream Makefile (create_Makefile.c_str());

      Makefile << "" << endl;
      Makefile << "ORCh= ../../../../ORCh" << endl;
      Makefile << "" << endl;
      Makefile << "all: computeAnalyticFlame" << endl;
      Makefile << "computeAnalyticFlame: computeAnalyticFlame.cpp" << endl;
      Makefile << "\t@echo \"Compiling files ...\"" << endl;
      Makefile << "\t@mpicxx -w -std=c++11 -L$(MPI_LIB) -I$(BOOSTPATH) -I$(GTCOMB_CT_HOME)/include -L$(GTCOMB_CT_HOME)/build/lib -L$(GTCOMB_CT_HOSTTYPE)  -I$(HDF5_INC)  ../../conditions.cpp $(ORCh)/read_write/read.cpp $(ORCh)/read_write/reaction.cpp $(ORCh)/tools/outputs.cpp $(ORCh)/tools/gnuplot.cpp  $(ORCh)/tools/fitness_criteria.cpp  $(ORCh)/read_write/species.cpp  $(ORCh)/cantera/computePremixedFlames.cpp $(ORCh)/cantera/particle.cpp $(ORCh)/cantera/flamemodel.cpp $(ORCh)/drgep/drgep.cpp $(ORCh)/read_write/QSSscenario.cpp $(ORCh)/optimisation/OptimScenario.cpp $(ORCh)/emst/emst_subs.o $(ORCh)/emst/emst.o -o computeAnalyticFlame computeAnalyticFlame.cpp -O3 -Wall -lcantera -luser -lm -lhdf5 -lhdf5_hl -lgfortran" << endl;
      Makefile << "" << endl;
      Makefile << "clean:" << endl;
      Makefile << "\t@rm -f computeAnalyticFlame" << endl;
      Makefile.close();
   }


}







