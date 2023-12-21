#include "computeMultipleInlet.h"
#include <fstream>	//Huu-Tri: "read from file" library
#include <iostream>	//Huu-Tri
#include <sys/stat.h> 	//Huu-Tri: Check file exist - stat function - 20 Nov 2019
#include <tensorflow/c/c_api.h> //Huu-Tri@20200723 : Add Tensorflow C API to use trained model from Python Keras

// Huu-Tri@20200724 : Add Tensorflow libraries "cppflow" - Use Tensorflow C-API to load ANN model and predict in C++
// /home/huutri/workdir/orch/ORCh/cppflow
// CppFlow : https://github.com/serizba/cppflow
#include "../cppflow/include/Model.h"
#include "../cppflow/include/Tensor.h"
#include "opencv2/core/core.hpp" // Lib for PCA
#include <numeric>
#include <iomanip>

#include <fstream>// Huu-Tri NGUYEN 20210212 - To read from text file
#include <iostream>
#include <sstream>
#include <string>
using namespace std;



/* EMST model */
extern "C" void emst_(int* mode_emst,int* np_emst,int* nc_emst,double* f_emst,double* state_emst,double* wt_emst,double* omdt_emst,double* fscale_emst,double* cvars_emst,int* info_emst); // EMST mixing model - Edited by Kaidi - Added by Huu-Tri Nguyen 10.12.2019


//---computeMultipleInlet---

computeMultipleInlet::computeMultipleInlet() //Constructeur
{}

double Get_random ()
{
    int random = rand() % 10000;
    double random_number = double(random)/10000.0;
    return random_number;
}

/* Check if file exists */
/* Return true if the file exists */
/* Huu-Tri Nguyen 20 Nov 2019 */
bool file_exists(char *filename)
{
	struct stat buffer;
	return(stat(filename, &buffer) == 0);	//if the file exits, return true
}


/* Print to file - Source term Deterministic */
void SaveToFile_wDeter(double *wDeter_T, double **wDeter_species, int nsp, int nbInlets, int i, double delta_t, int rank, IdealGasMix *mixture)
{
// Print the Deterministic source term to CFD_results/sourceTermDeter_Inlet*.txt

	if(rank==0)	//Parallel MPI stuff to prevent print multiple lines, rank = 0 <=> First processor
	{

		// Check if sourceTermDeter_Inlet*.txt file exists at the first step
		// If yes, clear the file content
		if((file_exists("CFD_results/sourceTermDeter_Inlet0.txt") || file_exists("CFD_results/sourceTermDeter_Inlet1.txt")) && i==0)
		{
		
			cout << " -------------- Warrning --------------" << endl;
			cout << "CFD_results/sourceTermDeter_Inlet*.txt exists. Clearing file ... " << endl;	
			
			ofstream sourceTermDeter_Inlet0_clear("CFD_results/sourceTermDeter_Inlet0.txt", ios::out | ios::trunc);	//open file in trunc mode to clear the content
			sourceTermDeter_Inlet0_clear.close(); //close the file

			ofstream sourceTermDeter_Inlet1_clear("CFD_results/sourceTermDeter_Inlet1.txt", ios::out | ios::trunc);	//open file in trunc mode to clear the content
			sourceTermDeter_Inlet1_clear.close(); //close the file
	
		}
		

		
		// Inlet 0
		ofstream sourceTermDeter_Inlet0("CFD_results/sourceTermDeter_Inlet0.txt",ios::app); //ios::app = append at the end of the file
		if(sourceTermDeter_Inlet0)
		{
			if(i==0)	//First step: Need to write the headline
			{
				// First line
				sourceTermDeter_Inlet0 << "#1:Time(s)	"; 
				sourceTermDeter_Inlet0 << "2:T	";
				for (int k=0; k<nsp; k++)
				{
					sourceTermDeter_Inlet0 << k+3 << ":" << mixture->speciesName(k) << "	";
				}
				sourceTermDeter_Inlet0 << endl;			

				// Data from ORCh
				sourceTermDeter_Inlet0 << i*delta_t << "	";
				sourceTermDeter_Inlet0 << wDeter_T[0] << "	";
				for (int k=0; k<nsp; k++)
				{
					sourceTermDeter_Inlet0 << wDeter_species[0][k] << "	";
				}
				sourceTermDeter_Inlet0 << endl;
			}
			else
			{		
				// Data from ORCh
				sourceTermDeter_Inlet0 << i*delta_t << "	";
				sourceTermDeter_Inlet0 << wDeter_T[0] << "	";
				for (int k=0; k<nsp; k++)
				{
					sourceTermDeter_Inlet0 << wDeter_species[0][k]  << "	";
					}
				sourceTermDeter_Inlet0 << endl;
			}
		}
		else
		{	
			cout << "ERROR: Impossible to write sourceTermDeter_Inlet0.txt" << endl;
			cout << "Please check computeMultipleInlet.cpp" << endl;
		}
	
		sourceTermDeter_Inlet0.close();

		// Inlet 1
		ofstream sourceTermDeter_Inlet1("CFD_results/sourceTermDeter_Inlet1.txt",ios::app); //ios::app = append at the end of the file
		if(sourceTermDeter_Inlet1)
		{
			if(i==0)	//First step: Need to write the headline
			{
				// First line
				sourceTermDeter_Inlet1 << "#1:Time(s)	"; 
				sourceTermDeter_Inlet1 << "2:T	";
				for (int k=0; k<nsp; k++)
				{
					sourceTermDeter_Inlet1 << k+3 << ":" << mixture->speciesName(k) << "	";
				}
				sourceTermDeter_Inlet1 << endl;			

				// Data from ORCh
				sourceTermDeter_Inlet1 << i*delta_t << "	";
				sourceTermDeter_Inlet1 << wDeter_T[1] << "	";
				for (int k=0; k<nsp; k++)
				{
					sourceTermDeter_Inlet1 << wDeter_species[1][k] << "	";
				}
				sourceTermDeter_Inlet1 << endl;
			}
			else
			{		
				// Data from ORCh
				sourceTermDeter_Inlet1 << i*delta_t << "	";
				sourceTermDeter_Inlet1 << wDeter_T[1] << "	";
				for (int k=0; k<nsp; k++)
				{
					sourceTermDeter_Inlet1 << wDeter_species[1][k]  << "	";
					}
				sourceTermDeter_Inlet1 << endl;
			}
		}
		else
		{	
			cout << "ERROR: Impossible to write sourceTermDeter_Inlet1.txt" << endl;
			cout << "Please check computeMultipleInlet.cpp" << endl;
		}
	
		sourceTermDeter_Inlet0.close();

	} // End if(rank==0)
} // END SaveToFile_wDeter



// Declare function for ANN - Huu-Tri 20210212

// END Declare function for ANN - Huu-Tri 20210212


void computeMultipleInlet::getMultipleInlet(
   string mech,
   string mech_desc,
   vector<MultipleInlet*> listInlets,
   vector<bool> Targets,
   bool new_mixing,
   string step,
   vector<vector<vector<double> > > &R_AD_Trajectories,
   vector<vector<double> > &max_j_on_Target,
   vector<vector<vector<double> > > &Ym_Trajectories_store,
   vector<vector<vector<double> > > &Production_Trajectories_ref,
   vector<vector<vector<double> > > &Consumption_Trajectories_ref,
   vector<vector<double> > &T_Trajectories_store,
   vector<double> &time_store,
   vector<bool> &SpeciesIntoReactants)
{

   mixture  = new IdealGasMix(mech,mech_desc);
   int nsp = mixture->nSpecies();
   int nbInlets = listInlets.size();

   getMixedGasesComposition(listInlets, step);

   double t = 0.0; //Initial computational time
   int nTot = 0; //Total number of particles

   //Treat parallel stuff
   int rank, nproc;

   	if (step != "Optimisation") {
		
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		MPI_Comm_size(MPI_COMM_WORLD, &nproc);
        }
	
		
      else
          rank = 0;

   double Particle_flowRate = 0.001; //elementary mass flow rate
   vector<int> nbParticles (nbInlets, 0); //Number of particles per inlet
   for (int n=0; n<nbInlets; n++)
   {
      nbParticles[n] = listInlets[n]->m_flowRate/Particle_flowRate;
      if (rank == 0)
      {
         cout << "Nb particles  " << n << "  " << nbParticles[n] << endl;
      }
      nTot += nbParticles[n];
   }

  //Treat parallel stuff

   	if (step != "Optimisation") {
		//
		// DAK : initialize Ifi_rank and Ila_rank (particles to be computed by this proc)
		//
		
		int Ifi = 0;
		int Ila = 0;
		RecvCounts = new int[nproc];
		Disp = new int[nproc];
		for (int r=0; r<nproc; r++) Disp[r] = 0;
		
		for (int r=0; r<nproc; r++)
		{
			if (r < nTot%nproc)
			{
				Ifi = (nTot/nproc)*r + r;
				Ila = (nTot/nproc)*r + (nTot/nproc) + r + 1;
			}
			else
			{
				Ifi = (nTot/nproc)*r + nTot%nproc;
				Ila = (nTot/nproc)*r + (nTot/nproc) + nTot%nproc;
			}
			
			RecvCounts[r] = (Ila-Ifi)*(nsp+2); //for Yks and H and T
			for (int rb=r; rb<nproc-1; rb++) Disp[rb+1] += (Ila-Ifi)*(nsp+2);
			
			if (rank == r)
			{
				Ifi_rank = Ifi;
				Ila_rank = Ila;
				nb_var_loc = (Ila-Ifi)*(nsp+2);
			}
		}
	}

	//Huu-Tri - 13.01.2020
	if(rank==0) cout << "nTot = " << nTot << endl;
   int AirPart = ceil(nbParticles[0]*0.4);

   //-----------------------------------------------------------------------------------------------------------------------//
   //   Randomly select the particles that will be mixed or read that into the "Selection.dat" file if new_mixing == false
   //-----------------------------------------------------------------------------------------------------------------------//

   vector<vector<int> > Particle_1;
   vector<vector<int> > Particle_2;

   double delta_t = dynamic_cast <Characteristics_MultipleInlet*> (listInlets[nbInlets-1])->m_delta_t;
   double tau_t = dynamic_cast <Characteristics_MultipleInlet*> (listInlets[nbInlets-1])->m_tau_t;
   int nbLines =  dynamic_cast <Characteristics_MultipleInlet*> (listInlets[nbInlets-1])->m_NbIterations;
	

   // Huu-Tri Nguyen - 20.01.2020
	if(rank==0) cout << "Time step = " << delta_t << " | Iterations = " << nbLines << " | Mixing time = " << tau_t << endl;

   double Pressure =  dynamic_cast <Characteristics_MultipleInlet*> (listInlets[nbInlets-1])->m_Pressure;
   double F =1;

   int Nmix = delta_t*nTot/tau_t;
   if (rank == 0)
      cout << "Nmix " << Nmix << endl;

   if (Nmix < 0) 
   {
      cout << "Problem with the delta_t and tau_t description, Nmix = " << Nmix << endl;
      getchar();
   }

   // by Kaidi@2019.11 - Huu-Tri Nguyen added - 10.12.2019
   if (Nmix*2 > nTot)
   {
      cout << "Problem with the delta_t and tau_t description, Nmix is too large!" << endl;
      getchar();
   }


   //seed the random number generator
   srand(time(NULL));

   if (new_mixing == false || step == "Optimisation_XML" || step == "Optimisation_Analytical")
   {
      ifstream Selection_read("Selection.dat");

      for (int nb=0; nb<nbLines; nb++)
      {
         Particle_1.push_back(vector<int> (Nmix));
         Particle_2.push_back(vector<int> (Nmix));

         double a;
         Selection_read >> a;
         for (int sp=0; sp<Nmix; sp++)
         {
            Selection_read >> Particle_1[nb][sp] >> Particle_2[nb][sp];
         }
      }
      Selection_read.close();
   }
   else
   {
      ofstream Selection("Selection.dat");

      int select;

      for (int nb=0; nb<nbLines; nb++)
      {


         Particle_1.push_back(vector<int> (Nmix));
         Particle_2.push_back(vector<int> (Nmix));

         Selection << nb << "  ";


         for (int sp=0; sp<Nmix; sp++)
         {
            //Select first and second particle to be mixed
            for (int fs=0; fs<2; fs++)
            {
               bool particle_found = false;
               while (!particle_found)
               {
//air dilution add 17/07/17 
       //           if (nb < nbLines/2) // On ne mélange que 60% de l'air entrant avec le reste
       //           {
       //  	     select = rand() % (nTot-AirPart) + AirPart; 
       //                
       //           }
       //           else // on rajoute les particules d'air restantes petit à petit jusqu à la fin
       //           {
       //              float c = nb;
       //              float d = nbLines;
       //              double b = abs(AirPart*(1 - ( 2*( c - d/2 )/d)));
       //              int a = b;

       // 	     select = rand() % (nTot-a) + a;
       // 
       //           }
                 select = rand() % nTot;

                  if (sp == 0 && fs == 0)
                     particle_found = true;
                  if (sp == 0 && fs == 1)
                  {
                     if (select != Particle_1[nb][0])
                        particle_found = true;
                  }

                  //Check that the particle hasn't been used yet

             	  /* Comment to add the edited code from Kaidi - Huu-Tri Nguyen 10.12.2019 */
	//HT	  for (int sp_test=0; sp_test<sp; sp_test++)
        //HT      {
        //HT             if (select != Particle_1[nb][sp_test] && select != Particle_2[nb][sp_test])
        //HT             {
        //HT                particle_found = true;
        //HT             }
        //HT      }
		/* End comment - Huu-Tri Nguyen 10.12.2019 */


		/* Add the corrected code from Kaidi -  Huu-Tri Nguyen 10.12.2019 */
		bool particle_used = false;
                for (int sp_test=0; sp_test<sp; sp_test++)
                {
                	if (select == Particle_1[nb][sp_test] || select == Particle_2[nb][sp_test])
                     	{
                       		particle_used = true;
                     	}
                }
                
		if (fs == 1)
                {
                	if (select == Particle_1[nb][sp])  
				particle_used = true;
                }
                
		if (particle_used)
                	particle_found = false;
                else
                	particle_found = true;
		/* End the corrected code from Kaidi -  Huu-Tri Nguyen 10.12.2019 */

               }

               if (fs == 0)
                  Particle_1[nb][sp] = select;
               else if (fs == 1)
                  Particle_2[nb][sp] = select;
            }
            Selection << Particle_1[nb][sp] << " " << Particle_2[nb][sp] << "   ";
         }
         Selection << endl;

      }
      Selection.close();
   }


   //---Trajectories store---
   vector<double> Hm_Trajectories(nbInlets, 0.0);
   vector<double> Zm_Trajectories(nbInlets, 0.0);
   vector<double> Hm_inletIni(nbInlets, 0.0); // Save initial enthalpy of each inlet - Huu-Tri Nguyen - 16.01.2020

   vector<double> Ym(nsp,0.0);
   double Hm = 0.0;
   double Zm = 0.0;
   double Tm = 0.0;
   double density = 0.0;

   //First create the Particles which will transport the gaseous and liquid phases
   vector<Particle*> listParticles;








   double Diameter_init; 
   double tau_vj;
   for (int n=0; n<nbInlets; n++)
   {      
      if (listInlets[n]->m_X_Species != "")
      {
         if (rank == 0)
         {
         cout << "Set the mole fraction of inlet " << n << endl;
         }
         mixture->setState_TPX(listInlets[n]->m_Temperature, listInlets[n]->m_Pressure, listInlets[n]->m_X_Species);
      }
      else if (listInlets[n]->m_Y_Species != "")
      {
         if (rank == 0)
         {
            cout << "Set the mass fraction of inlet " << n << endl;
         }
         mixture->setState_TPY(listInlets[n]->m_Temperature, listInlets[n]->m_Pressure, listInlets[n]->m_Y_Species);
      }

      mixture->getMassFractions(&Ym[0]);
      Hm  = mixture->enthalpy_mass();
      vector<double> Y_transfer (nsp, 0.0);
      for (int k=0; k<nsp; k++)
         Y_transfer[k] = Ym[k];
      density = mixture->density();

      for (int k=0; k<nsp; k++)
      {
         if (Ym[k] > 0.0)
         {
            SpeciesIntoReactants[k] = true;
         }
      }

      
	Hm_inletIni[n] = Hm; // Save initial enthalpy of each inlet - Huu-Tri Nguyen - 16.01.2020
	if(rank ==0) cout << "Hm initial inlet " << n << " = " << Hm_inletIni[n] << endl;

      if (n < nbInlets-1)  
      {
         //Composition space Lagrangian trajectories
         for (int k=0; k<nsp; k++)
            Ym_Trajectories_store[n][0][k] = Y_transfer[k];
 
         Hm_Trajectories[n] = Hm;
         T_Trajectories_store[n][0] = listInlets[n]->m_Temperature;

      }
	
      

      for (int i=0; i<nbParticles[n]; i++)
      {
         if (listInlets[n]->m_liquid)
         {
            double nbDroplets = Particle_flowRate/(listInlets[n]->m_density_liquid*(PI/6)*pow(listInlets[n]->m_DropletsDiameter,3.0));
            listParticles.push_back(new Particle(
                           Y_transfer, 
                           listInlets[n]->m_Temperature, 
                           Hm, 
                           1, 
                           0.0, 
                           nbDroplets, 
                           listInlets[n]->m_DropletsDiameter, 
                           0.001 /*0% gas, 100% liquid*/, 
			   Y_transfer, 
                           listInlets[n]->m_density_liquid, 
                           listInlets[n]->m_EvaporationLatentHeat));

            Diameter_init = listInlets[n]->m_DropletsDiameter;
            tau_vj = listInlets[n]->m_Tau_vj;
         }
         else
         {
            listParticles.push_back(new Particle(
                           Y_transfer, 
                           listInlets[n]->m_Temperature, 
                           Hm, 
                           0, 
                           density, 
                           0, 
                           0.0, 
                           1.0 /*100% gas, 0% liquid*/, 
                           vector<double> (nsp, 0.0), 
                           0.0, 
                           0.0));
         }
      }
   }





   //Table with the sensible enthalpy of species i at Tboiling
   vector<double> BoilingEnthalpy (nsp, 0.0);
   for (int k=0; k<nsp; k++)
   {
      double *Ymb = new double[nsp];
      for (int kbis=0; kbis<nsp; kbis++)
      {
         if (k != kbis)
            Ymb[kbis] = 0;
         else
            Ymb[kbis] = 1;
      }
      mixture->setState_TPY(listInlets[1]->m_Temperature, listInlets[1]->m_Pressure, Ymb);
      BoilingEnthalpy[k] = mixture->enthalpy_mass();
      delete[] Ymb;
   }





   double dt = dynamic_cast <Characteristics_MultipleInlet*> (listInlets[nbInlets-1])->m_delta_t;
   double gas_mass_p1;
   double gas_mass_p2;
   double Total_gas_mass;


   vector<double> Mean_Ym(nsp, 0.0);
   double Mean_Hm = 0.0;
   double Mean_Zm = 0.0;
   double Mean_Tm = 0.0;
	
   // Data for scatterplot
   // ofstream store ("outputs/data.dat"); // Store mean values
   // ofstream store_particles ("data_particles.dat");	// Store scatterplot
   // store_particles << "#1:time  2:Particle_number  3:Zfraction  4:T(K) 5:Y_CO2 6:Y_O2 7:particle type  8:ratio  9:Zst " << endl;
	

   /* Add EMST model initialization - Huu-Tri Nguyen - 10.12.2019 */
   // EMST mixing model initialization -- by Kaidi@2019.12
   int mode_emst = 1;		// 1:-the state variables of all particles are initialized. No mixing is performed (check emst.f)
   int np_emst = nTot;		// np_emst: number of particles
   int nc_emst = nsp+1;		// nc_emst: number of compositions variables (+1 to add Enthalpy at the end of array)
   double *f_emst = new double[np_emst*nc_emst];	// Particle composition
   double *state_emst = new double[np_emst];		// State (or age) variable
   double *wt_emst = new double[np_emst];		// Particle weights. wt(i)>0 is the numerical weight of the i-th particle
   double *fscale_emst = new double[nc_emst];		// Scale factors for the compositions (fscale(j)>0) - see explanation in emst.f
   double cvars_emst[6] = {0,0,0,0,0,0};		// Control variables for expert users. cvars_emst[6] = {0,0,0,0,0,0} is default
   int info_emst = 0;					// = 0 for successful execution; < 0 for error in input
   int i_emst;
   double Cphi = 2; 					// Mixing model constant, see explanation in emst.f
   double omdt_emst = delta_t/tau_t*Cphi;		// Normalized mixing time 
	if(rank ==0) cout << "omdt = " << omdt_emst << endl;
   i_emst = 0;
   for (int ic=0; ic<nc_emst; ic++)	// ic: running variable for number of composition nc_emst
   {
      for (int ip=0; ip<np_emst; ip++)	// ik: running variable for number of particle np_emst
      {
         if (ic<nsp)  
		f_emst[i_emst] = listParticles[ip]->m_Yk_gas[ic];
         else         
		f_emst[i_emst] = listParticles[ip]->m_H_gas;		// Ad enthalpy in the end of array to calculate
         i_emst++;
      }
   }

   for (int ip=0; ip<np_emst; ip++)
   {
      state_emst[ip] = 0;
      wt_emst[ip] = 1.0;
   }

   for (int ic=0; ic<nc_emst; ic++)
   {
      if (ic<nsp)   fscale_emst[ic] = 0.1;	// Intialization values - recommended by emst.f
      else          fscale_emst[ic] = 1e16;
   }

   emst_(&mode_emst,&np_emst,&nc_emst,f_emst,state_emst,wt_emst,&omdt_emst,fscale_emst,cvars_emst,&info_emst);
   if (info_emst != 0)
   {
      cout << "emst initialization failed" << endl;
      getchar();
   }
   /* END Add EMST model initialization - Huu-Tri Nguyen - 10.12.2019 */

   int ndil = 0;

    // ============================================================================
    // Declare flagANN parameters before loop- Huu-Tri Nguyen - 20210206 =========
    // ============================================================================


    // END Declare ANN parametes  ======================
    // =================================================
    
    
    
   /* =============================== */
   /* BEGIN BIG LOOP - Each time step */
   /* =============================== */
   for (int i=0; i<nbLines-1; i++)
   {
    /* ======================================================== */
    /* Read enthalpy CFD file - Huu-Tri NGUYEN - 30 August 2019 */

    int maxRow = 0; // Max row
    int maxCol = 0; // Max column
    vector<double> t_CFD;		// Use vector because it can easily resize
    vector<double> mean_hCFD;
    


    // Path file - Huu-Tri NGUYEN - 2020.01.10
    string dir ="/home/huutri/workdir/orch/Workcases/";
    string Case = "121_EMST_ORChData.FourUmons_NewVersionORChANN20210312_Stochastic_Qx10_NewConditions_outletGB_ANNmodels_SansPerte/"; // Be careful, add / at the end

    string inputCFD = "CFD_results/20190918_h_tau.txt";
    string dirCFDIn = dir + Case + inputCFD;

    string outputCFD ="CFD_results/Output_h_tau_Interpolated.txt";
    string dirCFDOut = dir + Case + outputCFD;

     bool flagCFD = false;		// flagCFD == true, use CFD result correction
     bool flagCFDDeter = false;		// flagCFDDeter == true, use CFD result correction for Deterministic closure
     bool flagCFDStochas = false;	// flagCFDStochas == true, use CFD result correction for Stochastic closure

    // HuuTri@20220107: flagHeatLossStochas == true, use a sink/source term to impose the heat loss (alpha_loss) for Stochastic closure
	// Sink term  = alpha_loss*(T_gas - T_wall)
	// Heat term  = beta_heat*(T_gas - T_wall)
	// The value of alpha_loss, beta_heat also depends on the time step delta_t
	// To calculate enthalpy: It should be multiplied by delta_t: 
		// H_gasAfterLoss = H_gas -  sinkTerm*delta_t =  H_ini - alpha_loss*(T_gas - T_wall)*delta_t
		// H_gasAfterHeated = H_gas -  sourceTerm*delta_t =  H_ini - beta_heat*(T_gas - T_wall)*delta_t
     bool flagHeatLossStochas = true;
     double T_wall = 1600;
     double alpha_loss = 1.5e+05;	//1.5e+05 // 3e+05 //4.5e+05 //1.5e+06 //0.1e+05;
     double beta_heat = 7.5e+05;
     if(flagHeatLossStochas)
     {
	if(rank==0 && i==0)
	{
	cout << "==== Heat loss sink/source term correction ====" << endl;
	cout << "alpha_loss = " << alpha_loss << " | beta_heat = " << beta_heat << " | Twall = " << T_wall << endl;
	}
     }

     if (flagCFD)
     {
	if(rank==0 && i==0) cout << "==== With CFD correction ====" << endl;

//	ifstream enthalpyCFDfile(fullPathCorrectionIn.c_str());	// Using Boost - ifstream need path in form of c_str, not only string

	ifstream enthalpyCFDfile(dirCFDIn.c_str());	// Using string - ifstream need path in form of c_str, not only string

	enthalpyCFDfile.ignore(500,'\n');	// To skip the first line: 500 characters or until getting down

	if(enthalpyCFDfile)
	{
		//Open file succesfully
		string lineFile, temp;
	  	stringstream ss;		// Create a stringstream to store every single string

				
		//Count line & col
		while(getline(enthalpyCFDfile,lineFile))
		{
			maxRow ++;	// Count the line
			
			ss.clear();
			ss << lineFile;	// Get the lineFile to ss (reference)
			while(ss >> temp)	// while ss still exist (temp is a string, ss is an adress)
			{
				maxCol++;	// Count every single string then divide by maxRow to obtain maxCol
			}
	
		}
		maxCol = maxCol/maxRow;		// Number of column  = Number of string / Number of row
//		cout << " row = " << maxRow << endl;
//		cout << " col = " << maxCol << endl;


		// Read the file again from the begining
		enthalpyCFDfile.clear();			//Clear all error flag
		enthalpyCFDfile.seekg(0, enthalpyCFDfile.beg);	//and move the cursor to the begining to read the file again
		enthalpyCFDfile.ignore(500,'\n');		//To skip the first line, cursor moves to second line
		double data[maxRow][maxCol];		// Array to store the data as double
                for(int row = 0; row < maxRow; row++)
               	{
                      for(int col = 0; col < maxCol; col++)
                      {
                          enthalpyCFDfile >> data[row][col];	// Store data
                        //  cout << "i = " << row << " j = " << col << " data = " << data[row][col] <<  endl;
			

		      }
		}

		for(int col = 0; col < maxCol; col++)
		{
			if(col==0)
			{

				for(int row = 0; row < maxRow; row++)
				{
					t_CFD.push_back(data[row][col]); // Get the first column into vector
				//	cout << "t_CFD = " << t_CFD[row] << endl;
				}
			}
			else
			{
				 for(int row = 0; row < maxRow; row++)
                                {
                                        mean_hCFD.push_back(data[row][col]);	// Get the second column into vector
				//	cout << "mean_hCFD = " << mean_hCFD[row] << endl;
                                }
			}
		}

	}
	else
	{
		cout << "ERREUR: Unable to open  enthalpyCFDfile" << endl;
		cout << "Please make sure that the folder CFD_results was created and h_tau.txt exits" << endl;
	}

				
	enthalpyCFDfile.close();
    } // flagCFD bracket
    else // flagCFD = false
    {
	if(rank==0 && i==0) cout << "Without CFD correction" << endl;
    }


/* END Read CFD */
/* ============ */




	double checkMean_Hm; //Huu-Tri NGUYEN 13 Nov 2019
      	//Mean values
      for (int k=0; k<nsp; k++)
         Mean_Ym[k] = 0.0;
      Mean_Hm = 0.0;
      Mean_Tm = 0.0;
      Total_gas_mass = 0.0;
      for (int p=ndil; p<nTot; p++)
      {
         Total_gas_mass += listParticles[p]->m_P_gas_liquid; // m_P_gas_liquid for 1p = 1
         for (int k=0; k<nsp; k++)
            Mean_Ym[k] += (listParticles[p]->m_P_gas_liquid*listParticles[p]->m_Yk_gas[k]);
         Mean_Hm += (listParticles[p]->m_P_gas_liquid*listParticles[p]->m_H_gas);
         Mean_Tm += (listParticles[p]->m_P_gas_liquid*listParticles[p]->m_T_gas);
      }
      for (int k=0; k<nsp; k++)
       Mean_Ym[k] /= Total_gas_mass;
      Mean_Hm /= Total_gas_mass;
      Mean_Tm /= Total_gas_mass;
	
	// Huu-Tri commented - 2020.03.06
	if(rank==0) 
	{
		cout << endl;		
		cout << " Mean_Hm at ite " << i << " = " << Mean_Hm <<endl;
		cout << " Mean_Tm at ite " << i << " = " << Mean_Tm << endl;
	}
	
	/* Save the mean data Stochastic to output file - Huu-Tri NGUYEN - 19 Nov 2019 */
	bool activateMeanData = false;
	if(activateMeanData)
	{
		if(rank==0)	//Parallel MPI stuff to prevent print multiple lines, rank = 0 <=> First processor
		{
			if(i==0) cout << "*Mean data saving is activated" << endl;
			// Check if meanDataStochastic.txt file exists at the first step
			// If yes, clear the file content
			if(file_exists("outputs/mean_dataParticle.dat") && i==0)
			{
				cout << " -------------- Warrning --------------" << endl;
				cout << " outputs/mean_dataParticle.dat exists. Clearing file ... " << endl;	
			
				ofstream mean_dataParticle_clear("outputs/mean_dataParticle.dat", ios::out | ios::trunc);	//open file in trunc mode to clear the content
				mean_dataParticle_clear.close(); //close the file
	
			}		


	
			ofstream mean_dataParticle("outputs/mean_dataParticle.dat",ios::app); //ios::app = append at the end of the file
			if(mean_dataParticle)
			{
				if(i==0)	//First step: Need to write the headline
				{
					// First line
					mean_dataParticle << "#Time(s)	"; 
					for (int k=0; k<nsp; k++)
					{
						mean_dataParticle << "Mean_Ym_" << mixture->speciesName(k) << "	";
					}
					mean_dataParticle << "Mean_Hm	";
					mean_dataParticle << "Mean_Tm(K)" << endl;

					// Data from ORCh
					mean_dataParticle << i*delta_t << "	";
	
					for (int k=0; k<nsp; k++)
					{
						mean_dataParticle << Mean_Ym[k] << "	";
					}
					mean_dataParticle << Mean_Hm << "	";
					mean_dataParticle << Mean_Tm << "	" << endl;
				}	
				else
				{		
					// Data from ORCh
					mean_dataParticle << i*delta_t << "	";
					for (int k=0; k<nsp; k++)
					{
						mean_dataParticle << Mean_Ym[k] << "	";
					}
					mean_dataParticle << Mean_Hm << "	";
					mean_dataParticle << Mean_Tm << "	" << endl;
				}
			}
			else
			{	
				cout << "ERROR: Impossible to write mean_dataParticle.dat" << endl;
				cout << "Please check computeMultipleInlet.cpp" << endl;
			}
	
			mean_dataParticle.close();
		} // End if(rank==0)
	} // End if(activateMeanData)
	/* END Save - Huu-Tri NGUYEN - 19 Nov 2019 */


  
	/* ==================== CORRECTION  DETERMINISTIC EQUATIONS ================== */
	/*CFD Enthalpy Correction - Huu-Tri NGUYEN - 30 Aout 2019*/
	/**Commented 13 Nov 2019 to implement the heat loss to STOCHASTIC EQUATIONS**/

    if (flagCFDDeter)	// flagCFD == true, use CFD result correction 
     {
//if(i<7) // Use heat loss Deterministic in  5 time steps to stabilize the first jump 
//{


	/* === Interpolate the data to match with time step of ORCh - Huu-Tri NGUYEN - 5 Sept. 2019 === */
        vector<double> timeORCh;           // Store all the time of ORCh = Iteration * delta_t
        double  mean_hCFDinterpo[nbLines]; // Store the mean CFD enthalpy interpolated - Size = nb of iterations
	double Mean_Hm_ini; // Huu-Tri NGUYEN - 14.01.2020 - Interpolation

	if(i==0)	// Huu-Tri NGUYEN - 14.01.2020
	{ 
		Mean_Hm_ini = Mean_Hm;	// Mean_Hm of iteration t, we are now t+1 because after Stochastic closure
		if(rank ==0) cout << " Mean ini = " << Mean_Hm_ini << endl;
	}


	mean_hCFDinterpo[0] = Mean_Hm_ini;  // No variation between CFD and ORCh >> mean_hCFDinterpo - Mean_Hm = 0

	for(int step = 0; step < nbLines; step++)	
	{
		timeORCh.push_back(step*delta_t); // size of timeORCh vector equals to nbLines (nb of iterations)
	}
	

	// Find Top, Bot position
	for(int step = 1; step < nbLines; step++)	// Start from second value (step = 1), mean_hCFDinterpo[0] = Mean_Hm
        {
		for(int row = 0; row < maxRow-1; row++)	// End before the last number
		{
			if(timeORCh[step] < t_CFD[0]) // Verify if out of range
			{
                       //         mean_hCFDinterpo[step] = Mean_Hm; // Smaller than the first CFD result
									// assume no loss                                        
				double tTop = 0.0;	// tTop < timeORCh < tBot
				double tBot = t_CFD[0];
				double hTop = mean_hCFDinterpo[0];
				double hBot = mean_hCFD[0];


				mean_hCFDinterpo[step] = hTop + ((hBot - hTop)/(tBot-tTop)*(timeORCh[step]-tTop));			


			}
			else if(timeORCh[step] > t_CFD[maxRow-1]) //  Bigger than the last CFD result, take the last value hm
			{

				mean_hCFDinterpo[step] = mean_hCFDinterpo[step-1];
			} 
			else	// In the range of interpolation
			{
				if(timeORCh[step] > t_CFD[row+1])	// Find the cursor position
				{	
					// Do not thing >> Move to next row
				}
				else if(timeORCh[step] == t_CFD[row])
				{	
					mean_hCFDinterpo[step] = mean_hCFD[row];
				}	
				else if(timeORCh[step] > t_CFD[row] && timeORCh[step] < t_CFD[row+1]) // Interpolate
				{
					double tTop = t_CFD[row];	// tTop < timeORCh < tBot
					double tBot = t_CFD[row+1];
					double hTop = mean_hCFD[row];
					double hBot = mean_hCFD[row+1];
					mean_hCFDinterpo[step] = hTop + ((hBot - hTop)/(tBot-tTop)*(timeORCh[step]-tTop));
				}		
			}
			 	
		}
	} 
	

	// Save the interpolated result to output 1
	ofstream mean_hCFD_interpolated(dirCFDOut.c_str());	// Using string - ifstream need path in form of c_str, not only string


	if(mean_hCFD_interpolated)
	{
		mean_hCFD_interpolated << "#Time(s)	"; 
		mean_hCFD_interpolated << "h_mean_interpolated(J/kg)" << endl;

		for(int step = 0; step < nbLines; step++)
		{
			mean_hCFD_interpolated << timeORCh[step] << "	";
			mean_hCFD_interpolated << mean_hCFDinterpo[step] << endl;
		}
	}
	else
	{
		cout << "ERROR: Unable to write h_tau_Interpolated.txt" << endl;
		cout << "Please check computeMultipleInlet.cpp" << endl;
	}

	mean_hCFD_interpolated.close();


	/* =================================== END Interpolate =================================== */
	

	if(timeORCh[i] == i*delta_t)	// i =nbIteration, timeORCh = i*delta_t >> Check if interpolation is OK
	{	
		
		Mean_Hm = Mean_Hm + (mean_hCFDinterpo[i] - Mean_Hm);	// Add the enthalpy variation
		if(rank ==0)
		{		
			cout << "Correction Deterministic at " << i*delta_t << "s || " << endl;
			cout << "Mean correction = " << Mean_Hm << endl;
		}		
		 cout << "Mean Hm after = " << Mean_Hm << endl;
	}
	else
	{
		cout << "Something is wrong with interpolation - check computeMultipleInlet.cpp" <<  endl;
	}	
// } //if(i<5)  bracket

   } // flagCFD bracket

	/* ================== END CORRECTION ================= */

     // store << t << "  ";
     // for (int k=0; k<nsp; k++)
     //    store << Mean_Ym[k] << "  ";
     // store << Mean_Tm << "  ";
     // store << endl;



     //richesse
     double phi[nTot];
	

  // =========== SCATTERPLOT ===========
     // In this section, there are 2 formulas in order to calculate Mixture fraction Z for each particle (Z_gas and Zn2)
     // Z_gas is based on the Bilger formula (atoms C,H,O - default in ORCh) but only for 1 fuel
     // Zn2 is based on N2 evolution (No reaction of NOx and  - by Huu-Tri NGUYEN - 07.01.2020
     // Zn2_Deter for the mixing Deterministic will be calculated after the correction (Find Stochastic correction)
     bool activateScatter = false;	//flag to activate Scatterplot (write data_particles.dat)
     if(activateScatter)
     {
	if(rank==0)	//Parallel MPI stuff to prevent print multiple lines, rank = 0 <=> First processor
	{
		if(i==0) cout << "*Scatterplot is activated" << endl;
		// Check if meanSpeciesProdRate.txt file exists at the first step
		// If yes, clear the file content
		if(file_exists("outputs/scatterplot_data.dat") && i==0)
		{
		
			cout << " -------------- Warrning --------------" << endl;
			cout << " outputs/scatterplot_data.dat exists. Clearing file ... " << endl;	
			
			ofstream store_particles_clear("outputs/scatterplot_data.dat", ios::out | ios::trunc);	//open file in trunc mode to clear the content
			store_particles_clear.close(); //close the file
	
		}		


		ofstream store_particles("outputs/scatterplot_data.dat",ios::app); //ios::app = append at the end of the file
		if(store_particles)
		{
			if(i==0)	// write the first line
			{  
				store_particles << "#1:time  2:Particle_number  3:Zfraction  4:T(K) 5:Y_CO2 6:Y_O2 7:particle type  8:ratio  9:Zst	10:Zc	11:Zo	12:Zh 13:Zn2	14:hp" << endl;
			}	

      			for (int z=ndil;z<nTot; z++)
      			{
      				double Yf, Yo, Yf_0, Yo2_0, s, Yco2, Ych4, Yco, Yh2o;
				
      				Yf_0 = 1; 
      				Yo2_0 = 0.232917;

				// Huu-Tri NGUYEN - 07.01.2020 - Calculate Z (mixture fraction) by N2
				double Yn2, Yn2_0, Yn2_f,Zn2;
				Zn2 = 0;
				Yn2_0 = 0.766990291; 	// UMONS case - Inlet air preheated 
				Yn2_f = 0.396761134;	// UMONS case - Inlet fuel
      
      				for (int k=0; k<nsp; k++)
      				{
         				if (mixture->speciesName(k) == "NC10H22")
            					Yf = listParticles[z]->m_Yk_gas[k]; 
         				if (mixture->speciesName(k) == "O2")
            					Yo = listParticles[z]->m_Yk_gas[k]; 
         				if (mixture->speciesName(k) == "CO2")
            					Yco2 = listParticles[z]->m_Yk_gas[k]; 
         				if (mixture->speciesName(k) == "CH4")
            					Ych4 = listParticles[z]->m_Yk_gas[k]; 
         				if (mixture->speciesName(k) == "H2O")
            					Yh2o = listParticles[z]->m_Yk_gas[k]; 
         				if (mixture->speciesName(k) == "CO")
            					Yco = listParticles[z]->m_Yk_gas[k]; 
         				if (mixture->speciesName(k) == "N2")		//Huu-Tri NGUYEN - 07.01.2020
            					Yn2 = listParticles[z]->m_Yk_gas[k];
				}      
       
				Zn2 = (Yn2-Yn2_0)/(Yn2_f-Yn2_0); //Huu-Tri NGUYEN - 07.01.2020: Mixture fraction based on N2
									// Zst_n2 is calculated in Excel file of case conditions

      				//mass stoichiometric coefficient
      				double W_NC10H22 = 142; // kg/kmol
      				double W_O2 = 32; // kg/kmol
      				double W_CH4 = 16;
      				double W_CO = 28;
      				double W_CO2 = 44;
      				double W_H2O = 18;
      				double W_C = 12;
      				double W_O = 16;
      				double W_H = 1;  
				double W_H2 = 2;	// Huu-Tri Nguyen - 16.01.2020
				double W_N2 = 28;  	// Huu-Tri Nguyen - 16.01.2020
 
      				int m = 6;	//22
      				int n = 2;     //10
      				int nuO = n + m/4;

      				double Zc = 0;
      				double Zo = 0;
      				double Zh = 0;

				double Y[nsp];
				for (int f=0;f<nsp;f++)
   					Y[f] = listParticles[z]->m_Yk_gas[f];


      				//Species mixture fraction - Bilger formula 
				// Zc = (No_atomsC_inSpecie * Weight_of_C * Y_species) / Weight_of_specie)
      				for (int k=0;k<nsp;k++)
      				{
         				if (mixture->species(k)->composition["C"] != 0)
            					Zc += mixture->species(k)->composition["C"]*W_C*Y[k]/( mixture->species(k)->composition["C"]*W_C +  mixture->species(k)->composition["H"]*W_H + mixture->species(k)->composition["O"]*W_O);
         				if (mixture->species(k)->composition["O"] != 0)
         					Zo += mixture->species(k)->composition["O"]*W_O*Y[k]/( mixture->species(k)->composition["C"]*W_C +  mixture->species(k)->composition["H"]*W_H +  mixture->species(k)->composition["O"]*W_O);
         				if (mixture->species(k)->composition["H"] != 0)
         					Zh += mixture->species(k)->composition["H"]*W_H*Y[k]/( mixture->species(k)->composition["C"]*W_C +  mixture->species(k)->composition["H"]*W_H +  mixture->species(k)->composition["O"]*W_O);
      				}				


     				//Species mixture fraction 
     				double Zc0 = n*W_C*Yf_0/(n*W_C + m*W_H);
     				double Zh0 = m*W_H*Yf_0/(n*W_C + m*W_H);
     				double Zo0 = 2*Yo2_0*W_O/(W_O2);



      				//Calcul de la fraction de mélange locale Z
      				listParticles[z]->m_Z_gas = (Zc/(n*W_C) + Zh/(m*W_H) + 2*(Yo2_0 - Zo)/(nuO*W_O2))/(Zc0/(n*W_C) + Zh0/(m*W_H)+ 2*Yo2_0/(nuO*W_O2));


				//Huu-Tri NGUYEN - check 
				//if(rank == 0)
				//	cout << " Particle " << z << " !!! delta = " << listParticles[z]->m_Z_gas-Zn2 << endl;

      				//Calcul de Zst
      				double Zst = (2*Yo2_0/(nuO*W_O2))/(Zc0/(n*W_C) + Zh0/(m*W_H)+ 2*Yo2_0/(nuO*W_O2));

      				//Calcul de la richesse locale
      				phi[z] =   listParticles[z]->m_Z_gas * (1 - Zst) / (Zst * (1 - listParticles[z]->m_Z_gas));

      				double particleType;
      				if (z < nbParticles[0])
         				particleType = 0;
      				else if (z > nbParticles[0]-1 && z < nbParticles[0] + nbParticles[1])
         				particleType = 1;
      				else if (z >= nbParticles[0] + nbParticles[1])
         				particleType = 2;
     

				// Store the result in data_particles.dat  
				store_particles << t << "  ";
      				store_particles << z+1 << "  ";
      				store_particles << listParticles[z]->m_Z_gas << "  ";
      				store_particles << listParticles[z]->m_T_gas  << "  ";
      				store_particles << Yf  << "  ";
      				store_particles << Yo  << "  ";
     				store_particles << particleType  << "  ";
     				store_particles << phi[z] << "  ";
	     			store_particles << Zst << "  ";    
	     			store_particles << Zc << "  ";    
     				store_particles << Zo << "  ";    
     				store_particles << Zh << "  ";  
  				store_particles << Zn2 << "  ";		//Huu-Tri NGUYEN - 07.01.2020
  				store_particles << listParticles[z]->m_H_gas << "  ";//Huu-Tri NGUYEN - 15.01.2020 - Enthalpy de particule
     				store_particles << endl; 
			
		

      			} //end of scatterplot part
		}
    		store_particles.close(); //store_particles for scatterplot

	
	// Huu-Tri Nguyen - Print enthalpy evolution of each inlet - 15.01.2020

	// Calculate Zn2 Deterministic
		double Yn2_0, Yn2_f;
		double Yn2_Deter0 = 0, Yn2_Deter1 = 0, Yn2_Deter2 = 0, Yn2_DeterMix = 0;;
		double Zn2_Deter0 = 0, Zn2_Deter1 = 0, Zn2_Deter2 = 0, Zn2_DeterMix = 0;
		Yn2_0 = 0.766990291; 	// UMONS case - Inlet air preheated 
		Yn2_f = 0.396761134;	// UMONS case - Inlet fuel
	    
		for (int k=0; k<nsp; k++)
		{	
			
			if (mixture->speciesName(k) == "N2")		//Huu-Tri NGUYEN - 07.01.2020
			{
				Yn2_DeterMix = Ym[k];
       		     		Yn2_Deter0 = Ym_Trajectories_store[0][i][k];
				Yn2_Deter1 = Ym_Trajectories_store[1][i][k];
			//	Yn2_Deter2 = Ym_Trajectories_store[2][i][k];	// Huu-Tri Commented - 2 inlets - 20.01.2020

			}
		}
			   
		Zn2_DeterMix = (Yn2_DeterMix-Yn2_0)/(Yn2_f-Yn2_0);
		Zn2_Deter0 = (Yn2_Deter0-Yn2_0)/(Yn2_f-Yn2_0);
		Zn2_Deter1 = (Yn2_Deter1-Yn2_0)/(Yn2_f-Yn2_0);
		Zn2_Deter2 = (Yn2_Deter2-Yn2_0)/(Yn2_f-Yn2_0);






		// Calculate enthalpy mix of inlets - UMONS case - Huu-Tri Nguyen 15.01.2020
		double h_gas = Hm_Trajectories[0];
		double h_air = Hm_Trajectories[1];
		double h_mixZ_Deter;
			h_mixZ_Deter = h_gas*Zn2_DeterMix + h_air*(1-Zn2_DeterMix);	// In the case of 2 inlets adidabatic
		double h_burntGas = Hm_Trajectories[2];


	
		if(file_exists("outputs/scatterplot_dataHDeter.dat") && i==0)
		{
			cout << " -------------- Warrning --------------" << endl;
			cout << " outputs/scatterplot_dataHDeter.dat exists. Clearing file ... " << endl;	
				
			ofstream dataEnthalpy_clear("outputs/scatterplot_dataHDeter.dat", ios::out | ios::trunc);	//open file in trunc mode to clear the content
			dataEnthalpy_clear.close(); //close the file
		} //end if(file_exists)		

		ofstream dataEnthalpy("outputs/scatterplot_dataHDeter.dat",ios::app); //ios::app = append at the end of the file
		if(dataEnthalpy)
		{
			if(i==0)	// write the first line
			{  
				dataEnthalpy << "#1:time	2:h_gas	3:h_air	4:h_GB	5:h_mixZ_Deter	6:Zn2_Mix	7:Zn2_gas	8:Zn2_air	9:Zn2_GB" << endl;
				dataEnthalpy << i << "	";
				dataEnthalpy << h_gas << "	";
				dataEnthalpy << h_air << "	";
				dataEnthalpy << h_burntGas << "	";
				dataEnthalpy << h_mixZ_Deter << "	";
				dataEnthalpy << Zn2_DeterMix << "	";
				dataEnthalpy << Zn2_Deter0 << "	";
				dataEnthalpy << Zn2_Deter1 << "	";
				dataEnthalpy << Zn2_Deter2 << "	";
				dataEnthalpy << endl;
			
			}
			else
			{
				dataEnthalpy << i << "	";
				dataEnthalpy << h_gas << "	";
				dataEnthalpy << h_air << "	";
				dataEnthalpy << h_burntGas << "	";
				dataEnthalpy << h_mixZ_Deter << "	";
				dataEnthalpy << Zn2_DeterMix << "	";
				dataEnthalpy << Zn2_Deter0 << "	";
				dataEnthalpy << Zn2_Deter1 << "	";
				dataEnthalpy << Zn2_Deter2 << "	";
				dataEnthalpy << endl;
			}
	
		} //end if(dataEnthalpy)	
		dataEnthalpy.close();
	} // end if(rank==0)	
    } //end if(activeScatter)	

     // =========== END SCATTERPLOT ===========	



     //store.close(); // for mean values >> already made a new file mean_dataParticle.txt- Huu-Tri Nguyen - 10.12.2019


     // =========== FULL DATA PARTICLES - Huu-Tri Nguyen 19.12.2019===========
     // Store Temperature & Mass fraction for each particle at each time step - Huu-Tri Nguyen 10.12.2019
     bool activateFullData = true;	//flag to activate Full data (write full_dataParticle.dat)
     if(activateFullData)
     {
	if(rank==0)	//Parallel MPI stuff to prevent print multiple lines, rank = 0 <=> First processor
	{
		if(i==0) 
			cout << "*Full data saving is activated" << endl;
		// Check if meanSpeciesProdRate.txt file exists at the first step
		// If yes, clear the file content
		if(file_exists("outputs/full_dataParticle.dat") && i==0)
		{
		
			cout << " -------------- Warrning --------------" << endl;
			cout << " outputs/full_dataParticle.dat exists. Clearing file ... " << endl;	
			
			ofstream full_dataParticle_clear("outputs/full_dataParticle.dat", ios::out | ios::trunc);	//open file in trunc mode to clear the content
			full_dataParticle_clear.close(); //close the file
	
		}	
		
		ofstream full_dataParticle("outputs/full_dataParticle.dat",ios::app); //ios::app = append at the end of the file
		if(full_dataParticle)
		{
			if(i==0)
			{
   				ofstream full_dataParticle ("outputs/full_dataParticle.dat"); 
   				full_dataParticle << "#Time	";
				full_dataParticle << "Particle_number	";
   				for (int k=0; k<nsp; k++)
      					full_dataParticle << mixture->speciesName(k) << "	";
				full_dataParticle << "Hm	";	//HT@2020.08.22 : Need to remove
				full_dataParticle << "Temperature	" << endl;

				for(int p=0; p<nTot; p++)
				{
     					full_dataParticle << t << "	";
     					full_dataParticle << p << "	";	// Particle starts from 1
     					for (int k=0; k<nsp; k++)
        					full_dataParticle << listParticles[p]->m_Yk_gas[k] << "	";
					full_dataParticle << listParticles[p]->m_H_gas << "	";	//HT2020.08.22 : Need to remove
					full_dataParticle << listParticles[p]->m_T_gas << "	" << endl;
				}
			}
			else
			{
				for(int p=0; p<nTot; p++)
				{
     					full_dataParticle << t << "	";
     					full_dataParticle << p << "	";	// Particle starts from 1
     					for (int k=0; k<nsp; k++)
        					full_dataParticle << listParticles[p]->m_Yk_gas[k] << "	";
                                        full_dataParticle << listParticles[p]->m_H_gas << "	";      //HT2020.08.22 : Need to remove
					full_dataParticle << listParticles[p]->m_T_gas << "	" << endl;
				}	
			}
		}

		full_dataParticle.close(); //close the file
	} //END if(rank==0)

     } //END if(activateFullData)
     // =========== END FULL DATA PARTICLES ===========

      // ====== LAGRANGIAN TRAJECTORIES - DETERMINISTIC ====== //	
      for (int n=0; n<nbInlets-1; n++)
      {
         for (int k=0; k<nsp; k++)
            Ym[k] = (Ym_Trajectories_store[n][i][k]-Mean_Ym[k])*exp(-(delta_t/(2*tau_t)))+Mean_Ym[k];

         
	
         Hm = (Hm_Trajectories[n]-Mean_Hm)*exp(-(delta_t/(2*tau_t)))+Mean_Hm;

         if (step == "DRGEP_Species" || step == "DRGEP_Reactions")
         {
            Next_Time_Step_with_drgep(Targets, Pressure, Ym, Hm, Tm, delta_t, R_AD_Trajectories[n], max_j_on_Target, step);
         }
         else if (step == "computeQSSCriteria")
         {
            Next_Time_Step(Pressure, Ym, Hm, Tm, delta_t, Production_Trajectories_ref, Consumption_Trajectories_ref, n, i);
         }
         else
         {

            Next_Time_Step(Pressure, Ym, Hm, Tm, delta_t);
             
             // Huu-Tri TEST - 20210206
            if(rank==0)
            {
                cout << " **** Advance Deterministic inlet " << n << " at " << i << " iterations" << endl;
            }
         }

	// Huu-Tri - After the reactor
         for (int k=0; k<nsp; k++)
            Ym_Trajectories_store[n][i+1][k] = (Ym[k]-Mean_Ym[k])*exp(-(delta_t/(2*tau_t)))+Mean_Ym[k];


         Hm_Trajectories[n] = (Hm-Mean_Hm)*exp(-(delta_t/(2*tau_t)))+Mean_Hm;
         T_Trajectories_store[n][i+1] = Tm;

      } // End "for" each inlet
      // ====== END LAGRANGIAN TRAJECTORIES - DETERMINISTIC ====== //



	//  Huu-Tri NGUYEN - Calculate source term Deterministic - 5 Dec 2019 
	//  1st step: Create the 2D array
	bool calculWDeter = false;
	if(calculWDeter)
	{
		double *wDeter_T = new double[nbInlets-1];
		double **wDeter_species  = new double*[nbInlets-1];	// 2D array ** with first dimension (number of inlets)
		for (int n=0; n<nbInlets-1; n++)	// Create second dimension 
		{
			wDeter_species[n] = new double[nsp];	// nsp = number of species has declared above
		}

		// 2nd step: Calculate source term
		for (int n=0; n<nbInlets-1; n++)
      		{	wDeter_T[n] = T_Trajectories_store[n][i+1] - T_Trajectories_store[n][i];
			wDeter_T[n] /= delta_t;

        		for (int k=0; k<nsp; k++)
			{
				wDeter_species[n][k] = (Ym_Trajectories_store[n][i+1][k] - Ym_Trajectories_store[n][i][k]);
				wDeter_species[n][k] /= delta_t;	
			}
		}
	
		// 3rd step: Save to file
		SaveToFile_wDeter(wDeter_T, wDeter_species, nsp, nbInlets, i, delta_t, rank, mixture);

	
		// Free w_species memory (array pointer should be deleted after use)		 
		delete wDeter_T;
		for (int n=0; n<nbInlets-1; n++)	// Create second dimension 
		{
			delete[] wDeter_species[n];
		}
		delete[] wDeter_species;
	} //End if(calculWDeter)	
	//  End Calculate source term deterministic - 5 Dec 2019
	

      //---Particles mixing---
      //
      //Add Nmix variation if dilution
     int Nmix2;
      if (i < nbLines/2)
	Nmix2 = delta_t*(nTot-AirPart+1)/tau_t;
      else 
         {
        
         float c = i;
         float d = nbLines;
         double b = abs(AirPart*(1 - ( 2*( c - d/2 )/d)));
         int a = b;
         Nmix2 = delta_t*(nTot - a)/tau_t;
         }

       /* ======== MIXING MODELS ======== */
       /* 2 options: Curl model or EMST model */
	bool activateCurl = false;	//true = Curl; false = EMST

	//Huu-Tri@20200922 : Save Y, T of particle before EMST to check  if EMST change particle - ANN
	double saveEMSTbefore[nTot][nsp];	
         for (int p=0; p<nTot; p++) //Copy mass fraction of each particle into saveEMSTbefore array
         {
            	for (int k=0; k<nsp; k++)
            	{
			saveEMSTbefore[p][k] = listParticles[p]->m_Yk_gas[k];	// Already checked cout 

            	}
         } 




	 /* CURL Mixing closure phi_p1 = phi_p2 = (phi_p1 + phi_p2)/2 */
	if(activateCurl)
	{
		if(rank ==0) cout << "Curl mixing model" << endl;
     	
      		for (int p=0; p<Nmix; p++)
      		{
          		double run_mixing = true;

          		if (run_mixing)
         		{
 	    		// Huu-Tri: Only gas case,  gas_mass_p1 = gas_mass_p2 = 0.01 = Particle_flowRate
            			gas_mass_p1 = listParticles[Particle_1[i][p]]->m_P_gas_liquid*Particle_flowRate;
             			gas_mass_p2 = listParticles[Particle_2[i][p]]->m_P_gas_liquid*Particle_flowRate; 

             			if (gas_mass_p1+gas_mass_p2 > 0.0)
           			{
                			for (int k=0; k<nsp; k++)
                			{
                   				listParticles[Particle_1[i][p]]->m_Yk_gas[k] = F*(gas_mass_p1*listParticles[Particle_1[i][p]]->m_Yk_gas[k]+gas_mass_p2*listParticles[Particle_2[i][p]]->m_Yk_gas[k])/(gas_mass_p1+gas_mass_p2);
                   				listParticles[Particle_2[i][p]]->m_Yk_gas[k] = listParticles[Particle_1[i][p]]->m_Yk_gas[k];
                			}

      		 				listParticles[Particle_1[i][p]]->m_H_gas = (gas_mass_p1*listParticles[Particle_1[i][p]]->m_H_gas+gas_mass_p2*listParticles[Particle_2[i][p]]->m_H_gas)/(gas_mass_p1+gas_mass_p2);//Original

               
//  Huu-Tri NGUYEN - Add a modified enthaply for heat loss
//                listParticles[Particle_1[i][p]]->m_H_gas = (gas_mass_p1*listParticles[Particle_1[i][p]]->m_H_gas+gas_mass_p2*listParticles[Particle_2[i][p]]->m_H_gas)/(gas_mass_p1+gas_mass_p2) + varEnthalpyCFD;// Huu-Tri Stochastic heat loss 14 Nov 2019 


						listParticles[Particle_2[i][p]]->m_H_gas = listParticles[Particle_1[i][p]]->m_H_gas;
	

             			}
         		}
       		}
	}

      /* END Commented Curl model to add EMST model - Huu-Tri Nguyen -10.12.2019 */
	else // EMST model
	{
		if(rank ==0) 
		{
			cout << "EMST mixing model" << endl;
		}
		/* Add EMST model - Huu-Tri Nguyen -10.12.2019 */
     		//EMST mixing model -- by Kaidi@2019.12
      		double run_mixing = true;

     		if (run_mixing==true)
      		{
         		mode_emst = 2;	// 2 -mixing is performed and the state variables are incremented.

         		i_emst = 0;
         		for (int ic=0; ic<nc_emst; ic++)
         		{
            			for (int ip=0; ip<np_emst; ip++)
            			{
               				if (ic<nsp)  
						f_emst[i_emst] = listParticles[ip]->m_Yk_gas[ic];
               				else         
						f_emst[i_emst] = listParticles[ip]->m_H_gas;
               			i_emst++;
            			}
         		}
    
		// Mixing with EMST model
		bool byPassEMST = false; // By pass EMST to use only ANN - HuuTri@20200919
		if (byPassEMST == false)
		{
			if(rank==0 & i==0) {
                                cout << " Use EMST - Not by pass " << endl;
                        }
			emst_(&mode_emst,&np_emst,&nc_emst,f_emst,state_emst,wt_emst,&omdt_emst,fscale_emst,cvars_emst,&info_emst);
			
			int EMST_mixedParticles =0;			
			for (int ip=0; ip<np_emst; ip++) //HuuTri@20201029: Print number of mixed particles
  			{
				if(state_emst[ip] >0)
				{
      					EMST_mixedParticles++; // Particle mixed state_emst>0, non-mixed state_emst<0
				}
      				// wt_emst[ip];
   			}
			if(rank==0)
			{ 
				cout << "Mixed particles = " << EMST_mixedParticles << endl;
			}        	
		}
		else
		{
			if(rank==0 & i==0) {
				cout << " By pass EMST for ANN !!!! " << endl;
			}
		}
			if (info_emst != 0)
         		{
            			cout << "emst failed" << endl;
            			getchar();
         		}

         		i_emst = 0;


         		for (int ic=0; ic<nc_emst; ic++)
         		{
            			for (int ip=0; ip<np_emst; ip++)
            			{
               				if (ic<nsp)
               				{
                  				listParticles[ip]->m_Yk_gas[ic] = f_emst[i_emst];
                  				//listParticles[ip]->m_Yk_gas[ic] = listParticles[ip]->m_Yk_gas[ic] * (((rand() / double(RAND_MAX))*2.0-1.0)*1.0/100.0 + 1.0);
               				}
               				else
               				{
                  				listParticles[ip]->m_H_gas = f_emst[i_emst];
                  				//T_o = 353;
                  				//if (listParticles[ip]->m_T_gas > T_o) listParticles[ip]->m_H_gas = listParticles[ip]->m_H_gas * (((rand() / double(RAND_MAX))*2.0-1.0)*1.0/100.0 + 1.0);
                  				//if (i > nbLines/3)
						//Heatloss Camille - Commented by Huu-Tri Nguyen -10.12.2019
				                  			
//HT						if (true)
//HT                  				{
//HT                    				T_o = 353; 
//HT                     				if (listParticles[ip]->m_T_gas > T_o) 
//HT							listParticles[ip]->m_H_gas = listParticles[ip]->m_H_gas - K*alpha_loss*(listParticles[ip]->m_T_gas - T_o)*delta_t;
//HT              				}


						// Heatloss Camille - Corrected by Huu-Tri Nguyen - 20220107
    						// HuuTri@20220107: flagHeatLossStochas == true, use a sink/source term to impose the heat loss (alpha_loss) for Stochastic closure
						// Sink term  = alpha_loss*(T_gas - T_wall)
						// Heat term  = beta_heat*(T_gas - T_wall)
						// The value of alpha_loss, beta_heat also depends on the time step delta_t
						// To calculate enthalpy: It should be multiplied by delta_t: 
						// H_gasAfterLoss = H_gas -  sinkTerm*delta_t =  H_ini - alpha_loss*(T_gas - T_wall)*delta_t
						// H_gasAfterHeated = H_gas -  sourceTerm*delta_t =  H_ini - beta_heat*(T_gas - T_wall)*delta_t
						if (flagHeatLossStochas)
						{
							// HuuTri@20220107: alpha_loss and T_wall should be declared above (find "flagHeatLossStochas" keyword)
							// HuuTri@20220107: If Tgas > Twall, the heat will be lost >> alpha_loss
							if (listParticles[ip]->m_T_gas > T_wall)
							{
								listParticles[ip]->m_H_gas = listParticles[ip]->m_H_gas - alpha_loss*(listParticles[ip]->m_T_gas - T_wall)*delta_t;	
							}
							else // When Tgas < Twall, the gas is heated by walls
							{
								listParticles[ip]->m_H_gas = listParticles[ip]->m_H_gas - beta_heat*(listParticles[ip]->m_T_gas - T_wall)*delta_t;
							}


						} // end if (flagHeatLossStochas)
              				}
               				i_emst++;
            			}
         		}
    	
		}
	} // End else EMST model
      // End of EMST mixing model

	/* END Add EMST model - Huu-Tri Nguyen -10.12.2019 */

	



	// ===== SCATTERPLOT Particle 2 - AFTER the correction Enthalpy - Huu-Tri Nguyen - 16.01.2020 
    if(activateScatter)
    {
	if(rank==0)	//Parallel MPI stuff to prevent print multiple lines, rank = 0 <=> First processor
	{

		if(file_exists("outputs/scatterplot_dataAfterCorrection.dat") && i==0)
		{
		
			cout << " -------------- Warrning --------------" << endl;
			cout << " outputs/scatterplot_dataAfterCorrection.dat exists. Clearing file ... " << endl;	
			
			ofstream store_particlesAfter_clear("outputs/scatterplot_dataAfterCorrection.dat", ios::out | ios::trunc);	//open file in trunc mode to clear the content
			store_particlesAfter_clear.close(); //close the file

		}		

		ofstream store_particlesAfter("outputs/scatterplot_dataAfterCorrection.dat",ios::app); //ios::app = append at the end of the file
		if(store_particlesAfter)
		{
			if(i==0)	// write the first line
			{  
				store_particlesAfter << "#1:time  2:Particle_number  3:particle type  4:Zn2	5:hp" << endl;
			}	
		
      			for (int z=ndil;z<nTot; z++)
      			{

				// Huu-Tri NGUYEN - 07.01.2020 - Calculate Z (mixture fraction) by N2
				double Yn2, Yn2_0, Yn2_f,Zn2;
				Zn2 = 0;
				Yn2_0 = 0.766990291; 	// UMONS case - Inlet air preheated 
				Yn2_f = 0.396761134;	// UMONS case - Inlet fuel
      
      				for (int k=0; k<nsp; k++)
      				{
         			 
         				if (mixture->speciesName(k) == "N2")		//Huu-Tri NGUYEN - 07.01.2020
            						Yn2 = listParticles[z]->m_Yk_gas[k];
				}      
       
				Zn2 = (Yn2-Yn2_0)/(Yn2_f-Yn2_0); //Huu-Tri NGUYEN - 07.01.2020: Mixture fraction based on N2
									// Zst_n2 is calculated in Excel file of case conditions

      				

      				double particleType;
      				if (z < nbParticles[0])
         				particleType = 0;
      				else if (z > nbParticles[0]-1 && z < nbParticles[0] + nbParticles[1])
         				particleType = 1;
      				else if (z >= nbParticles[0] + nbParticles[1])
         				particleType = 2;
     

				// Store the result in data_particles.dat  
				store_particlesAfter << t << "  ";
      				store_particlesAfter << z+1 << "  ";
     				store_particlesAfter << particleType  << "  "; 
  				store_particlesAfter << Zn2 << "  ";		//Huu-Tri NGUYEN - 07.01.2020
  				store_particlesAfter << listParticles[z]->m_H_gas << "  ";//Huu-Tri NGUYEN - 15.01.2020 - Enthalpy de particule
     				store_particlesAfter << endl; 
			
		
      			} //end of scatterplot part
		}
    		store_particlesAfter.close(); //store_particlesAfter for scatterplot after the correction
		
	} //End(rank==0)
    }// End if(activateScatter)



     //Evaporation
      for (int p=ndil; p<nTot; p++)
      {
         if (listParticles[p]->m_P_gas_liquid < 1.0)
         {
            double Diameter;
            if ((t+dt)/tau_vj < 1)
               Diameter = Diameter_init*pow((1-((t+dt)/tau_vj)),0.5);
            else
               Diameter = 0.0;

            double old_Diameter = listParticles[p]->m_droplets_diameter;
            listParticles[p]->m_droplets_diameter = Diameter;

            double gas_mass = listParticles[p]->m_P_gas_liquid*Particle_flowRate;
            double liquid_mass = (1-listParticles[p]->m_P_gas_liquid)*Particle_flowRate;
            double ReleasedVapor = listParticles[p]->m_N_droplets*listParticles[p]->m_density_liquid*(PI/6)*(pow(old_Diameter,3)-pow(Diameter,3));
            listParticles[p]->m_P_gas_liquid = (gas_mass+ReleasedVapor)/Particle_flowRate;

            for (int k=0; k<nsp; k++)
            {
               double Yk_gas_init = listParticles[p]->m_Yk_gas[k];
               listParticles[p]->m_Yk_gas[k] = ((gas_mass*Yk_gas_init)+(ReleasedVapor*listParticles[p]->m_Yk_liquid[k]))/(gas_mass+ReleasedVapor);
            }

            double H_gas_init = listParticles[p]->m_H_gas;
            double ReleasedEnthalpy = 0.0;
            for (int k=0; k<nsp; k++)
            {
               ReleasedEnthalpy += ReleasedVapor*listParticles[p]->m_Yk_liquid[k]*(BoilingEnthalpy[k]-listParticles[p]->m_EvaporationLatentHeat);
            }
            listParticles[p]->m_H_gas = ((gas_mass*H_gas_init)+(ReleasedEnthalpy))/(gas_mass+ReleasedVapor);
         }


      }


    // Calculate right hand-side of Y(=mixing + source term) and T (=source term) - Huu-Tri NGUYEN - 2019.12.05
    bool activateSourceTermParticle = false;
	// Initialization
	double *Tm_gas_before = new double[nTot];	// 1D array temperature for each particle
	double **Yk_gas_before  = new double*[nTot];	// 2D array ** with first dimension (number of particles)
	for (int p=0; p<nTot; p++)	// Create second dimension 
	{
		Yk_gas_before[p] = new double[nsp];	// nsp = number of species has declared above
	}


	for (int p=ndil; p<nTot; p++)
      	{

		Tm_gas_before[p] = listParticles[p]->m_T_gas;		

 		for (int k=0; k<nsp; k++)
        	{
              	 	Yk_gas_before[p][k] = listParticles[p]->m_Yk_gas[k];
		}	
	}
   // End Calculate right hand-side of Y(=mixing + source term) and T (=source term) - Huu-Tri NGUYEN - 2019.12.05

     // =========== FULL DATA PARTICLES  After EMST - Huu-Tri Nguyen 17.09.2020===========
     // Store Temperature & Mass fraction for each particle at each time step BUT AFTER EMST- Huu-Tri Nguyen 17.09.2020
     // Move from 1205 to Here
     // ORCh : dY/dt = EMST then this one then dY/dt = wdot
     bool activateFullDataAftEMST = true;	//flag to activate Full data (write full_dataParticle.dat)
     if(activateFullDataAftEMST)
     {
	if(rank==0)	//Parallel MPI stuff to prevent print multiple lines, rank = 0 <=> First processor
	{
		cout << "Save after EMST" << endl;
		if(i==0) 
			cout << "*Full data AFTER EMST saving is activated" << endl;
		// Check if meanSpeciesProdRate.txt file exists at the first step
		// If yes, clear the file content
		if(file_exists("outputs/full_dataParticleAftEMST.dat") && i==0)
		{
		
			cout << " -------------- Warrning --------------" << endl;
			cout << " outputs/full_dataParticleAftEMST.dat exists. Clearing file ... " << endl;	
			
			ofstream full_dataParticle_clear("outputs/full_dataParticleAftEMST.dat", ios::out | ios::trunc);	//open file in trunc mode to clear the content
			full_dataParticle_clear.close(); //close the file
	
		}	
		
		ofstream full_dataParticle("outputs/full_dataParticleAftEMST.dat",ios::app); //ios::app = append at the end of the file
		if(full_dataParticle)
		{
			if(i==0)
			{
   				ofstream full_dataParticle ("outputs/full_dataParticleAftEMST.dat"); 
   				full_dataParticle << "#Time	";
				full_dataParticle << "Particle_number	";
   				for (int k=0; k<nsp; k++)
      					full_dataParticle << mixture->speciesName(k) << "	";
				full_dataParticle << "Hm	";	//HT@2020.08.22 : Need to remove
				full_dataParticle << "Temperature	" << endl;

				for(int p=0; p<nTot; p++)
				{
     					full_dataParticle << t << "	";
     					full_dataParticle << p << "	";	// Particle starts from 1
     					for (int k=0; k<nsp; k++)
        					full_dataParticle << listParticles[p]->m_Yk_gas[k] << "	";
					full_dataParticle << listParticles[p]->m_H_gas << "	";	//HT2020.08.22 : Need to remove
					full_dataParticle << listParticles[p]->m_T_gas << "	" << endl;
				}
			}
			else
			{
				for(int p=0; p<nTot; p++)
				{
     					full_dataParticle << t << "	";
     					full_dataParticle << p << "	";	// Particle starts from 1
     					for (int k=0; k<nsp; k++)
        					full_dataParticle << listParticles[p]->m_Yk_gas[k] << "	";
                                        full_dataParticle << listParticles[p]->m_H_gas << "	";      //HT2020.08.22 : Need to remove
					full_dataParticle << listParticles[p]->m_T_gas << "	" << endl;
				}	
			}
		}

		full_dataParticle.close(); //close the file
	} //END if(rank==0)

     } //END if(activateFullDataAftEMST)
     // =========== END FULL DATA PARTICLES AFTER EMST ===========


    // Huu-Tri@20200724 : Load model with a path to the .pb file. (Using cppflow)
    bool flagANN = false;


	// HuuTri@20220111: Commented
	//if(rank==0)	//HT@2020.08.20 check after test
        //{
        //        cout << "Particle 700 before - T= " << listParticles[700]->m_T_gas << " |Enthalpy = " <<  listParticles[700]->m_H_gas << endl;
        //}


	/* ==================================== */
	/* ============ STOCHASTIC ============ */
	/* ==================================== */
      if (step == "Optimisation")
         Reacting(listParticles, nsp, dt, Pressure);
       
      else
	{
	 // Use ANN to predict the source term, then calculate the next-iteration state of each particle - Huu-Tri@20200729
	 	if (flagANN==true)  // Huu-Tri@20200724
		{
			if (rank==0) cout << "Use ANN to advance Y and T" << endl;
		}
	 	else	// Use Cantera to integrate - obtain the next step
		{
			if(i==0)
			{
				if(rank==0)
				{	
					cout << "***Use Cantera ConstPressureReactor to advance" << " |flagANN  = " << flagANN << endl;
         			}
			}
			ReactingParallel(listParticles, nsp, dt, Pressure);
		}
	}	
	/* ==================================== */
	/* ==================================== */

	// =========== FULL DATA PARTICLES  After Reactor - Huu-Tri Nguyen 17.12.2020 ===========
     // Store Temperature & Mass fraction for each particle at each time step BUT AFTER REACTOR- Huu-Tri Nguyen 17.12.2020
     // This is the label of ANN Regression
     // ORCh : dY/dt = EMST then AftEMST then dY/dt = wdot then AftREACTOR
     bool activateFullDataAftREACTOR = true;	//flag to activate Full data (write full_dataParticle.dat)
     if(activateFullDataAftREACTOR)
     {
	if(rank==0)	//Parallel MPI stuff to prevent print multiple lines, rank = 0 <=> First processor
	{
		cout << "Save after Reactor" << endl;
		if(i==0) 
			cout << "*Full data AFTER REACTOR saving is activated" << endl;
		// Check if meanSpeciesProdRate.txt file exists at the first step
		// If yes, clear the file content
		if(file_exists("outputs/full_dataParticleAftREACTOR.dat") && i==0)
		{
		
			cout << " -------------- Warrning --------------" << endl;
			cout << " outputs/full_dataParticleAftREACTOR.dat exists. Clearing file ... " << endl;	
			
			ofstream full_dataParticleAftREACTOR_clear("outputs/full_dataParticleAftREACTOR.dat", ios::out | ios::trunc);	//open file in trunc mode to clear the content
			full_dataParticleAftREACTOR_clear.close(); //close the file
	
		}	
		
		ofstream full_dataParticleAftREACTOR("outputs/full_dataParticleAftREACTOR.dat",ios::app); //ios::app = append at the end of the file
		if(full_dataParticleAftREACTOR)
		{
			if(i==0)
			{
   				ofstream full_dataParticleAftREACTOR("outputs/full_dataParticleAftREACTOR.dat"); 
   				full_dataParticleAftREACTOR << "#Time	";
				full_dataParticleAftREACTOR << "Particle_number	";
   				for (int k=0; k<nsp; k++)
      					full_dataParticleAftREACTOR << mixture->speciesName(k) << "	";
				full_dataParticleAftREACTOR << "Hm	";	//HT@2020.08.22 : Need to remove
				full_dataParticleAftREACTOR << "Temperature	" << endl;

				for(int p=0; p<nTot; p++)
				{
     					full_dataParticleAftREACTOR << t << "	";
     					full_dataParticleAftREACTOR << p << "	";	// Particle starts from 1
     					for (int k=0; k<nsp; k++)
        					full_dataParticleAftREACTOR << listParticles[p]->m_Yk_gas[k] << "	";
					full_dataParticleAftREACTOR << listParticles[p]->m_H_gas << "	";	//HT2020.08.22 : Need to remove
					full_dataParticleAftREACTOR << listParticles[p]->m_T_gas << "	" << endl;
				}
			}
			else
			{
				for(int p=0; p<nTot; p++)
				{
     					full_dataParticleAftREACTOR << t << "	";
     					full_dataParticleAftREACTOR << p << "	";	// Particle starts from 1
     					for (int k=0; k<nsp; k++)
        					full_dataParticleAftREACTOR << listParticles[p]->m_Yk_gas[k] << "	";
                                        full_dataParticleAftREACTOR << listParticles[p]->m_H_gas << "	";      //HT2020.08.22 : Need to remove
					full_dataParticleAftREACTOR << listParticles[p]->m_T_gas << "	" << endl;
				}	
			}
		}

		full_dataParticleAftREACTOR.close(); //close the  file
	} //END if(rank==0)

     } //END if(activateFullDataAftEMST)
     // =========== END FULL DATA PARTICLES AFTER EMST ===========



    double varEnthalpyCFD = 0.0; 	// Correction CFD Stochastic - Huu-Tri NGUYEN - 13 Nov 2019
   
     if (flagCFDStochas)	// flagCFD == true, use CFD result correction 
     {
	
	/* === Interpolate the data to match with time step of ORCh - Huu-Tri NGUYEN - 5 Sept. 2019 === */
        vector<double> timeORCh;           // Store all the time of ORCh = Iteration * delta_t
        double  mean_hCFDinterpo[nbLines]; // Store the mean CFD enthalpy interpolated - Size = nb of iterations
	double Mean_Hm_ini; // Huu-Tri NGUYEN - 14.01.2020 - Interpolation

	if(i==0)	// Huu-Tri NGUYEN - 14.01.2020
	{ 
		Mean_Hm_ini = Mean_Hm;	// Mean_Hm of iteration t, we are now t+1 because after Stochastic closure
		if(rank ==0) cout << " Correction Stochastic - Mean ini = " << Mean_Hm_ini << endl;
	}


	mean_hCFDinterpo[0] = Mean_Hm_ini;  // No variation between CFD and ORCh >> mean_hCFDinterpo - Mean_Hm = 0

	for(int step = 0; step < nbLines; step++)	
	{
		timeORCh.push_back(step*delta_t); // size of timeORCh vector equals to nbLines (nb of iterations)
	}
	

	// Find Top, Bot position
	for(int step = 1; step < nbLines; step++)	// Start from second value (step = 1), mean_hCFDinterpo[0] = Mean_Hm
        {
		for(int row = 0; row < maxRow-1; row++)	// End before the last number
		{
			if(timeORCh[step] < t_CFD[0]) // Verify if out of range
			{
                       //         mean_hCFDinterpo[step] = Mean_Hm; // Smaller than the first CFD result
									// assume no loss                                        
				double tTop = 0.0;	// tTop < timeORCh < tBot
				double tBot = t_CFD[0];
				double hTop = mean_hCFDinterpo[0];
				double hBot = mean_hCFD[0];


				mean_hCFDinterpo[step] = hTop + ((hBot - hTop)/(tBot-tTop)*(timeORCh[step]-tTop));			


			}
			else if(timeORCh[step] > t_CFD[maxRow-1]) //  Bigger than the last CFD result, take the last value hm
			{

				mean_hCFDinterpo[step] = mean_hCFDinterpo[step-1];
			} 
			else	// In the range of interpolation
			{
				if(timeORCh[step] > t_CFD[row+1])	// Find the cursor position
				{	
					// Do not thing >> Move to next row
				}
				else if(timeORCh[step] == t_CFD[row])
				{	
					mean_hCFDinterpo[step] = mean_hCFD[row];
				}	
				else if(timeORCh[step] > t_CFD[row] && timeORCh[step] < t_CFD[row+1]) // Interpolate
				{
					double tTop = t_CFD[row];	// tTop < timeORCh < tBot
					double tBot = t_CFD[row+1];
					double hTop = mean_hCFD[row];
					double hBot = mean_hCFD[row+1];
					mean_hCFDinterpo[step] = hTop + ((hBot - hTop)/(tBot-tTop)*(timeORCh[step]-tTop));
				}		
			}
			 	
		}
	} 
	

	// Save the interpolated result to output 1
	ofstream mean_hCFD_interpolated(dirCFDOut.c_str());	// Using string - ifstream need path in form of c_str, not only string


	if(mean_hCFD_interpolated)
	{
		mean_hCFD_interpolated << "#Time(s)	"; 
		mean_hCFD_interpolated << "h_mean_interpolated(J/kg)" << endl;

		for(int step = 0; step < nbLines; step++)
		{
			mean_hCFD_interpolated << timeORCh[step] << "	";
			mean_hCFD_interpolated << mean_hCFDinterpo[step] << endl;
		}
	}
	else
	{
		cout << "ERROR: Unable to write h_tau_Interpolated.txt" << endl;
		cout << "Please check computeMultipleInlet.cpp" << endl;
	}

	mean_hCFD_interpolated.close();


	/* =================================== END Interpolate =================================== */
	
	/* =================== CORRECTION STOCHASITIC EQUATIONS ==================== */
        /* CFD Enthalpy Correction - Huu-Tri NGUYEN - 13 Nov 2019 */
        /* Calculate the variation between hCFD and Mean_Hm for each particle */
        /* Calculate heat-loss weight for each particle */
        /* See more at Stochastic equations part below */

	double Mean_Hm_NextStep = 0.0;

 	for (int p=ndil; p<nTot; p++)
      	{	
         	Mean_Hm_NextStep += (listParticles[p]->m_P_gas_liquid*listParticles[p]->m_H_gas);
      	}
      	Mean_Hm_NextStep /= Total_gas_mass;




       if(timeORCh[i] == i*delta_t)
        {
 		/* Main equation */
		/* hp_new = hp_old + weight*varEnthalpyCFD */

		/** varEnthalpyCFd = hCFD - Mean_Hm **/
		/** weight = hp_old / Mean_Hm **/ //Calculate this later (below) due to Curl mixing closure


                varEnthalpyCFD = mean_hCFDinterpo[i+1] - Mean_Hm_NextStep;	// After Stochastic closure so we have calculate for next time step Mean_Hm_NextStep[i+1]
		
		// Huu-Tri Commented check Correction Stochastic - 2020.03.05
	/*	if(rank==0)	//Print only 1 time on screen
		{
			cout << " ======================== Correction Stochastic ========================= " << endl;
			cout << " ----- At " << i*delta_t << "s -step " << i << " -----" << endl; 
			cout << "hCFD at t+1 " << (i+1)*delta_t << " = " << mean_hCFDinterpo[i+1] << endl;
			cout << "Mean_Hm_NextStep at t+1 " << (i+1)*delta_t << "s = " << Mean_Hm_NextStep << endl;
			cout << "Variation = hCFD[i+1] - Mean_Hm_NextStep = " << varEnthalpyCFD << endl;
		}
	*/
		// Averaged loss on nTot particles
		// varEnthalpyCFD *= nTot;
		//if(rank ==0)
		//	cout << "Variation*nTot = " << varEnthalpyCFD << endl;


        }
        else
        {
                cout << "Something is wrong with interpolation - check computeMultipleInlet.cpp" << endl;
        }




 
        /* ================== END CORRECTION ================= */
	

	/*====================== CORRECTION STOCHASTIC EQUATIONS ====================== */
	/* CFD Enthalpy Correction STOCHASTIC EQUATIONS - Huu-Tri NGUYEN - 13 Novembre 2019*/
	/* Start at i==1 (i==0 is initial condition) */	
	
	double checkEnthalpy = 0.0;
	double sumVar = 0.0;
	double sumH_gas = 0.0;
	int mixingIte = -1; // Number of mixing iterations before enthalpy correction  

	if(i<mixingIte) // Wait "mixingIte" steps before correct enthalpy	//i=0, initial state
	{
		// Without correction			if(rank==0) 
		cout << "The particles mix in " << mixingIte << " iterations = " << (mixingIte*delta_t)*1000 << "ms before correction" << endl << endl;
	}
	else // i> mixingIte
	{
	
	// Commented by Huu-Tri Nguyen - 2020.03.05
	//	if(rank==0) 
	//	{	
	//		cout << "Correction Enthalpy" << endl;

	//	}
	
		// Define Particle type
	//	vector<double> particleType;
	//	for (int p=0; p<nTot; p++)
	//	{
      	//		if (p < nbParticles[0])
        //			particleType[p] = 0;
      	//		else if (p > nbParticles[0]-1 && p < nbParticles[0] + nbParticles[1])
        // 			particleType[p] = 1;
      	//		else if (p >= nbParticles[0] + nbParticles[1])
        // 			particleType[p] = 2;
	//	}


		// Calculate hTot_abs = sum(abs(h_particle))
		double hTot_abs = 0.0;
		double sum_Tm_C = 0.0;	// Huu-Tri Nguyen - 26.01.2020 - Calculate weight
		for (int p=0; p<nTot; p++)
		{

			hTot_abs += abs(listParticles[p]->m_H_gas);
			sum_Tm_C += (listParticles[p]->m_T_gas - 273.15);

		}
	
		// check enthalpy particle
		double Mean_Hm_NextStep_checkBefore = 0.0;
        	for (int p=ndil; p<nTot; p++)
        	{
                	Mean_Hm_NextStep_checkBefore += (listParticles[p]->m_P_gas_liquid*listParticles[p]->m_H_gas);
		
        	}
        
		Mean_Hm_NextStep_checkBefore /= Total_gas_mass;
	
		// Commented by Huu-Tri Nguyen - 2020.03.05
	/*	if(rank ==0) cout << " Mean_Hm_NextStep_check = " << Mean_Hm_NextStep_checkBefore << " :: ndil = " << ndil << " :: Total gas mass = " << Total_gas_mass << " :: nTot = " << nTot << endl; 
		

		if(rank ==0) cout << "sum_Tm_C = " << sum_Tm_C <<  " Mean_Tm_C next step = " << sum_Tm_C/nTot << endl;
	*/

		double sum_weight = 0.0; // HT 26.01.2020
		
		
		for (int p=0; p<nTot; p++)
		{	
			double weight_particle = 0;
			sumVar += varEnthalpyCFD; // Should = varEnthalpy before /nTot
			sumH_gas += listParticles[p]->m_P_gas_liquid*listParticles[p]->m_H_gas;	// Should = meanHm (before correction)
	
			//Correction
		//HT	weight_particle = listParticles[p]->m_H_gas/Mean_Hm_NextStep;

		//HT	weight_particle = abs(listParticles[p]->m_H_gas)/hTot_abs*nTot;	// New formulation - Huu-Tri - 22.01.2020

			weight_particle = (listParticles[p]->m_T_gas-273.15)/sum_Tm_C*nTot; 
			sum_weight += weight_particle;
		//	weight_particle = 1.0/nTot;
		//	if(rank==0) 	cout << "var = " << varEnthalpyCFD << " Weight = " << weight_particle << " * = " << weight_particle*varEnthalpyCFD << endl;

//		if(rank==0) cout << "Particle " << p << " before correction = " << listParticles[p]->m_H_gas << endl; 	
 //		if(p >= nbParticles[0])
 			listParticles[p]->m_H_gas = listParticles[p]->m_H_gas +  weight_particle*varEnthalpyCFD;	// Add abs(weight_particle) - Huu-Tri Nguyen - 17.01.2020	


//		if(rank==0) cout << " 	after correction = " << listParticles[p]->m_H_gas << endl;

			checkEnthalpy += listParticles[p]->m_H_gas; // Should = mean_CFDinterpo (after correction)
//			if(rank==0) cout << " check Enthaly particle " << p << " = " << checkEnthalpy << endl;
		
		}


		sumH_gas /= Total_gas_mass;		
		checkEnthalpy /= nTot;
		
		// Check enthalpy after correction
		double Mean_Hm_NextStep_checkAfter = 0.0;
                for (int p=ndil; p<nTot; p++)
                {
                        Mean_Hm_NextStep_checkAfter += (listParticles[p]->m_P_gas_liquid*listParticles[p]->m_H_gas);

                }
                Mean_Hm_NextStep_checkAfter /= Total_gas_mass;

		double varBeforeAfter = 0.0;
		varBeforeAfter = Mean_Hm_NextStep_checkAfter - Mean_Hm_NextStep_checkBefore;


	/* Commented by Huu-Tri Nguyen - 2020.03.05 
               if(rank ==0) 
		{
			cout << " Mean_Hm_NextStep_check AFTER = " << Mean_Hm_NextStep_checkAfter  << endl;
			cout << " varBeforeAfter = Mean_Hm_NextStep_checkAfter - Mean_Hm_NextStep_checkBefore  = " << varBeforeAfter << "  >> should = varEnthalpyCFD " << endl;
		} 

		if(rank==0)
		{
		cout << " --------- Correction Check ---------- " << endl;
		cout << " sum_weight = " << sum_weight << endl;
		cout << "Before correction, sumH_gas (should = Mean_Hm) =  " << sumH_gas << endl;
		cout << "After correction, checkEnthalpy (should = hCFD) = " << checkEnthalpy << endl;	
		cout << "Variation = checkEnthalpy - mean_hCFDinterpo[i+1] = " <<  checkEnthalpy - mean_hCFDinterpo[i+1] << endl;
		cout << endl;
		}
	*/
	} //End else


	// ==== Check correction of particles enthalpy  - Huu-Tri Nguyen - 16.01.2020 ==== //
	bool hmindataAfter = false;
	if(rank==0 && hmindataAfter)
	{	
		if(file_exists("outputs/Z_hmindataAfterCorrection.dat") && i==0)
		{
			cout << " -------------- Warrning --------------" << endl;
			cout << " outputs/Z_hmindataAfterCorrection.dat exists. Clearing file ... " << endl;	
			
			ofstream Z_hmindataAfter_clear("outputs/Z_hmindataAfterCorrection.dat", ios::out | ios::trunc);	//open file in trunc mode to clear the content
			Z_hmindataAfter_clear.close(); //close the file

		}		

		ofstream Z_hmindataAfter("outputs/Z_hmindataAfterCorrection.dat",ios::app); //ios::app = append at the end of the file
		if(Z_hmindataAfter)
		{
			if(i==0)	// write the first line
			{  
				Z_hmindataAfter << "#1:time  2:Particle_number  3:particle type  4:Zn2p	5:hminp	6:hp_befo	7:hp_after	8:OutOfBound" << endl;
			}		


			double Yn2_0, Yn2_f;
			double Yn2_gbIni = 0, Zn2_gbIni = 0  ; // Huu-Tri Nguyen 16.01.2020 - To calculate Zn2_gbIni >> hmin(Z)
			Yn2_0 = 0.766990291; 	// UMONS case - Inlet air preheated 
			Yn2_f = 0.396761134;	// UMONS case - Inlet fuel
	    
			// Find initial mixture fraction Zn2_gbIni of inlet GB 
			for (int k=0; k<nsp; k++)
			{	
			
				if (mixture->speciesName(k) == "N2")		//Huu-Tri NGUYEN - 07.01.2020
				{
			//GB		Yn2_gbIni = Ym_Trajectories_store[2][0][k];	// Only take the first step i=0
				}
			}

			Zn2_gbIni = ( Yn2_gbIni-Yn2_0)/(Yn2_f-Yn2_0);	// Initial mixture fraction of Burnt gas

	
			// Calculate enthalpy mix of inlets - UMONS case - Huu-Tri Nguyen 15.01.2020
			double h_gasIni = Hm_inletIni[0];
			double h_airIni = Hm_inletIni[1];
			double h_burntGasIni = Hm_inletIni[2];
			double hminParticle = 0;			// Save hmin(Z) of each particle
			double hmaxParticle =0; 			// Save hmax(Z) of each particle - 20.01.2020
			double Zn2Particle = 0;
			//double h_mixZ_Deter;
			//	h_mixZ_Deter = h_gas*Zn2_DeterMix + h_air*(1-Zn2_DeterMix);	// In the case of 2 inlets adidabatic


			// Huu-Tri NGUYEN - 07.01.2020 - Calculate Z (mixture fraction) by N2
			for (int p=0; p<nTot; p++)
      			{	

				
				// Define Particle type
      				double particleType;
      				if (p < nbParticles[0])
         				particleType = 0;
      				else if (p > nbParticles[0]-1 && p < nbParticles[0] + nbParticles[1])
         				particleType = 1;
      				else if (p >= nbParticles[0] + nbParticles[1])
         				particleType = 2;

				double Yn2, Yn2_0, Yn2_f;
				Zn2Particle = 0;
				Yn2_0 = 0.766990291; 	// UMONS case - Inlet air preheated 
				Yn2_f = 0.396761134;	// UMONS case - Inlet fuel
				int OutOfBound = 0; 	// To know if the particle gets out the triangle (=1) or Not (=0)
      
      				for (int k=0; k<nsp; k++)
      				{
	
         				if (mixture->speciesName(k) == "N2")		//Huu-Tri NGUYEN - 07.01.2020
            					Yn2 = listParticles[p]->m_Yk_gas[k];
				}      
       	
				Zn2Particle = (Yn2-Yn2_0)/(Yn2_f-Yn2_0); //Huu-Tri NGUYEN - 07.01.2020: Mixture fraction based on N2 of 1 particle
									// Zst_n2 is calculated in Excel file of case conditions
				// Calculate hmin(Z) of this particle
				if(Zn2Particle <= Zn2_gbIni)		// Left branch air-burnt gas on (h,z) space
				{
					hminParticle = (h_burntGasIni - h_airIni)/Zn2_gbIni*Zn2Particle + h_airIni; 
				}
				else // Zn2Particle > Zn2_gbIni 	// Right branch burnt gas-fuel on (h,z) space
				{
					hminParticle = (h_gasIni - h_burntGasIni)/(1 - Zn2_gbIni)*(Zn2Particle - Zn2_gbIni) + h_burntGasIni;
				}



				// Calculate hmax - 20.01.2020
				hmaxParticle = (h_gasIni - h_airIni)*Zn2Particle + h_airIni;


				// Print to file
				Z_hmindataAfter << t << "	";
				Z_hmindataAfter << p+1 << "	";
				Z_hmindataAfter << particleType << "	";
				Z_hmindataAfter << Zn2Particle << "	";
				Z_hmindataAfter << hminParticle << "	";
				Z_hmindataAfter << listParticles[p]->m_H_gas << "	";	// Particle enthalpy BEFORE the comparison with hmin

				// Check if Enthalpy after correction of this particle < or > hminParticle
				// If < hminParticle, particle go outside the triangle (h_gasIni, h_airIni, h_burntGasIni) 
				// Enthalpy of this particle should be = hminParticle
			//	if(listParticles[p]->m_H_gas < hminParticle)
			//	{
			//		OutOfBound = 1; 	// Particle is out of boudary hmin
	//HT				listParticles[p]->m_H_gas =  hminParticle;
			//	}





				Z_hmindataAfter << listParticles[p]->m_H_gas << "	";	// Particle enthalpy AFTER the comparison with hmin

				Z_hmindataAfter << OutOfBound << "	";
				Z_hmindataAfter << endl;

			} //end for(p = 0>pTot)
		} //End if(Z_hmindataAfter)

	Z_hmindataAfter.close();
	// End Check correction of particles enthalpy  - Huu-Tri Nguyen - 16.01.2020 //
	} // End if(rank==0)

  } //End if(flagCFD)
	/*==================== END CORRECTION ====================*/



	// h_particle = hmin (After Stochastic closure) - Huu-Tri Nguyen 20.01.2020
	// For iteration i+1 >> Calculate Mean_Hm before go into Lagrangian closure
	// This part should be place here because when i==0, all particles has the initial state so h_particle[i==0] = hmin[i==0] = hInletIni 

			double Yn2_0, Yn2_f;
			double Yn2_gbIni = 0, Zn2_gbIni = 0  ; // Huu-Tri Nguyen 16.01.2020 - To calculate Zn2_gbIni >> hmin(Z)
			Yn2_0 = 0.766990291; 	// UMONS case - Inlet air preheated 
			Yn2_f = 0.396761134;	// UMONS case - Inlet fuel
	    
			// Find initial mixture fraction Zn2_gbIni of inlet GB 
			for (int k=0; k<nsp; k++)
			{	
			
				if (mixture->speciesName(k) == "N2")		//Huu-Tri NGUYEN - 07.01.2020
				{
					Yn2_gbIni = Ym_Trajectories_store[2][0][k];	// Only take the first step i=0
				}
			}

			Zn2_gbIni = ( Yn2_gbIni-Yn2_0)/(Yn2_f-Yn2_0);	// Initial mixture fraction of Burnt gas

	
			// Calculate enthalpy mix of inlets - UMONS case - Huu-Tri Nguyen 15.01.2020
			double h_gasIni = Hm_inletIni[0];
			double h_airIni = Hm_inletIni[1];
			double h_burntGasIni = Hm_inletIni[2];
			double hminParticle = 0;			// Save hmin(Z) of each particle
			double hmaxParticle =0; 			// Save hmax(Z) of each particle - 20.01.2020
			double Zn2Particle = 0;
			//double h_mixZ_Deter;
			//	h_mixZ_Deter = h_gas*Zn2_DeterMix + h_air*(1-Zn2_DeterMix);	// In the case of 2 inlets adidabatic


			// Huu-Tri NGUYEN - 07.01.2020 - Calculate Z (mixture fraction) by N2
			for (int p=0; p<nTot; p++)
      			{
				
				double Yn2, Yn2_0, Yn2_f;
				Zn2Particle = 0;
				Yn2_0 = 0.766990291; 	// UMONS case - Inlet air preheated 
				Yn2_f = 0.396761134;	// UMONS case - Inlet fuel
				int OutOfBound = 0; 	// To know if the particle gets out the triangle (=1) or Not (=0)
      
      				for (int k=0; k<nsp; k++)
      				{
	
         				if (mixture->speciesName(k) == "N2")		//Huu-Tri NGUYEN - 07.01.2020
            					Yn2 = listParticles[p]->m_Yk_gas[k];
				}      
       	
				Zn2Particle = (Yn2-Yn2_0)/(Yn2_f-Yn2_0); //Huu-Tri NGUYEN - 07.01.2020: Mixture fraction based on N2 of 1 particle
									// Zst_n2 is calculated in Excel file of case conditions
				// Calculate hmin(Z) of this particle
				if(Zn2Particle <= Zn2_gbIni)		// Left branch air-burnt gas on (h,z) space
				{
					hminParticle = (h_burntGasIni - h_airIni)/Zn2_gbIni*Zn2Particle + h_airIni; 
				}
				else // Zn2Particle > Zn2_gbIni 	// Right branch burnt gas-fuel on (h,z) space
				{
					hminParticle = (h_gasIni - h_burntGasIni)/(1 - Zn2_gbIni)*(Zn2Particle - Zn2_gbIni) + h_burntGasIni;
				}



				// Calculate hmax - 20.01.2020
				hmaxParticle = (h_gasIni - h_airIni)*Zn2Particle + h_airIni;


				// Check if Enthalpy after correction of this particle < or > hminParticle
				// If < hminParticle, particle go outside the triangle (h_gasIni, h_airIni, h_burntGasIni) 
				// Enthalpy of this particle should be = hminParticle
			//	if(listParticles[p]->m_H_gas < hminParticle)
			//	{
			//		OutOfBound = 1; 	// Particle is out of boudary hmin
	//HT				listParticles[p]->m_H_gas =  hminParticle;
			//	}


			} //end for(p = 0>pTot)

	//  END h_particle = hmin (After Stochastic closure) - Huu-Tri Nguyen 20.01.2020

    // Calculate right hand-side of Y(=mixing + source term) and T (=source term) - Huu-Tri NGUYEN - 2019.12.05
   	// Initialization
	double *Tm_gas_after = new double[nTot];	// 1D array temperature for each particle
	double *rightSide_T  = new double[nTot];	// 1D array right-hand-side for each particle

	double **Yk_gas_after  = new double*[nTot];	// 2D array ** with first dimension (number of particles)
	for (int p=0; p<nTot; p++)	// Create second dimension 
	{
		Yk_gas_after[p] = new double[nsp];	// nsp = number of species has declared above

	}
	double **rightSide  = new double*[nTot];	// 2D array ** with first dimension (number of particles)
	for (int p=0; p<nTot; p++)	// Create second dimension 
	{
		rightSide[p] = new double[nsp];	// nsp = number of species has declared above
	}

	double *meanSourceTerm_Stochas = new double[nsp];
	double meanSourceTerm_T = 0.0;


    if(activateSourceTermParticle)
    { 
	// Calculation
	// For each particle
	for (int p=ndil; p<nTot; p++)
      	{
		Tm_gas_after[p] = listParticles[p]->m_T_gas;		
		rightSide_T[p] = Tm_gas_after[p] - Tm_gas_before[p];
		rightSide_T[p] /= delta_t;
	
 		for (int k=0; k<nsp; k++)
        	{
              	 	Yk_gas_after[p][k] = listParticles[p]->m_Yk_gas[k];
			rightSide[p][k] = Yk_gas_after[p][k] - Yk_gas_before[p][k];	// HT@2020.09.17: Correct after - before
			rightSide[p][k] /= delta_t;
	
		}	
		
	}
	
	// Sum and Averaged on all particles
	for (int p=ndil; p<nTot; p++)
      	{
		meanSourceTerm_T += rightSide_T[p];
		for (int k=0; k<nsp; k++)
        	{
			meanSourceTerm_Stochas[k] += rightSide[p][k]; // rightSide = Mix + w (Must subtract Mix) 
		}
	}
	

	meanSourceTerm_T /= nTot;
	
	for (int k=0; k<nsp; k++)
	{
		meanSourceTerm_Stochas[k] /= nTot;
			//if(rank==0) // commented by Huu-Tri@20200731
				//cout << "meanSource term of " << mixture->speciesName(k) << " = " << meanSourceTerm_Stochas [k]  << endl;
	}
	

	// Step 3: Print the MEAN species production rate to CFD_results/meanSpeciesProdRate.txt

	if(rank==0)	//Parallel MPI stuff to prevent print multiple lines, rank = 0 <=> First processor
	{

		// Check if meanSpeciesProdRate.txt file exists at the first step
		// If yes, clear the file content
		if(file_exists("outputs/meanSpeciesProdRate.txt") && i==0)
		{
		
			cout << " -------------- Warrning --------------" << endl;
			cout << " outputs/meanSpeciesProdRate.txt exists. Clearing file ... " << endl;	
			
			ofstream meanSpeciesProdRate_clear("outputs/meanSpeciesProdRate.txt", ios::out |  ios::trunc);	//open file in trunc mode to clear the content
			meanSpeciesProdRate_clear.close(); //close the file
	
		}		


	
		ofstream meanSpeciesProdRate("outputs/meanSpeciesProdRate.txt",ios::app); //ios::app = append at the end of the file
		if(meanSpeciesProdRate)
		{
			if(i==0)	//First step: Need to write the headline
			{
				// First line
				meanSpeciesProdRate << "Time	"; 
					//meanSpeciesProdRate << "T	"; //HT@2020.09.08
				for (int k=0; k<nsp; k++)
				{
				//	meanSpeciesProdRate << k+3 << ":" << mixture->speciesName(k) << "	";
				meanSpeciesProdRate << mixture->speciesName(k) << "        "; //Added by Huu-Tri@2020.09.08
				}	
				meanSpeciesProdRate << "T	"; 
				meanSpeciesProdRate << endl;			

				// Data from ORCh
				meanSpeciesProdRate << i*delta_t << "	";
					//meanSpeciesProdRate << meanSourceTerm_T << "	";
				for (int k=0; k<nsp; k++)
				{
					meanSpeciesProdRate << meanSourceTerm_Stochas[k] << "	";
				}	
				meanSpeciesProdRate << meanSourceTerm_T << "    "; //HT@2020.09.08
				meanSpeciesProdRate << endl;
			}
			else
			{		
				// Data from ORCh
				meanSpeciesProdRate << i*delta_t << "	";
					//meanSpeciesProdRate << meanSourceTerm_T << "	"; //HT@2020.09.08
				for (int k=0; k<nsp; k++)
				{
					meanSpeciesProdRate << meanSourceTerm_Stochas[k]  << "	";
				}
				meanSpeciesProdRate << meanSourceTerm_T << "    ";
				meanSpeciesProdRate << endl;
			}
		}
		else
		{	
			cout << "ERROR: Impossible to write meanSpeciesProdRate.txt" << endl;
			cout << "Please check computeMultipleInlet.cpp" << endl;
		}
	
		meanSpeciesProdRate.close();
	}

	// Step 4: Print the source term of each particle to CFD_results/dataSourceTerm_Particle.txt
	if(rank==0)	//Parallel MPI stuff to prevent print multiple lines, rank = 0 <=> First processor
	{

		// Check if dataSourceTerm_Particle.txt file exists at the first step
		// If yes, clear the file content
		if(file_exists("outputs/dataSourceTerm_Particle.txt") && i==0)
		{
		
			cout << " -------------- Warrning --------------" << endl;
			cout << " outputs/dataSourceTerm_Particle.txt exists. Clearing file ... " << endl;	
			
			ofstream dataSourceTerm_Particle_clear("outputs/dataSourceTerm_Particle.txt", ios::out |  ios::trunc);	//open file in trunc mode to clear the content
			dataSourceTerm_Particle_clear.close(); //close the file
	
		}		


	
		ofstream dataSourceTerm_Particle("outputs/dataSourceTerm_Particle.txt",ios::app); //ios::app = append at the end of the file
		if(dataSourceTerm_Particle)
		{
			if(i==0)	//First step: Need to write the headline
			{
				// First line
				dataSourceTerm_Particle << "#Time	"; 
				dataSourceTerm_Particle << "Particle_number	"; 
					//dataSourceTerm_Particle << "T	"; //HT@2020.09.08
				for (int k=0; k<nsp; k++)
				{
				//	dataSourceTerm_Particle << k+4 << ":" << mixture->speciesName(k) << "	";
				dataSourceTerm_Particle << mixture->speciesName(k) << "	"; // Added by Huu-Tri@2020.09.08
				}
				dataSourceTerm_Particle << "T";
				dataSourceTerm_Particle << endl;			

				// Data from ORCh
				for(int p=0; p<nTot; p++)
				{
					dataSourceTerm_Particle << (i+1)*delta_t << "	";	//rightSide is calculated on (i+1)-i
					dataSourceTerm_Particle << p << "	";
						//dataSourceTerm_Particle << rightSide_T[p] << "	"; //HT@2020.09.08
					for (int k=0; k<nsp; k++)
					{
						dataSourceTerm_Particle << rightSide[p][k] << "	";
					}
					dataSourceTerm_Particle << rightSide_T[p] << "  ";
					dataSourceTerm_Particle << endl;
				}
			}
			else
			{		
				// Data from ORCh
				for(int p=0; p<nTot; p++)
				{
					dataSourceTerm_Particle << (i+1)*delta_t << "	";
					dataSourceTerm_Particle << p << "	";
						//dataSourceTerm_Particle << rightSide_T[p] << "	"; //HT@2020.09.08
					for (int k=0; k<nsp; k++)
					{
						dataSourceTerm_Particle << rightSide[p][k] << "	";
					}
					dataSourceTerm_Particle << rightSide_T[p] << "  ";
					dataSourceTerm_Particle << endl;
				}
			}
		}
		else
		{	
			cout << "ERROR: Impossible to write dataSourceTerm_Particle.txt" << endl;
			cout << "Please check computeMultipleInlet.cpp" << endl;
		}
	
		dataSourceTerm_Particle.close();
	}

   } // End if(activateSourceTermParticle)


	// Free w_species memory (array pointer should be deleted after use)	
	delete Tm_gas_before;
	delete Tm_gas_after;
	delete rightSide_T; 

	for (int p=0; p<nTot; p++)	// Create second dimension 
	{
		delete[] Yk_gas_before[p];
	}
	delete[] Yk_gas_before;

	for (int p=0; p<nTot; p++)	// Create second dimension 
	{
		delete[] Yk_gas_after[p];
	}
	delete[] Yk_gas_after;

	for (int p=0; p<nTot; p++)	// Create second dimension 
	{
		delete[] rightSide[p];
	}
	delete[] rightSide;
	
	delete meanSourceTerm_Stochas ;
	
   // End Calculate right hand-side of Y(=mixing + source term) and T (=source term) - Huu-Tri NGUYEN - 2019.12.05


      // Store time step	
      t = t + dt;
      time_store[i+1] = t;

      // Print to screen
       if (rank == 0) cout << '.' << flush;
       //if (rank ==0)  cout << t << " |*| " << flush;
   } // END of big loop i

    // FLAGANN: Release OPENCV MAT >> If declare before for(i), it should release after loop
    
    // END FLAGANN: Release OPENCV MAT


} //end of main computeMultipleInlet::getMultipleInlet







void computeMultipleInlet::Reacting(vector<Particle*> &listParticles, int nsp, double dt, double Pressure)
{
   int nTot = listParticles.size();
   vector<double> Ym(nsp,0.0);
   double Hm;
   double Zm;
   double Tm;

   for (int p=0; p<nTot; p++)
   {
      for (int k=0; k<nsp; k++)
         Ym[k] = listParticles[p]->m_Yk_gas[k];

      Hm = listParticles[p]->m_H_gas;

      Next_Time_Step(Pressure, Ym, Hm, Tm, dt);

      for (int k=0; k<nsp; k++)
         listParticles[p]->m_Yk_gas[k] = Ym[k];

      listParticles[p]->m_H_gas = Hm;
      listParticles[p]->m_T_gas = Tm;
   }
}

void computeMultipleInlet::ReactingParallel(vector<Particle*> &listParticles, int nsp, double dt, double Pressure)
{
   int nTot = listParticles.size();
   vector<double> Ym(nsp,0.0);
   double Hm;
   double Tm;

   for (int p=Ifi_rank; p<Ila_rank; p++)
   {
      for (int k=0; k<nsp; k++)
         Ym[k] = listParticles[p]->m_Yk_gas[k];

      Hm = listParticles[p]->m_H_gas;

      Next_Time_Step(Pressure, Ym, Hm, Tm, dt);

      for (int k=0; k<nsp; k++)
         listParticles[p]->m_Yk_gas[k] = Ym[k];

      listParticles[p]->m_H_gas = Hm;
      listParticles[p]->m_T_gas = Tm;
   }

   double Data_Proc[nb_var_loc];
   double Data_All[nTot*(nsp+2)];

   int count = 0;
   for (int p=Ifi_rank; p<Ila_rank; p++)
   {
      for (int k=0; k<nsp; k++)
      {
         Data_Proc[count] = listParticles[p]->m_Yk_gas[k];
         count += 1;
      }
      Data_Proc[count] = listParticles[p]->m_H_gas;
      count += 1;
       
      Data_Proc[count] = listParticles[p]->m_T_gas;
      count += 1;
   }

   MPI_Allgatherv(Data_Proc, nb_var_loc, MPI_DOUBLE, Data_All, RecvCounts, Disp, MPI_DOUBLE, MPI_COMM_WORLD);

   count = 0;
   for (int p=0; p<Ifi_rank; p++)
   {
 
      for (int k=0; k<nsp; k++)
      {
         listParticles[p]->m_Yk_gas[k] = Data_All[count];
         count += 1;
      }
      listParticles[p]->m_H_gas = Data_All[count];
      count += 1;

      listParticles[p]->m_T_gas = Data_All[count];
      count += 1;
   }
    //
	count += nb_var_loc;
	//
	for (int p=Ila_rank; p<nTot; p++)
	{
		for (int k=0; k<nsp; k++)
		{
			listParticles[p]->m_Yk_gas[k] = Data_All[count];
			count += 1;
		}
      listParticles[p]->m_H_gas = Data_All[count];
      count += 1;

		listParticles[p]->m_T_gas = Data_All[count];
		count += 1;
	}
}

void computeMultipleInlet::getMixedGasesComposition(vector<MultipleInlet*> listInlets, string step)
{
	int rank, nproc;
   if (step != "Optimisation")
   {
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	} else {
		rank = -1;
   }

   int nsp = mixture->nSpecies();
   int nbInlets = listInlets.size(); 

	vector<double> Compo_Yk_mixed(nsp,0.0);
   double Compo_H_mixed = 0.0;
   double Total_flowRate = 0.0;

   //Composition of the burned gases
   for (int n=0; n<nbInlets; n++)
   {
      if (listInlets[n]->m_X_Species != "")
      {
         if (rank == 0)
          //  cout << "Set the mole fraction of inlet " << n << endl;
			mixture->setState_TPX(listInlets[n]->m_Temperature, listInlets[n]->m_Pressure, listInlets[n]->m_X_Species);
		} else if (listInlets[n]->m_Y_Species != "") {
          //  cout << "Set the mass fraction of inlet " << n << endl;
			mixture->setState_TPY(listInlets[n]->m_Temperature, listInlets[n]->m_Pressure, listInlets[n]->m_Y_Species);
      }

		vector<double> Ym(nsp,0.0);
      double Hm = 0.0;
		mixture->getMassFractions(&Ym[0]);
		Hm = mixture->enthalpy_mass();

		for (int k=0; k<nsp; k++) Compo_Yk_mixed[k] += listInlets[n]->m_flowRate*Ym[k];
      Compo_H_mixed += listInlets[n]->m_flowRate*Hm;
      Total_flowRate += listInlets[n]->m_flowRate;
   }
   
   for (int k=0; k<nsp; k++)
   {
      Compo_Yk_mixed[k] /= Total_flowRate;
   }
   Compo_H_mixed /= Total_flowRate;
   if (rank == 0)
   {
      cout << endl <<  "Composition to enter for the equilibrium computation to get the Burned gases" << endl;
      cout << "Compo_H_mixed " << Compo_H_mixed << endl;
   }

	mixture->setMassFractions(&Compo_Yk_mixed[0]);
	mixture->setState_HP(Compo_H_mixed, listInlets[0]->m_Pressure);

	vector<double> Xm(nsp);
   double T_mixed = 0.0;
	mixture->getMoleFractions(&Xm[0]);
	T_mixed = mixture->temperature();

   if (rank == 0)
   {
      for (int k=0; k<nsp; k++)
      {
         if (Xm[k] != 0.0)
            cout << "X_" << mixture->speciesName(k) << ": " << Xm[k] << endl;
      }
      cout << "T_mixed " << T_mixed << endl;
   }
}

void computeMultipleInlet::Next_Time_Step_with_drgep(vector<bool> Targets, double P, vector<double> &Ym, double &Hm, double &Tm, double delta_t, 
   vector<vector<double> > &R_AD_Trajectories, vector<vector<double> > &max_j_on_Target, string step)
{
	try {
		mixture->setMassFractions(&Ym[0]);
   mixture->setState_HP(Hm, P);

   ConstPressureReactor reac;
   reac.insert(*mixture);
   ReactorNet sim;
   sim.addReactor(reac);

   sim.advance(delta_t);
	} catch (...) {
		return;
	}

   Hm  = mixture->enthalpy_mass();
   Tm = mixture->temperature();
	mixture->getMassFractions(&Ym[0]);

   drgep *species_relations = new drgep();
	species_relations->drgep_0D_species(mixture, Targets, R_AD_Trajectories, 0, 0.0);

  // Comment the duplicated part (???) - Huu-TriNGUYEN@2020.10.14  	
   /*if (step == "DRGEP_Reactions")
   drgep *species_relations = new drgep();
   species_relations->drgep_0D_species(mixture, Targets, R_AD_Trajectories, n, time);

   if (step == "DRGEP_Reactions")
   drgep *species_relations = new drgep();
   species_relations->drgep_0D_species(mixture, Targets, R_AD_Trajectories, n, time);

   if (step == "DRGEP_Reactions")
   drgep *species_relations = new drgep();
   species_relations->drgep_0D_species(mixture, Targets, R_AD_Trajectories, n, time);*/

   if (step == "DRGEP_Reactions")
   {
      int nsp = mixture->nSpecies();
      int nreac = mixture->nReactions();

      vector<vector<double> > rj_for_k (nsp, vector<double> (nreac,0.0));
      species_relations->drgep_0D_reactions(mixture, rj_for_k);

      for (int ka=0; ka<nsp; ka++)
      {
         for (int kb=0; kb<nsp; kb++)
         {
            for (int j=0; j<nreac; j++)
            {
               if (max_j_on_Target[ka][j] < R_AD_Trajectories[ka][kb]*rj_for_k[kb][j])
                  max_j_on_Target[ka][j] = R_AD_Trajectories[ka][kb]*rj_for_k[kb][j];
            }
         }
      }
   }
}

//Next_Time_Step without drgep analysis (to use for the computations with optimisation)
void computeMultipleInlet::Next_Time_Step(double P, vector<double> &Ym, double &Hm, double &Tm, double delta_t)
{
	try {
		mixture->setMassFractions(&Ym[0]);
   mixture->setState_HP(Hm, P);

   ConstPressureReactor reac;
   reac.insert(*mixture);
   ReactorNet sim;
   sim.addReactor(reac);

 
//   cout << "68 " << sim.componentName(68) << endl;	// Huu-Tri Nguyen - To find out the stiff component that CVODE cannot find the time step to integrate >> Cantera error: Components with largest weighted error estimates: 48: -0.250213
//   cout << "71 " << sim.componentName(71) << endl;

   sim.advance(delta_t);
	} catch (...) {
		return;
	}

   // Function to get the internal step number during advance() - Huu-Tri@2020.09.15
   //int numStep = 0;
   //sim.getNumInternalStepReactor(numStep);
   //cout << "Internal step number = " << numStep << endl;

   Hm  = mixture->enthalpy_mass();
   Tm = mixture->temperature();
	mixture->getMassFractions(&Ym[0]);
}

//Next_Time_Step with QSS analysis 
void computeMultipleInlet::Next_Time_Step(double P, vector<double> &Ym, double &Hm, double &Tm, double delta_t,
                    vector<vector<vector<double> > > &Production_Trajectories_ref, vector<vector<vector<double> > > &Consumption_Trajectories_ref, int nInlet, int nLine)
{
	mixture->setMassFractions(&Ym[0]);
   mixture->setState_HP(Hm, P);

   ConstPressureReactor reac;
   reac.insert(*mixture);
   ReactorNet sim;
   sim.addReactor(reac);

   sim.advance(delta_t);

   Hm  = mixture->enthalpy_mass();
   Tm = mixture->temperature();
	mixture->getMassFractions(&Ym[0]);

   int nreac = mixture->nReactions();
   int nsp = mixture->nSpecies();

   double* fwdRates = new double[nreac];
   double* revRates = new double[nreac];
   mixture->getFwdRatesOfProgress(fwdRates);
   mixture->getRevRatesOfProgress(revRates);

   for (int k=0; k<nsp; k++)
   {
      double omega_k_prod = 0.0;
      double omega_k_cons = 0.0;
      for (int j=0; j<nreac; j++)
      {
         omega_k_prod += mixture->productStoichCoeff(k,j)*fwdRates[j]
                        +mixture->reactantStoichCoeff(k,j)*revRates[j];
         omega_k_cons += mixture->reactantStoichCoeff(k,j)*fwdRates[j]
                        +mixture->productStoichCoeff(k,j)*revRates[j];
      }
      Production_Trajectories_ref[nInlet][nLine][k] = omega_k_prod;
      Consumption_Trajectories_ref[nInlet][nLine][k] = omega_k_cons;
   }
   delete mixture;
}


computeMultipleInlet::~computeMultipleInlet() //Destructeur
{}

