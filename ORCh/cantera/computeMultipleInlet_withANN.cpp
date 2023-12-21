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

bool getFileContent(std::string fileName, std::vector<std::string> &vecOfStrs) // Read from file to vector<string>
{
    // Open the File
    std::ifstream in(fileName.c_str());
    // Check if object is valid
    if(!in)
    {
        std::cerr << "Cannot open the File : "<<fileName<<std::endl;
        return false;
    }
    std::string str;
    // Read the next line from File untill it reaches the end.
    while(std::getline(in, str,',')) // Comma delimited
    {
        // Line contains string of length > 0 then save it in vector
        if(str.size() > 0)
            vecOfStrs.push_back(str);
    }
    //Close The File
    in.close();
    return true;
}

void readFromCommaDelimitedFile_Float(std::string fileName, std::vector<float> &floatVector)   // For "row type seperated by comma" file : Read from file to vector<float>
{
    // Open the File
    std::ifstream in(fileName.c_str());
    // Check if object is valid
    if(!in)
    {
        std::cerr << "Cannot open the File : "<< fileName <<std::endl;
        std::cout << "Cannot open the File : "<< fileName <<std::endl;
    }
    std::string str;
    std::vector<std::string> vecOfStrs;
    // Read the next line from File untill it reaches the end.
    while(std::getline(in, str,',')) // Comma delimited
    {
        // Line contains string of length > 0 then save it in vector
        if(str.size() > 0)
            vecOfStrs.push_back(str);
    }
    
    // Convert string to Float
    floatVector.reserve(vecOfStrs.size()); // Secure the size of vector
    std::transform(vecOfStrs.begin(), vecOfStrs.end(), back_inserter(floatVector),
              [](string const& val) {return stod(val);});
    
    //Close The File
    in.close();
    
}

void readFromCommaDelimitedFile_Double(std::string fileName, std::vector<double> &doubleVector)   // For "row type seperated by comma" file : Read from file to vector<double>
{
    // Open the File
    std::ifstream in(fileName.c_str());
    // Check if object is valid
    if(!in)
    {
        std::cerr << "Cannot open the File : "<< fileName <<std::endl;
        std::cout << "Cannot open the File : "<< fileName <<std::endl;
    }
    std::string str;
    std::vector<std::string> vecOfStrs;
    // Read the next line from File untill it reaches the end.
    while(std::getline(in, str,',')) // Comma delimited
    {
        // Line contains string of length > 0 then save it in vector
        if(str.size() > 0)
            vecOfStrs.push_back(str);
    }
    
    // Convert string to Double
    doubleVector.reserve(vecOfStrs.size()); // Secure the size of vector
    std::transform(vecOfStrs.begin(), vecOfStrs.end(), back_inserter(doubleVector),
                   [](string const& val) {return stod(val);});
    
    //Close The File
    in.close();
   
}

void advanceANN(std::vector<float> &child_MeanINPUTLoad, std::vector<float> &child_StdINPUTLoad,
                std::vector<float> &child_MeanOUTPUTLoad, std::vector<float> &child_StdOUTPUTLoad,
                std::vector<float> &child_MaxOUTPUTLoad, std::vector<float> &child_MinOUTPUTLoad,
                std::vector<double> &maxAftREACTORLoad, std::vector<double>  &minAftREACTORLoad,
                cv::PCA pcaANN,
                Tensor &inputPCAANNRegLoad, Model &modelANNRegLoad, Tensor &outputStandardizeANNRegLoad,
                double &Tm, std::vector<double> &Ym,
                int &numVarANN, int &nsp,
                std::vector<float> &input_childLocalStandardized, cv::Mat input_childLocalStandardized_Mat,
                cv::Mat inputPCA, std::vector<float> &inputPCA_Vec,
                std::vector<float> &outputStandardizeANN_Vec,
                IdealGasMix *mixture)
{
    // LOCAL INPUT Standardize
    // Standardize T
    input_childLocalStandardized[numVarANN-1] = Tm - child_MeanINPUTLoad[numVarANN-1];
    input_childLocalStandardized[numVarANN-1] /= child_StdINPUTLoad[numVarANN-1];
    // Standardize Y
    for (int kANN=0; kANN<nsp-1; kANN++) // Don't take N2 (last species)
    {
        input_childLocalStandardized[kANN] = Ym[kANN] - child_MeanINPUTLoad[kANN];
        input_childLocalStandardized[kANN] /= child_StdINPUTLoad[kANN];
        //Cross-check if there is N2 inside, if yes, should move N2 to the end of the scheme
        if(mixture->speciesName(kANN) == "N2") {cout << "Warning : N2 in the ANN!!!" << endl;}
    }
    
    
    // PCA Transform
    // Load to matrix
    input_childLocalStandardized_Mat = cv::Mat(1, numVarANN, CV_32FC1, input_childLocalStandardized.data());
    // Project to LPCA space: inputPCA = input of ANN in LPCA space
    pcaANN.project(input_childLocalStandardized_Mat, inputPCA);
    // Copy matrix inputPCA to inputPCA_Vec vector - cppflow cannot communicate with cv::Mat of OpenCV
    if (inputPCA.isContinuous()) // Continuous memory block
    {
        inputPCA_Vec.assign((float*)inputPCA.data,(float*)inputPCA.data + inputPCA.total()*inputPCA.channels());
    } // General formular for multiple channels matrix, here we have only 1 channel (image RGB = 3 channels)
    
    
    // ANN regression to predict Standardized variation
    // Set input tensor for cppflow (Tensorflow)
    inputPCAANNRegLoad.set_data(inputPCA_Vec);
    // Predict the  Standardized variation - saved as outputANN tensor
    modelANNRegLoad.run(inputPCAANNRegLoad,outputStandardizeANNRegLoad);
    // Get value of outputPCAANNReg tensor (ANN output) to Vector outputStandardizeANN_Vec
    int kRun=0;
    for (float f : outputStandardizeANNRegLoad.get_data<float>())
    {
        outputStandardizeANN_Vec[kRun] = f;
        kRun++;
    }
    
    // INVERSE outputStandardizeANN to outputANN (variation of Y,T in delta_t)
    //  X = standardizedX*standardDeviationScale + meanScale
    // T
    outputStandardizeANN_Vec[numVarANN-1] = outputStandardizeANN_Vec[numVarANN-1] * child_StdOUTPUTLoad[numVarANN-1];
    outputStandardizeANN_Vec[numVarANN-1] = outputStandardizeANN_Vec[numVarANN-1] + child_MeanOUTPUTLoad[numVarANN-1];
    // Verify if output (Composition space) out-of-bound
    outputStandardizeANN_Vec[numVarANN-1] = std::max(outputStandardizeANN_Vec[numVarANN-1],child_MinOUTPUTLoad[numVarANN-1]);
    outputStandardizeANN_Vec[numVarANN-1] = std::min(outputStandardizeANN_Vec[numVarANN-1],child_MaxOUTPUTLoad[numVarANN-1]);
    // Calculate Temperature and bound on Max Min GLOBAL data
    Tm = Tm + outputStandardizeANN_Vec[numVarANN-1];
    Tm = std::min(std::max(Tm,minAftREACTORLoad[numVarANN-1]),maxAftREACTORLoad[numVarANN-1]); // Check if <minValue or >maxValue on GLOBAL
    
    //  Y
    for (int kANN=0; kANN<nsp-1; kANN++) // Don't take N2 (last species)
    {
        outputStandardizeANN_Vec[kANN] = outputStandardizeANN_Vec[kANN] * child_StdOUTPUTLoad[kANN];
        outputStandardizeANN_Vec[kANN] = outputStandardizeANN_Vec[kANN] + child_MeanOUTPUTLoad[kANN];
        // Verify if output (Composition space) out-of-bound
        outputStandardizeANN_Vec[kANN] = std::max(outputStandardizeANN_Vec[kANN],child_MinOUTPUTLoad[kANN]);
        outputStandardizeANN_Vec[kANN] = std::min(outputStandardizeANN_Vec[kANN],child_MaxOUTPUTLoad[kANN]);
        // Calculate Y and bound on Max Min GLOBAL data
        Ym[kANN] = Ym[kANN] + outputStandardizeANN_Vec[kANN];
        Ym[kANN] = std::min(std::max(Ym[kANN],minAftREACTORLoad[kANN]),maxAftREACTORLoad[kANN]); // Check if <minValue or >maxValue
        //Cross-check if there is N2 inside, if yes, should move N2 to the end of the scheme
        if(mixture->speciesName(kANN) == "N2")
            cout << "Warning : N2 in the ANN!!!" << endl;
    }
}
// END Declare function for ANN - Huu-Tri 20210212

// HuuTri@20211006: advanceANN without PCA test
void advanceANN_withoutPCA(std::vector<float> &child_MeanINPUTLoad, std::vector<float> &child_StdINPUTLoad,
                std::vector<float> &child_MeanOUTPUTLoad, std::vector<float> &child_StdOUTPUTLoad,
                std::vector<float> &child_MaxOUTPUTLoad, std::vector<float> &child_MinOUTPUTLoad,
                std::vector<double> &maxAftREACTORLoad, std::vector<double>  &minAftREACTORLoad,
                Tensor &inputPCAANNRegLoad, Model &modelANNRegLoad, Tensor &outputStandardizeANNRegLoad,
                double &Tm, std::vector<double> &Ym,
                int &numVarANN, int &nsp,
                std::vector<float> &input_childLocalStandardized, 
                std::vector<float> &outputStandardizeANN_Vec,
                IdealGasMix *mixture)
{
    // LOCAL INPUT Standardize
    // Standardize T
    input_childLocalStandardized[numVarANN-1] = Tm - child_MeanINPUTLoad[numVarANN-1];
    input_childLocalStandardized[numVarANN-1] /= child_StdINPUTLoad[numVarANN-1];
    // Standardize Y
    for (int kANN=0; kANN<nsp-1; kANN++) // Don't take N2 (last species)
    {
        input_childLocalStandardized[kANN] = Ym[kANN] - child_MeanINPUTLoad[kANN];
        input_childLocalStandardized[kANN] /= child_StdINPUTLoad[kANN];
        //Cross-check if there is N2 inside, if yes, should move N2 to the end of the scheme
        if(mixture->speciesName(kANN) == "N2") {cout << "Warning : N2 in the ANN!!!" << endl;}
    }
    
    
    // PCA Transform
    // Load to matrix
    //HT@20211006 input_childLocalStandardized_Mat = cv::Mat(1, numVarANN, CV_32FC1, input_childLocalStandardized.data());
    // Project to LPCA space: inputPCA = input of ANN in LPCA space
    //HT@20211006 pcaANN.project(input_childLocalStandardized_Mat, inputPCA);
    // Copy matrix inputPCA to inputPCA_Vec vector - cppflow cannot communicate with cv::Mat of OpenCV
    //HT@20211006 if (inputPCA.isContinuous()) // Continuous memory block
    //HT@20211006 {
    //HT@20211006    inputPCA_Vec.assign((float*)inputPCA.data,(float*)inputPCA.data + inputPCA.total()*inputPCA.channels());
    //HT@20211006 } // General formular for multiple channels matrix, here we have only 1 channel (image RGB = 3 channels)
    



    
    // ANN regression to predict Standardized variation
    // Set input tensor for cppflow (Tensorflow)
    //HT@20211006 inputPCAANNRegLoad.set_data(inputPCA_Vec);
    inputPCAANNRegLoad.set_data(input_childLocalStandardized); // HT@20211006: without PCA

    // Predict the  Standardized variation - saved as outputANN tensor
    modelANNRegLoad.run(inputPCAANNRegLoad,outputStandardizeANNRegLoad);
    // Get value of outputPCAANNReg tensor (ANN output) to Vector outputStandardizeANN_Vec
    int kRun=0;
    for (float f : outputStandardizeANNRegLoad.get_data<float>())
    {
        outputStandardizeANN_Vec[kRun] = f;
        kRun++;
    }
    
    // INVERSE outputStandardizeANN to outputANN (variation of Y,T in delta_t)
    //  X = standardizedX*standardDeviationScale + meanScale
    // T
    outputStandardizeANN_Vec[numVarANN-1] = outputStandardizeANN_Vec[numVarANN-1] * child_StdOUTPUTLoad[numVarANN-1];
    outputStandardizeANN_Vec[numVarANN-1] = outputStandardizeANN_Vec[numVarANN-1] + child_MeanOUTPUTLoad[numVarANN-1];
    // Verify if output (Composition space) out-of-bound
    outputStandardizeANN_Vec[numVarANN-1] = std::max(outputStandardizeANN_Vec[numVarANN-1],child_MinOUTPUTLoad[numVarANN-1]);
    outputStandardizeANN_Vec[numVarANN-1] = std::min(outputStandardizeANN_Vec[numVarANN-1],child_MaxOUTPUTLoad[numVarANN-1]);
    // Calculate Temperature and bound on Max Min GLOBAL data
    Tm = Tm + outputStandardizeANN_Vec[numVarANN-1];
    Tm = std::min(std::max(Tm,minAftREACTORLoad[numVarANN-1]),maxAftREACTORLoad[numVarANN-1]); // Check if <minValue or >maxValue on GLOBAL
    
    //  Y
    for (int kANN=0; kANN<nsp-1; kANN++) // Don't take N2 (last species)
    {
        outputStandardizeANN_Vec[kANN] = outputStandardizeANN_Vec[kANN] * child_StdOUTPUTLoad[kANN];
        outputStandardizeANN_Vec[kANN] = outputStandardizeANN_Vec[kANN] + child_MeanOUTPUTLoad[kANN];
        // Verify if output (Composition space) out-of-bound
        outputStandardizeANN_Vec[kANN] = std::max(outputStandardizeANN_Vec[kANN],child_MinOUTPUTLoad[kANN]);
        outputStandardizeANN_Vec[kANN] = std::min(outputStandardizeANN_Vec[kANN],child_MaxOUTPUTLoad[kANN]);
        // Calculate Y and bound on Max Min GLOBAL data
        Ym[kANN] = Ym[kANN] + outputStandardizeANN_Vec[kANN];
        Ym[kANN] = std::min(std::max(Ym[kANN],minAftREACTORLoad[kANN]),maxAftREACTORLoad[kANN]); // Check if <minValue or >maxValue
        //Cross-check if there is N2 inside, if yes, should move N2 to the end of the scheme
        if(mixture->speciesName(kANN) == "N2")
            cout << "Warning : N2 in the ANN!!!" << endl;
    }
}
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
   int rank, nb_processus;

   	if (step != "Optimisation") {
                MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	}

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


    
    int numVarANN = 9; // CasKaidi: 10 species (remove N2),  T (remove H)
    int numComponentPCA = 9;   // Global number of components of each cluster
    if(rank==0)
    {
        cout << "Number of variables: " << numVarANN;
        cout << "| Number of components: " << numComponentPCA;
        cout << endl;
    }
    
    // Overall case path
    std::string casePath = "/home/huutri/workdir/ML-DL_ORCh/dataORChCases/Hydrogen/cas1/Run9-casKissel-exact/";

    if(rank==0)
    {
      cout << " ====== WITHOUT PCA with Heat loss ====== "  << endl;
      cout << "ANN case path: " <<  casePath << endl;
    }
    
    //Declare GLOBAL Standardize for K-means (or ANN Classifier): standardizedX = (X - meanScale) / standardDeviationScale
    std::string GLOBAL_meanScalePATH = casePath + "0.treatedFile/scalerGlobal_OnlyMean.txt";
    vector<float> GLOBAL_meanScale;
    readFromCommaDelimitedFile_Float(GLOBAL_meanScalePATH, GLOBAL_meanScale); // readFromCommaDelimitedFile_Float(path, vectorFloat)
    
    std::string GLOBAL_stdScalePATH = casePath + "0.treatedFile/scalerGlobal_OnlyStandardDeviation.txt";
    vector<float> GLOBAL_stdScale;
    readFromCommaDelimitedFile_Float(GLOBAL_stdScalePATH, GLOBAL_stdScale);
    
    
    // Declare Max Min Label after REACTOR (maxAftREACTOR Y,T not Var) of each child/grandChild cluster (who have ANNReg)
    // Child of parent 0
    std::string child0_0_maxValPATH = casePath + "0.treatedFile/parentCluster0/maxMin_LabelClustered/label_OnlyMax_child0_0.txt";
    std::string child0_0_minValPATH = casePath + "0.treatedFile/parentCluster0/maxMin_LabelClustered/label_OnlyMin_child0_0.txt";
    vector<double> child0_0_maxVal;
    vector<double> child0_0_minVal;
    readFromCommaDelimitedFile_Double(child0_0_maxValPATH, child0_0_maxVal);
    readFromCommaDelimitedFile_Double(child0_0_minValPATH, child0_0_minVal);
    
    std::string child0_1_maxValPATH = casePath + "0.treatedFile/parentCluster0/maxMin_LabelClustered/label_OnlyMax_child0_1.txt";
    std::string child0_1_minValPATH = casePath + "0.treatedFile/parentCluster0/maxMin_LabelClustered/label_OnlyMin_child0_1.txt";
    vector<double> child0_1_maxVal;
    vector<double> child0_1_minVal;
    readFromCommaDelimitedFile_Double(child0_1_maxValPATH, child0_1_maxVal);
    readFromCommaDelimitedFile_Double(child0_1_minValPATH, child0_1_minVal);
    
    std::string child0_2_maxValPATH = casePath + "0.treatedFile/parentCluster0/maxMin_LabelClustered/label_OnlyMax_child0_2.txt";
    std::string child0_2_minValPATH = casePath + "0.treatedFile/parentCluster0/maxMin_LabelClustered/label_OnlyMin_child0_2.txt";
    vector<double> child0_2_maxVal;
    vector<double> child0_2_minVal;
    readFromCommaDelimitedFile_Double(child0_2_maxValPATH, child0_2_maxVal);
    readFromCommaDelimitedFile_Double(child0_2_minValPATH, child0_2_minVal);
    
    std::string child0_3_maxValPATH = casePath + "0.treatedFile/parentCluster0/maxMin_LabelClustered/label_OnlyMax_child0_3.txt";
    std::string child0_3_minValPATH = casePath + "0.treatedFile/parentCluster0/maxMin_LabelClustered/label_OnlyMin_child0_3.txt";
    vector<double> child0_3_maxVal;
    vector<double> child0_3_minVal;
    readFromCommaDelimitedFile_Double(child0_3_maxValPATH, child0_3_maxVal);
    readFromCommaDelimitedFile_Double(child0_3_minValPATH, child0_3_minVal);
    
    std::string child0_4_maxValPATH = casePath + "0.treatedFile/parentCluster0/maxMin_LabelClustered/label_OnlyMax_child0_4.txt";
    std::string child0_4_minValPATH = casePath + "0.treatedFile/parentCluster0/maxMin_LabelClustered/label_OnlyMin_child0_4.txt";
    vector<double> child0_4_maxVal;
    vector<double> child0_4_minVal;
    readFromCommaDelimitedFile_Double(child0_4_maxValPATH, child0_4_maxVal);
    readFromCommaDelimitedFile_Double(child0_4_minValPATH, child0_4_minVal);
    
    std::string child0_5_maxValPATH = casePath + "0.treatedFile/parentCluster0/maxMin_LabelClustered/label_OnlyMax_child0_5.txt";
    std::string child0_5_minValPATH = casePath + "0.treatedFile/parentCluster0/maxMin_LabelClustered/label_OnlyMin_child0_5.txt";
    vector<double> child0_5_maxVal;
    vector<double> child0_5_minVal;
    readFromCommaDelimitedFile_Double(child0_5_maxValPATH, child0_5_maxVal);
    readFromCommaDelimitedFile_Double(child0_5_minValPATH, child0_5_minVal);
    
    std::string child0_6_maxValPATH = casePath + "0.treatedFile/parentCluster0/maxMin_LabelClustered/label_OnlyMax_child0_6.txt";
    std::string child0_6_minValPATH = casePath + "0.treatedFile/parentCluster0/maxMin_LabelClustered/label_OnlyMin_child0_6.txt";
    vector<double> child0_6_maxVal;
    vector<double> child0_6_minVal;
    readFromCommaDelimitedFile_Double(child0_6_maxValPATH, child0_6_maxVal);
    readFromCommaDelimitedFile_Double(child0_6_minValPATH, child0_6_minVal);
    
    std::string child0_7_maxValPATH = casePath + "0.treatedFile/parentCluster0/maxMin_LabelClustered/label_OnlyMax_child0_7.txt";
    std::string child0_7_minValPATH = casePath + "0.treatedFile/parentCluster0/maxMin_LabelClustered/label_OnlyMin_child0_7.txt";
    vector<double> child0_7_maxVal;
    vector<double> child0_7_minVal;
    readFromCommaDelimitedFile_Double(child0_7_maxValPATH, child0_7_maxVal);
    readFromCommaDelimitedFile_Double(child0_7_minValPATH, child0_7_minVal);
    
    std::string child0_8_maxValPATH = casePath + "0.treatedFile/parentCluster0/maxMin_LabelClustered/label_OnlyMax_child0_8.txt";
    std::string child0_8_minValPATH = casePath + "0.treatedFile/parentCluster0/maxMin_LabelClustered/label_OnlyMin_child0_8.txt";
    vector<double> child0_8_maxVal;
    vector<double> child0_8_minVal;
    readFromCommaDelimitedFile_Double(child0_8_maxValPATH, child0_8_maxVal);
    readFromCommaDelimitedFile_Double(child0_8_minValPATH, child0_8_minVal);
    
    std::string child0_9_maxValPATH = casePath + "0.treatedFile/parentCluster0/maxMin_LabelClustered/label_OnlyMax_child0_9.txt";
    std::string child0_9_minValPATH = casePath + "0.treatedFile/parentCluster0/maxMin_LabelClustered/label_OnlyMin_child0_9.txt";
    vector<double> child0_9_maxVal;
    vector<double> child0_9_minVal;
    readFromCommaDelimitedFile_Double(child0_9_maxValPATH, child0_9_maxVal);
    readFromCommaDelimitedFile_Double(child0_9_minValPATH, child0_9_minVal);
    
    std::string child0_10_maxValPATH = casePath + "0.treatedFile/parentCluster0/maxMin_LabelClustered/label_OnlyMax_child0_10.txt";
    std::string child0_10_minValPATH = casePath + "0.treatedFile/parentCluster0/maxMin_LabelClustered/label_OnlyMin_child0_10.txt";
    vector<double> child0_10_maxVal;
    vector<double> child0_10_minVal;
    readFromCommaDelimitedFile_Double(child0_10_maxValPATH, child0_10_maxVal);
    readFromCommaDelimitedFile_Double(child0_10_minValPATH, child0_10_minVal);
    
    std::string child0_11_maxValPATH = casePath + "0.treatedFile/parentCluster0/maxMin_LabelClustered/label_OnlyMax_child0_11.txt";
    std::string child0_11_minValPATH = casePath + "0.treatedFile/parentCluster0/maxMin_LabelClustered/label_OnlyMin_child0_11.txt";
    vector<double> child0_11_maxVal;
    vector<double> child0_11_minVal;
    readFromCommaDelimitedFile_Double(child0_11_maxValPATH, child0_11_maxVal);
    readFromCommaDelimitedFile_Double(child0_11_minValPATH, child0_11_minVal);
    
    std::string child0_12_maxValPATH = casePath + "0.treatedFile/parentCluster0/maxMin_LabelClustered/label_OnlyMax_child0_12.txt";
    std::string child0_12_minValPATH = casePath + "0.treatedFile/parentCluster0/maxMin_LabelClustered/label_OnlyMin_child0_12.txt";
    vector<double> child0_12_maxVal;
    vector<double> child0_12_minVal;
    readFromCommaDelimitedFile_Double(child0_12_maxValPATH, child0_12_maxVal);
    readFromCommaDelimitedFile_Double(child0_12_minValPATH, child0_12_minVal);
    
    std::string child0_13_maxValPATH = casePath + "0.treatedFile/parentCluster0/maxMin_LabelClustered/label_OnlyMax_child0_13.txt";
    std::string child0_13_minValPATH = casePath + "0.treatedFile/parentCluster0/maxMin_LabelClustered/label_OnlyMin_child0_13.txt";
    vector<double> child0_13_maxVal;
    vector<double> child0_13_minVal;
    readFromCommaDelimitedFile_Double(child0_13_maxValPATH, child0_13_maxVal);
    readFromCommaDelimitedFile_Double(child0_13_minValPATH, child0_13_minVal);
    
    std::string child0_14_maxValPATH = casePath + "0.treatedFile/parentCluster0/maxMin_LabelClustered/label_OnlyMax_child0_14.txt";
    std::string child0_14_minValPATH = casePath + "0.treatedFile/parentCluster0/maxMin_LabelClustered/label_OnlyMin_child0_14.txt";
    vector<double> child0_14_maxVal;
    vector<double> child0_14_minVal;
    readFromCommaDelimitedFile_Double(child0_14_maxValPATH, child0_14_maxVal);
    readFromCommaDelimitedFile_Double(child0_14_minValPATH, child0_14_minVal);
    
    std::string child0_15_maxValPATH = casePath + "0.treatedFile/parentCluster0/maxMin_LabelClustered/label_OnlyMax_child0_15.txt";
    std::string child0_15_minValPATH = casePath + "0.treatedFile/parentCluster0/maxMin_LabelClustered/label_OnlyMin_child0_15.txt";
    vector<double> child0_15_maxVal;
    vector<double> child0_15_minVal;
    readFromCommaDelimitedFile_Double(child0_15_maxValPATH, child0_15_maxVal);
    readFromCommaDelimitedFile_Double(child0_15_minValPATH, child0_15_minVal);
    
    std::string child0_16_maxValPATH = casePath + "0.treatedFile/parentCluster0/maxMin_LabelClustered/label_OnlyMax_child0_16.txt";
    std::string child0_16_minValPATH = casePath + "0.treatedFile/parentCluster0/maxMin_LabelClustered/label_OnlyMin_child0_16.txt";
    vector<double> child0_16_maxVal;
    vector<double> child0_16_minVal;
    readFromCommaDelimitedFile_Double(child0_16_maxValPATH, child0_16_maxVal);
    readFromCommaDelimitedFile_Double(child0_16_minValPATH, child0_16_minVal);
    
    std::string child0_17_maxValPATH = casePath + "0.treatedFile/parentCluster0/maxMin_LabelClustered/label_OnlyMax_child0_17.txt";
    std::string child0_17_minValPATH = casePath + "0.treatedFile/parentCluster0/maxMin_LabelClustered/label_OnlyMin_child0_17.txt";
    vector<double> child0_17_maxVal;
    vector<double> child0_17_minVal;
    readFromCommaDelimitedFile_Double(child0_17_maxValPATH, child0_17_maxVal);
    readFromCommaDelimitedFile_Double(child0_17_minValPATH, child0_17_minVal);
    
    std::string child0_18_maxValPATH = casePath + "0.treatedFile/parentCluster0/maxMin_LabelClustered/label_OnlyMax_child0_18.txt";
    std::string child0_18_minValPATH = casePath + "0.treatedFile/parentCluster0/maxMin_LabelClustered/label_OnlyMin_child0_18.txt";
    vector<double> child0_18_maxVal;
    vector<double> child0_18_minVal;
    readFromCommaDelimitedFile_Double(child0_18_maxValPATH, child0_18_maxVal);
    readFromCommaDelimitedFile_Double(child0_18_minValPATH, child0_18_minVal);
    
    std::string child0_19_maxValPATH = casePath + "0.treatedFile/parentCluster0/maxMin_LabelClustered/label_OnlyMax_child0_19.txt";
    std::string child0_19_minValPATH = casePath + "0.treatedFile/parentCluster0/maxMin_LabelClustered/label_OnlyMin_child0_19.txt";
    vector<double> child0_19_maxVal;
    vector<double> child0_19_minVal;
    readFromCommaDelimitedFile_Double(child0_19_maxValPATH, child0_19_maxVal);
    readFromCommaDelimitedFile_Double(child0_19_minValPATH, child0_19_minVal);

    // grandChild of parent 1
    // grandChild of child1_0
    std::string grandChild1_0_0_maxValPATH = casePath + "0.treatedFile/parentCluster1/childCluster0/maxMin_LabelClustered/label_OnlyMax_grandChild1_0_0.txt";
    std::string grandChild1_0_0_minValPATH = casePath + "0.treatedFile/parentCluster1/childCluster0/maxMin_LabelClustered/label_OnlyMin_grandChild1_0_0.txt";
    vector<double> grandChild1_0_0_maxVal;
    vector<double> grandChild1_0_0_minVal;
    readFromCommaDelimitedFile_Double(grandChild1_0_0_maxValPATH, grandChild1_0_0_maxVal);
    readFromCommaDelimitedFile_Double(grandChild1_0_0_minValPATH, grandChild1_0_0_minVal);
    
    std::string grandChild1_0_1_maxValPATH = casePath + "0.treatedFile/parentCluster1/childCluster0/maxMin_LabelClustered/label_OnlyMax_grandChild1_0_1.txt";
    std::string grandChild1_0_1_minValPATH = casePath + "0.treatedFile/parentCluster1/childCluster0/maxMin_LabelClustered/label_OnlyMin_grandChild1_0_1.txt";
    vector<double> grandChild1_0_1_maxVal;
    vector<double> grandChild1_0_1_minVal;
    readFromCommaDelimitedFile_Double(grandChild1_0_1_maxValPATH, grandChild1_0_1_maxVal);
    readFromCommaDelimitedFile_Double(grandChild1_0_1_minValPATH, grandChild1_0_1_minVal);
    
    std::string grandChild1_0_2_maxValPATH = casePath + "0.treatedFile/parentCluster1/childCluster0/maxMin_LabelClustered/label_OnlyMax_grandChild1_0_2.txt";
    std::string grandChild1_0_2_minValPATH = casePath + "0.treatedFile/parentCluster1/childCluster0/maxMin_LabelClustered/label_OnlyMin_grandChild1_0_2.txt";
    vector<double> grandChild1_0_2_maxVal;
    vector<double> grandChild1_0_2_minVal;
    readFromCommaDelimitedFile_Double(grandChild1_0_2_maxValPATH, grandChild1_0_2_maxVal);
    readFromCommaDelimitedFile_Double(grandChild1_0_2_minValPATH, grandChild1_0_2_minVal);
    
    std::string grandChild1_0_3_maxValPATH = casePath + "0.treatedFile/parentCluster1/childCluster0/maxMin_LabelClustered/label_OnlyMax_grandChild1_0_3.txt";
    std::string grandChild1_0_3_minValPATH = casePath + "0.treatedFile/parentCluster1/childCluster0/maxMin_LabelClustered/label_OnlyMin_grandChild1_0_3.txt";
    vector<double> grandChild1_0_3_maxVal;
    vector<double> grandChild1_0_3_minVal;
    readFromCommaDelimitedFile_Double(grandChild1_0_3_maxValPATH, grandChild1_0_3_maxVal);
    readFromCommaDelimitedFile_Double(grandChild1_0_3_minValPATH, grandChild1_0_3_minVal);
    
    std::string grandChild1_0_4_maxValPATH = casePath + "0.treatedFile/parentCluster1/childCluster0/maxMin_LabelClustered/label_OnlyMax_grandChild1_0_4.txt";
    std::string grandChild1_0_4_minValPATH = casePath + "0.treatedFile/parentCluster1/childCluster0/maxMin_LabelClustered/label_OnlyMin_grandChild1_0_4.txt";
    vector<double> grandChild1_0_4_maxVal;
    vector<double> grandChild1_0_4_minVal;
    readFromCommaDelimitedFile_Double(grandChild1_0_4_maxValPATH, grandChild1_0_4_maxVal);
    readFromCommaDelimitedFile_Double(grandChild1_0_4_minValPATH, grandChild1_0_4_minVal);
    
    std::string grandChild1_0_5_maxValPATH = casePath + "0.treatedFile/parentCluster1/childCluster0/maxMin_LabelClustered/label_OnlyMax_grandChild1_0_5.txt";
    std::string grandChild1_0_5_minValPATH = casePath + "0.treatedFile/parentCluster1/childCluster0/maxMin_LabelClustered/label_OnlyMin_grandChild1_0_5.txt";
    vector<double> grandChild1_0_5_maxVal;
    vector<double> grandChild1_0_5_minVal;
    readFromCommaDelimitedFile_Double(grandChild1_0_5_maxValPATH, grandChild1_0_5_maxVal);
    readFromCommaDelimitedFile_Double(grandChild1_0_5_minValPATH, grandChild1_0_5_minVal);
    
    std::string grandChild1_0_6_maxValPATH = casePath + "0.treatedFile/parentCluster1/childCluster0/maxMin_LabelClustered/label_OnlyMax_grandChild1_0_6.txt";
    std::string grandChild1_0_6_minValPATH = casePath + "0.treatedFile/parentCluster1/childCluster0/maxMin_LabelClustered/label_OnlyMin_grandChild1_0_6.txt";
    vector<double> grandChild1_0_6_maxVal;
    vector<double> grandChild1_0_6_minVal;
    readFromCommaDelimitedFile_Double(grandChild1_0_6_maxValPATH, grandChild1_0_6_maxVal);
    readFromCommaDelimitedFile_Double(grandChild1_0_6_minValPATH, grandChild1_0_6_minVal);
    
    std::string grandChild1_0_7_maxValPATH = casePath + "0.treatedFile/parentCluster1/childCluster0/maxMin_LabelClustered/label_OnlyMax_grandChild1_0_7.txt";
    std::string grandChild1_0_7_minValPATH = casePath + "0.treatedFile/parentCluster1/childCluster0/maxMin_LabelClustered/label_OnlyMin_grandChild1_0_7.txt";
    vector<double> grandChild1_0_7_maxVal;
    vector<double> grandChild1_0_7_minVal;
    readFromCommaDelimitedFile_Double(grandChild1_0_7_maxValPATH, grandChild1_0_7_maxVal);
    readFromCommaDelimitedFile_Double(grandChild1_0_7_minValPATH, grandChild1_0_7_minVal);
    
    std::string grandChild1_0_8_maxValPATH = casePath + "0.treatedFile/parentCluster1/childCluster0/maxMin_LabelClustered/label_OnlyMax_grandChild1_0_8.txt";
    std::string grandChild1_0_8_minValPATH = casePath + "0.treatedFile/parentCluster1/childCluster0/maxMin_LabelClustered/label_OnlyMin_grandChild1_0_8.txt";
    vector<double> grandChild1_0_8_maxVal;
    vector<double> grandChild1_0_8_minVal;
    readFromCommaDelimitedFile_Double(grandChild1_0_8_maxValPATH, grandChild1_0_8_maxVal);
    readFromCommaDelimitedFile_Double(grandChild1_0_8_minValPATH, grandChild1_0_8_minVal);
    
    std::string grandChild1_0_9_maxValPATH = casePath + "0.treatedFile/parentCluster1/childCluster0/maxMin_LabelClustered/label_OnlyMax_grandChild1_0_9.txt";
    std::string grandChild1_0_9_minValPATH = casePath + "0.treatedFile/parentCluster1/childCluster0/maxMin_LabelClustered/label_OnlyMin_grandChild1_0_9.txt";
    vector<double> grandChild1_0_9_maxVal;
    vector<double> grandChild1_0_9_minVal;
    readFromCommaDelimitedFile_Double(grandChild1_0_9_maxValPATH, grandChild1_0_9_maxVal);
    readFromCommaDelimitedFile_Double(grandChild1_0_9_minValPATH, grandChild1_0_9_minVal);
    
    std::string grandChild1_0_10_maxValPATH = casePath + "0.treatedFile/parentCluster1/childCluster0/maxMin_LabelClustered/label_OnlyMax_grandChild1_0_10.txt";
    std::string grandChild1_0_10_minValPATH = casePath + "0.treatedFile/parentCluster1/childCluster0/maxMin_LabelClustered/label_OnlyMin_grandChild1_0_10.txt";
    vector<double> grandChild1_0_10_maxVal;
    vector<double> grandChild1_0_10_minVal;
    readFromCommaDelimitedFile_Double(grandChild1_0_10_maxValPATH, grandChild1_0_10_maxVal);
    readFromCommaDelimitedFile_Double(grandChild1_0_10_minValPATH, grandChild1_0_10_minVal);
    
    std::string grandChild1_0_11_maxValPATH = casePath + "0.treatedFile/parentCluster1/childCluster0/maxMin_LabelClustered/label_OnlyMax_grandChild1_0_11.txt";
    std::string grandChild1_0_11_minValPATH = casePath + "0.treatedFile/parentCluster1/childCluster0/maxMin_LabelClustered/label_OnlyMin_grandChild1_0_11.txt";
    vector<double> grandChild1_0_11_maxVal;
    vector<double> grandChild1_0_11_minVal;
    readFromCommaDelimitedFile_Double(grandChild1_0_11_maxValPATH, grandChild1_0_11_maxVal);
    readFromCommaDelimitedFile_Double(grandChild1_0_11_minValPATH, grandChild1_0_11_minVal);
    
    std::string grandChild1_0_12_maxValPATH = casePath + "0.treatedFile/parentCluster1/childCluster0/maxMin_LabelClustered/label_OnlyMax_grandChild1_0_12.txt";
    std::string grandChild1_0_12_minValPATH = casePath + "0.treatedFile/parentCluster1/childCluster0/maxMin_LabelClustered/label_OnlyMin_grandChild1_0_12.txt";
    vector<double> grandChild1_0_12_maxVal;
    vector<double> grandChild1_0_12_minVal;
    readFromCommaDelimitedFile_Double(grandChild1_0_12_maxValPATH, grandChild1_0_12_maxVal);
    readFromCommaDelimitedFile_Double(grandChild1_0_12_minValPATH, grandChild1_0_12_minVal);
    
    std::string grandChild1_0_13_maxValPATH = casePath + "0.treatedFile/parentCluster1/childCluster0/maxMin_LabelClustered/label_OnlyMax_grandChild1_0_13.txt";
    std::string grandChild1_0_13_minValPATH = casePath + "0.treatedFile/parentCluster1/childCluster0/maxMin_LabelClustered/label_OnlyMin_grandChild1_0_13.txt";
    vector<double> grandChild1_0_13_maxVal;
    vector<double> grandChild1_0_13_minVal;
    readFromCommaDelimitedFile_Double(grandChild1_0_13_maxValPATH, grandChild1_0_13_maxVal);
    readFromCommaDelimitedFile_Double(grandChild1_0_13_minValPATH, grandChild1_0_13_minVal);
    
    std::string grandChild1_0_14_maxValPATH = casePath + "0.treatedFile/parentCluster1/childCluster0/maxMin_LabelClustered/label_OnlyMax_grandChild1_0_14.txt";
    std::string grandChild1_0_14_minValPATH = casePath + "0.treatedFile/parentCluster1/childCluster0/maxMin_LabelClustered/label_OnlyMin_grandChild1_0_14.txt";
    vector<double> grandChild1_0_14_maxVal;
    vector<double> grandChild1_0_14_minVal;
    readFromCommaDelimitedFile_Double(grandChild1_0_14_maxValPATH, grandChild1_0_14_maxVal);
    readFromCommaDelimitedFile_Double(grandChild1_0_14_minValPATH, grandChild1_0_14_minVal);
    
    std::string grandChild1_0_15_maxValPATH = casePath + "0.treatedFile/parentCluster1/childCluster0/maxMin_LabelClustered/label_OnlyMax_grandChild1_0_15.txt";
    std::string grandChild1_0_15_minValPATH = casePath + "0.treatedFile/parentCluster1/childCluster0/maxMin_LabelClustered/label_OnlyMin_grandChild1_0_15.txt";
    vector<double> grandChild1_0_15_maxVal;
    vector<double> grandChild1_0_15_minVal;
    readFromCommaDelimitedFile_Double(grandChild1_0_15_maxValPATH, grandChild1_0_15_maxVal);
    readFromCommaDelimitedFile_Double(grandChild1_0_15_minValPATH, grandChild1_0_15_minVal);
    
    std::string grandChild1_0_16_maxValPATH = casePath + "0.treatedFile/parentCluster1/childCluster0/maxMin_LabelClustered/label_OnlyMax_grandChild1_0_16.txt";
    std::string grandChild1_0_16_minValPATH = casePath + "0.treatedFile/parentCluster1/childCluster0/maxMin_LabelClustered/label_OnlyMin_grandChild1_0_16.txt";
    vector<double> grandChild1_0_16_maxVal;
    vector<double> grandChild1_0_16_minVal;
    readFromCommaDelimitedFile_Double(grandChild1_0_16_maxValPATH, grandChild1_0_16_maxVal);
    readFromCommaDelimitedFile_Double(grandChild1_0_16_minValPATH, grandChild1_0_16_minVal);
    
    std::string grandChild1_0_17_maxValPATH = casePath + "0.treatedFile/parentCluster1/childCluster0/maxMin_LabelClustered/label_OnlyMax_grandChild1_0_17.txt";
    std::string grandChild1_0_17_minValPATH = casePath + "0.treatedFile/parentCluster1/childCluster0/maxMin_LabelClustered/label_OnlyMin_grandChild1_0_17.txt";
    vector<double> grandChild1_0_17_maxVal;
    vector<double> grandChild1_0_17_minVal;
    readFromCommaDelimitedFile_Double(grandChild1_0_17_maxValPATH, grandChild1_0_17_maxVal);
    readFromCommaDelimitedFile_Double(grandChild1_0_17_minValPATH, grandChild1_0_17_minVal);
    
    std::string grandChild1_0_18_maxValPATH = casePath + "0.treatedFile/parentCluster1/childCluster0/maxMin_LabelClustered/label_OnlyMax_grandChild1_0_18.txt";
    std::string grandChild1_0_18_minValPATH = casePath + "0.treatedFile/parentCluster1/childCluster0/maxMin_LabelClustered/label_OnlyMin_grandChild1_0_18.txt";
    vector<double> grandChild1_0_18_maxVal;
    vector<double> grandChild1_0_18_minVal;
    readFromCommaDelimitedFile_Double(grandChild1_0_18_maxValPATH, grandChild1_0_18_maxVal);
    readFromCommaDelimitedFile_Double(grandChild1_0_18_minValPATH, grandChild1_0_18_minVal);
    
    std::string grandChild1_0_19_maxValPATH = casePath + "0.treatedFile/parentCluster1/childCluster0/maxMin_LabelClustered/label_OnlyMax_grandChild1_0_19.txt";
    std::string grandChild1_0_19_minValPATH = casePath + "0.treatedFile/parentCluster1/childCluster0/maxMin_LabelClustered/label_OnlyMin_grandChild1_0_19.txt";
    vector<double> grandChild1_0_19_maxVal;
    vector<double> grandChild1_0_19_minVal;
    readFromCommaDelimitedFile_Double(grandChild1_0_19_maxValPATH, grandChild1_0_19_maxVal);
    readFromCommaDelimitedFile_Double(grandChild1_0_19_minValPATH, grandChild1_0_19_minVal);
    
    
    // grandChild of child1_1
    std::string grandChild1_1_0_maxValPATH = casePath + "0.treatedFile/parentCluster1/childCluster1/maxMin_LabelClustered/label_OnlyMax_grandChild1_1_0.txt";
    std::string grandChild1_1_0_minValPATH = casePath + "0.treatedFile/parentCluster1/childCluster1/maxMin_LabelClustered/label_OnlyMin_grandChild1_1_0.txt";
    vector<double> grandChild1_1_0_maxVal;
    vector<double> grandChild1_1_0_minVal;
    readFromCommaDelimitedFile_Double(grandChild1_1_0_maxValPATH, grandChild1_1_0_maxVal);
    readFromCommaDelimitedFile_Double(grandChild1_1_0_minValPATH, grandChild1_1_0_minVal);
    
    std::string grandChild1_1_1_maxValPATH = casePath + "0.treatedFile/parentCluster1/childCluster1/maxMin_LabelClustered/label_OnlyMax_grandChild1_1_1.txt";
    std::string grandChild1_1_1_minValPATH = casePath + "0.treatedFile/parentCluster1/childCluster1/maxMin_LabelClustered/label_OnlyMin_grandChild1_1_1.txt";
    vector<double> grandChild1_1_1_maxVal;
    vector<double> grandChild1_1_1_minVal;
    readFromCommaDelimitedFile_Double(grandChild1_1_1_maxValPATH, grandChild1_1_1_maxVal);
    readFromCommaDelimitedFile_Double(grandChild1_1_1_minValPATH, grandChild1_1_1_minVal);
    
    std::string grandChild1_1_2_maxValPATH = casePath + "0.treatedFile/parentCluster1/childCluster1/maxMin_LabelClustered/label_OnlyMax_grandChild1_1_2.txt";
    std::string grandChild1_1_2_minValPATH = casePath + "0.treatedFile/parentCluster1/childCluster1/maxMin_LabelClustered/label_OnlyMin_grandChild1_1_2.txt";
    vector<double> grandChild1_1_2_maxVal;
    vector<double> grandChild1_1_2_minVal;
    readFromCommaDelimitedFile_Double(grandChild1_1_2_maxValPATH, grandChild1_1_2_maxVal);
    readFromCommaDelimitedFile_Double(grandChild1_1_2_minValPATH, grandChild1_1_2_minVal);
    
    std::string grandChild1_1_3_maxValPATH = casePath + "0.treatedFile/parentCluster1/childCluster1/maxMin_LabelClustered/label_OnlyMax_grandChild1_1_3.txt";
    std::string grandChild1_1_3_minValPATH = casePath + "0.treatedFile/parentCluster1/childCluster1/maxMin_LabelClustered/label_OnlyMin_grandChild1_1_3.txt";
    vector<double> grandChild1_1_3_maxVal;
    vector<double> grandChild1_1_3_minVal;
    readFromCommaDelimitedFile_Double(grandChild1_1_3_maxValPATH, grandChild1_1_3_maxVal);
    readFromCommaDelimitedFile_Double(grandChild1_1_3_minValPATH, grandChild1_1_3_minVal);
    
    std::string grandChild1_1_4_maxValPATH = casePath + "0.treatedFile/parentCluster1/childCluster1/maxMin_LabelClustered/label_OnlyMax_grandChild1_1_4.txt";
    std::string grandChild1_1_4_minValPATH = casePath + "0.treatedFile/parentCluster1/childCluster1/maxMin_LabelClustered/label_OnlyMin_grandChild1_1_4.txt";
    vector<double> grandChild1_1_4_maxVal;
    vector<double> grandChild1_1_4_minVal;
    readFromCommaDelimitedFile_Double(grandChild1_1_4_maxValPATH, grandChild1_1_4_maxVal);
    readFromCommaDelimitedFile_Double(grandChild1_1_4_minValPATH, grandChild1_1_4_minVal);
    
    std::string grandChild1_1_5_maxValPATH = casePath + "0.treatedFile/parentCluster1/childCluster1/maxMin_LabelClustered/label_OnlyMax_grandChild1_1_5.txt";
    std::string grandChild1_1_5_minValPATH = casePath + "0.treatedFile/parentCluster1/childCluster1/maxMin_LabelClustered/label_OnlyMin_grandChild1_1_5.txt";
    vector<double> grandChild1_1_5_maxVal;
    vector<double> grandChild1_1_5_minVal;
    readFromCommaDelimitedFile_Double(grandChild1_1_5_maxValPATH, grandChild1_1_5_maxVal);
    readFromCommaDelimitedFile_Double(grandChild1_1_5_minValPATH, grandChild1_1_5_minVal);
    
    std::string grandChild1_1_6_maxValPATH = casePath + "0.treatedFile/parentCluster1/childCluster1/maxMin_LabelClustered/label_OnlyMax_grandChild1_1_6.txt";
    std::string grandChild1_1_6_minValPATH = casePath + "0.treatedFile/parentCluster1/childCluster1/maxMin_LabelClustered/label_OnlyMin_grandChild1_1_6.txt";
    vector<double> grandChild1_1_6_maxVal;
    vector<double> grandChild1_1_6_minVal;
    readFromCommaDelimitedFile_Double(grandChild1_1_6_maxValPATH, grandChild1_1_6_maxVal);
    readFromCommaDelimitedFile_Double(grandChild1_1_6_minValPATH, grandChild1_1_6_minVal);
    
    std::string grandChild1_1_7_maxValPATH = casePath + "0.treatedFile/parentCluster1/childCluster1/maxMin_LabelClustered/label_OnlyMax_grandChild1_1_7.txt";
    std::string grandChild1_1_7_minValPATH = casePath + "0.treatedFile/parentCluster1/childCluster1/maxMin_LabelClustered/label_OnlyMin_grandChild1_1_7.txt";
    vector<double> grandChild1_1_7_maxVal;
    vector<double> grandChild1_1_7_minVal;
    readFromCommaDelimitedFile_Double(grandChild1_1_7_maxValPATH, grandChild1_1_7_maxVal);
    readFromCommaDelimitedFile_Double(grandChild1_1_7_minValPATH, grandChild1_1_7_minVal);
    
    std::string grandChild1_1_8_maxValPATH = casePath + "0.treatedFile/parentCluster1/childCluster1/maxMin_LabelClustered/label_OnlyMax_grandChild1_1_8.txt";
    std::string grandChild1_1_8_minValPATH = casePath + "0.treatedFile/parentCluster1/childCluster1/maxMin_LabelClustered/label_OnlyMin_grandChild1_1_8.txt";
    vector<double> grandChild1_1_8_maxVal;
    vector<double> grandChild1_1_8_minVal;
    readFromCommaDelimitedFile_Double(grandChild1_1_8_maxValPATH, grandChild1_1_8_maxVal);
    readFromCommaDelimitedFile_Double(grandChild1_1_8_minValPATH, grandChild1_1_8_minVal);
    
    std::string grandChild1_1_9_maxValPATH = casePath + "0.treatedFile/parentCluster1/childCluster1/maxMin_LabelClustered/label_OnlyMax_grandChild1_1_9.txt";
    std::string grandChild1_1_9_minValPATH = casePath + "0.treatedFile/parentCluster1/childCluster1/maxMin_LabelClustered/label_OnlyMin_grandChild1_1_9.txt";
    vector<double> grandChild1_1_9_maxVal;
    vector<double> grandChild1_1_9_minVal;
    readFromCommaDelimitedFile_Double(grandChild1_1_9_maxValPATH, grandChild1_1_9_maxVal);
    readFromCommaDelimitedFile_Double(grandChild1_1_9_minValPATH, grandChild1_1_9_minVal);
    
    
    
    // grandChild of child1_2
    std::string grandChild1_2_0_maxValPATH = casePath + "0.treatedFile/parentCluster1/childCluster2/maxMin_LabelClustered/label_OnlyMax_grandChild1_2_0.txt";
    std::string grandChild1_2_0_minValPATH = casePath + "0.treatedFile/parentCluster1/childCluster2/maxMin_LabelClustered/label_OnlyMin_grandChild1_2_0.txt";
    vector<double> grandChild1_2_0_maxVal;
    vector<double> grandChild1_2_0_minVal;
    readFromCommaDelimitedFile_Double(grandChild1_2_0_maxValPATH, grandChild1_2_0_maxVal);
    readFromCommaDelimitedFile_Double(grandChild1_2_0_minValPATH, grandChild1_2_0_minVal);
    
    std::string grandChild1_2_1_maxValPATH = casePath + "0.treatedFile/parentCluster1/childCluster2/maxMin_LabelClustered/label_OnlyMax_grandChild1_2_1.txt";
    std::string grandChild1_2_1_minValPATH = casePath + "0.treatedFile/parentCluster1/childCluster2/maxMin_LabelClustered/label_OnlyMin_grandChild1_2_1.txt";
    vector<double> grandChild1_2_1_maxVal;
    vector<double> grandChild1_2_1_minVal;
    readFromCommaDelimitedFile_Double(grandChild1_2_1_maxValPATH, grandChild1_2_1_maxVal);
    readFromCommaDelimitedFile_Double(grandChild1_2_1_minValPATH, grandChild1_2_1_minVal);
    
    std::string grandChild1_2_2_maxValPATH = casePath + "0.treatedFile/parentCluster1/childCluster2/maxMin_LabelClustered/label_OnlyMax_grandChild1_2_2.txt";
    std::string grandChild1_2_2_minValPATH = casePath + "0.treatedFile/parentCluster1/childCluster2/maxMin_LabelClustered/label_OnlyMin_grandChild1_2_2.txt";
    vector<double> grandChild1_2_2_maxVal;
    vector<double> grandChild1_2_2_minVal;
    readFromCommaDelimitedFile_Double(grandChild1_2_2_maxValPATH, grandChild1_2_2_maxVal);
    readFromCommaDelimitedFile_Double(grandChild1_2_2_minValPATH, grandChild1_2_2_minVal);
    
    std::string grandChild1_2_3_maxValPATH = casePath + "0.treatedFile/parentCluster1/childCluster2/maxMin_LabelClustered/label_OnlyMax_grandChild1_2_3.txt";
    std::string grandChild1_2_3_minValPATH = casePath + "0.treatedFile/parentCluster1/childCluster2/maxMin_LabelClustered/label_OnlyMin_grandChild1_2_3.txt";
    vector<double> grandChild1_2_3_maxVal;
    vector<double> grandChild1_2_3_minVal;
    readFromCommaDelimitedFile_Double(grandChild1_2_3_maxValPATH, grandChild1_2_3_maxVal);
    readFromCommaDelimitedFile_Double(grandChild1_2_3_minValPATH, grandChild1_2_3_minVal);
    
    std::string grandChild1_2_4_maxValPATH = casePath + "0.treatedFile/parentCluster1/childCluster2/maxMin_LabelClustered/label_OnlyMax_grandChild1_2_4.txt";
    std::string grandChild1_2_4_minValPATH = casePath + "0.treatedFile/parentCluster1/childCluster2/maxMin_LabelClustered/label_OnlyMin_grandChild1_2_4.txt";
    vector<double> grandChild1_2_4_maxVal;
    vector<double> grandChild1_2_4_minVal;
    readFromCommaDelimitedFile_Double(grandChild1_2_4_maxValPATH, grandChild1_2_4_maxVal);
    readFromCommaDelimitedFile_Double(grandChild1_2_4_minValPATH, grandChild1_2_4_minVal);
    
    std::string grandChild1_2_5_maxValPATH = casePath + "0.treatedFile/parentCluster1/childCluster2/maxMin_LabelClustered/label_OnlyMax_grandChild1_2_5.txt";
    std::string grandChild1_2_5_minValPATH = casePath + "0.treatedFile/parentCluster1/childCluster2/maxMin_LabelClustered/label_OnlyMin_grandChild1_2_5.txt";
    vector<double> grandChild1_2_5_maxVal;
    vector<double> grandChild1_2_5_minVal;
    readFromCommaDelimitedFile_Double(grandChild1_2_5_maxValPATH, grandChild1_2_5_maxVal);
    readFromCommaDelimitedFile_Double(grandChild1_2_5_minValPATH, grandChild1_2_5_minVal);
    
    std::string grandChild1_2_6_maxValPATH = casePath + "0.treatedFile/parentCluster1/childCluster2/maxMin_LabelClustered/label_OnlyMax_grandChild1_2_6.txt";
    std::string grandChild1_2_6_minValPATH = casePath + "0.treatedFile/parentCluster1/childCluster2/maxMin_LabelClustered/label_OnlyMin_grandChild1_2_6.txt";
    vector<double> grandChild1_2_6_maxVal;
    vector<double> grandChild1_2_6_minVal;
    readFromCommaDelimitedFile_Double(grandChild1_2_6_maxValPATH, grandChild1_2_6_maxVal);
    readFromCommaDelimitedFile_Double(grandChild1_2_6_minValPATH, grandChild1_2_6_minVal);
    
    std::string grandChild1_2_7_maxValPATH = casePath + "0.treatedFile/parentCluster1/childCluster2/maxMin_LabelClustered/label_OnlyMax_grandChild1_2_7.txt";
    std::string grandChild1_2_7_minValPATH = casePath + "0.treatedFile/parentCluster1/childCluster2/maxMin_LabelClustered/label_OnlyMin_grandChild1_2_7.txt";
    vector<double> grandChild1_2_7_maxVal;
    vector<double> grandChild1_2_7_minVal;
    readFromCommaDelimitedFile_Double(grandChild1_2_7_maxValPATH, grandChild1_2_7_maxVal);
    readFromCommaDelimitedFile_Double(grandChild1_2_7_minValPATH, grandChild1_2_7_minVal);
    
    std::string grandChild1_2_8_maxValPATH = casePath + "0.treatedFile/parentCluster1/childCluster2/maxMin_LabelClustered/label_OnlyMax_grandChild1_2_8.txt";
    std::string grandChild1_2_8_minValPATH = casePath + "0.treatedFile/parentCluster1/childCluster2/maxMin_LabelClustered/label_OnlyMin_grandChild1_2_8.txt";
    vector<double> grandChild1_2_8_maxVal;
    vector<double> grandChild1_2_8_minVal;
    readFromCommaDelimitedFile_Double(grandChild1_2_8_maxValPATH, grandChild1_2_8_maxVal);
    readFromCommaDelimitedFile_Double(grandChild1_2_8_minValPATH, grandChild1_2_8_minVal);
    
    std::string grandChild1_2_9_maxValPATH = casePath + "0.treatedFile/parentCluster1/childCluster2/maxMin_LabelClustered/label_OnlyMax_grandChild1_2_9.txt";
    std::string grandChild1_2_9_minValPATH = casePath + "0.treatedFile/parentCluster1/childCluster2/maxMin_LabelClustered/label_OnlyMin_grandChild1_2_9.txt";
    vector<double> grandChild1_2_9_maxVal;
    vector<double> grandChild1_2_9_minVal;
    readFromCommaDelimitedFile_Double(grandChild1_2_9_maxValPATH, grandChild1_2_9_maxVal);
    readFromCommaDelimitedFile_Double(grandChild1_2_9_minValPATH, grandChild1_2_9_minVal);
    
    
    
    //Declare centroids for 2 PARENT clusters to calculate distance >> Classifier indexParentCluster
    std::string centroid_parentCluster0PATH = casePath + "0.treatedFile/centroids_GlobalKmeans_parentCluster0.txt";
    std::string centroid_parentCluster1PATH = casePath + "0.treatedFile/centroids_GlobalKmeans_parentCluster1.txt";
    vector<float> centroid_parentCluster0;
    vector<float> centroid_parentCluster1;
    readFromCommaDelimitedFile_Float(centroid_parentCluster0PATH, centroid_parentCluster0);
    readFromCommaDelimitedFile_Float(centroid_parentCluster1PATH, centroid_parentCluster1);
    // Transform to openCV MAT
    cv::Mat r0 = cv::Mat(numVarANN, 1, CV_32F, centroid_parentCluster0.data());
    cv::Mat r1 = cv::Mat(numVarANN, 1, CV_32F, centroid_parentCluster1.data());
    
    //Declare PARENT LOCAL INPUT Standardize for each PARENT cluster: standardizedX = (X - meanScale) / standardDeviationScale
    std::string parentCluster0_meanScaleINPUTPATH = casePath + "0.treatedFile/parentCluster0/scalerParentCluster0_OnlyMean.txt";
    vector<float> parentCluster0_meanScaleINPUT;
    readFromCommaDelimitedFile_Float(parentCluster0_meanScaleINPUTPATH, parentCluster0_meanScaleINPUT);
    std::string parentCluster0_stdScaleINPUTPATH = casePath + "0.treatedFile/parentCluster0/scalerParentCluster0_OnlyStandardDeviation.txt";
    vector<float> parentCluster0_stdScaleINPUT;
    readFromCommaDelimitedFile_Float(parentCluster0_stdScaleINPUTPATH, parentCluster0_stdScaleINPUT);

    std::string parentCluster1_meanScaleINPUTPATH = casePath + "0.treatedFile/parentCluster1/scalerParentCluster1_OnlyMean.txt";
    vector<float> parentCluster1_meanScaleINPUT;
    readFromCommaDelimitedFile_Float(parentCluster1_meanScaleINPUTPATH, parentCluster1_meanScaleINPUT);
    std::string parentCluster1_stdScaleINPUTPATH = casePath + "0.treatedFile/parentCluster1/scalerParentCluster1_OnlyStandardDeviation.txt";
    vector<float> parentCluster1_stdScaleINPUT;
    readFromCommaDelimitedFile_Float(parentCluster1_stdScaleINPUTPATH, parentCluster1_stdScaleINPUT);
    
    //Declare centroids for CHILD cluster to calculate distance >> Classifier >> indexChildCluster
    // Child of parentCluster0
    std::string centroid_childCluster0_0PATH = casePath + "0.treatedFile/parentCluster0/centroids_parentCluster0_child0.txt";
    vector<float> centroid_childCluster0_0;
    readFromCommaDelimitedFile_Float(centroid_childCluster0_0PATH, centroid_childCluster0_0);

    std::string centroid_childCluster0_1PATH = casePath + "0.treatedFile/parentCluster0/centroids_parentCluster0_child1.txt";
    vector<float> centroid_childCluster0_1;
    readFromCommaDelimitedFile_Float(centroid_childCluster0_1PATH, centroid_childCluster0_1);
    
    std::string centroid_childCluster0_2PATH = casePath + "0.treatedFile/parentCluster0/centroids_parentCluster0_child2.txt";
    vector<float> centroid_childCluster0_2;
    readFromCommaDelimitedFile_Float(centroid_childCluster0_2PATH, centroid_childCluster0_2);
    
    std::string centroid_childCluster0_3PATH = casePath + "0.treatedFile/parentCluster0/centroids_parentCluster0_child3.txt";
    vector<float> centroid_childCluster0_3;
    readFromCommaDelimitedFile_Float(centroid_childCluster0_3PATH, centroid_childCluster0_3);
    
    std::string centroid_childCluster0_4PATH = casePath + "0.treatedFile/parentCluster0/centroids_parentCluster0_child4.txt";
    vector<float> centroid_childCluster0_4;
    readFromCommaDelimitedFile_Float(centroid_childCluster0_4PATH, centroid_childCluster0_4);
    
    std::string centroid_childCluster0_5PATH = casePath + "0.treatedFile/parentCluster0/centroids_parentCluster0_child5.txt";
    vector<float> centroid_childCluster0_5;
    readFromCommaDelimitedFile_Float(centroid_childCluster0_5PATH, centroid_childCluster0_5);
    
    std::string centroid_childCluster0_6PATH = casePath + "0.treatedFile/parentCluster0/centroids_parentCluster0_child6.txt";
    vector<float> centroid_childCluster0_6;
    readFromCommaDelimitedFile_Float(centroid_childCluster0_6PATH, centroid_childCluster0_6);
    
    std::string centroid_childCluster0_7PATH = casePath + "0.treatedFile/parentCluster0/centroids_parentCluster0_child7.txt";
    vector<float> centroid_childCluster0_7;
    readFromCommaDelimitedFile_Float(centroid_childCluster0_7PATH, centroid_childCluster0_7);
    
    std::string centroid_childCluster0_8PATH = casePath + "0.treatedFile/parentCluster0/centroids_parentCluster0_child8.txt";
    vector<float> centroid_childCluster0_8;
    readFromCommaDelimitedFile_Float(centroid_childCluster0_8PATH, centroid_childCluster0_8);
    
    std::string centroid_childCluster0_9PATH = casePath + "0.treatedFile/parentCluster0/centroids_parentCluster0_child9.txt";
    vector<float> centroid_childCluster0_9;
    readFromCommaDelimitedFile_Float(centroid_childCluster0_9PATH, centroid_childCluster0_9);
    
    std::string centroid_childCluster0_10PATH = casePath + "0.treatedFile/parentCluster0/centroids_parentCluster0_child10.txt";
    vector<float> centroid_childCluster0_10;
    readFromCommaDelimitedFile_Float(centroid_childCluster0_10PATH, centroid_childCluster0_10);
    
    std::string centroid_childCluster0_11PATH = casePath + "0.treatedFile/parentCluster0/centroids_parentCluster0_child11.txt";
    vector<float> centroid_childCluster0_11;
    readFromCommaDelimitedFile_Float(centroid_childCluster0_11PATH, centroid_childCluster0_11);
    
    std::string centroid_childCluster0_12PATH = casePath + "0.treatedFile/parentCluster0/centroids_parentCluster0_child12.txt";
    vector<float> centroid_childCluster0_12;
    readFromCommaDelimitedFile_Float(centroid_childCluster0_12PATH, centroid_childCluster0_12);
    
    std::string centroid_childCluster0_13PATH = casePath + "0.treatedFile/parentCluster0/centroids_parentCluster0_child13.txt";
    vector<float> centroid_childCluster0_13;
    readFromCommaDelimitedFile_Float(centroid_childCluster0_13PATH, centroid_childCluster0_13);
    
    std::string centroid_childCluster0_14PATH = casePath + "0.treatedFile/parentCluster0/centroids_parentCluster0_child14.txt";
    vector<float> centroid_childCluster0_14;
    readFromCommaDelimitedFile_Float(centroid_childCluster0_14PATH, centroid_childCluster0_14);
    
    std::string centroid_childCluster0_15PATH = casePath + "0.treatedFile/parentCluster0/centroids_parentCluster0_child15.txt";
    vector<float> centroid_childCluster0_15;
    readFromCommaDelimitedFile_Float(centroid_childCluster0_15PATH, centroid_childCluster0_15);
    
    std::string centroid_childCluster0_16PATH = casePath + "0.treatedFile/parentCluster0/centroids_parentCluster0_child16.txt";
    vector<float> centroid_childCluster0_16;
    readFromCommaDelimitedFile_Float(centroid_childCluster0_16PATH, centroid_childCluster0_16);
    
    std::string centroid_childCluster0_17PATH = casePath + "0.treatedFile/parentCluster0/centroids_parentCluster0_child17.txt";
    vector<float> centroid_childCluster0_17;
    readFromCommaDelimitedFile_Float(centroid_childCluster0_17PATH, centroid_childCluster0_17);
    
    std::string centroid_childCluster0_18PATH = casePath + "0.treatedFile/parentCluster0/centroids_parentCluster0_child18.txt";
    vector<float> centroid_childCluster0_18;
    readFromCommaDelimitedFile_Float(centroid_childCluster0_18PATH, centroid_childCluster0_18);
    
    std::string centroid_childCluster0_19PATH = casePath + "0.treatedFile/parentCluster0/centroids_parentCluster0_child19.txt";
    vector<float> centroid_childCluster0_19;
    readFromCommaDelimitedFile_Float(centroid_childCluster0_19PATH, centroid_childCluster0_19);
    
    // Transform to openCV MAT
    cv::Mat r0_0 = cv::Mat(numVarANN, 1, CV_32F, centroid_childCluster0_0.data());
    cv::Mat r0_1 = cv::Mat(numVarANN, 1, CV_32F, centroid_childCluster0_1.data());
    cv::Mat r0_2 = cv::Mat(numVarANN, 1, CV_32F, centroid_childCluster0_2.data());
    cv::Mat r0_3 = cv::Mat(numVarANN, 1, CV_32F, centroid_childCluster0_3.data());
    cv::Mat r0_4 = cv::Mat(numVarANN, 1, CV_32F, centroid_childCluster0_4.data());
    cv::Mat r0_5 = cv::Mat(numVarANN, 1, CV_32F, centroid_childCluster0_5.data());
    cv::Mat r0_6 = cv::Mat(numVarANN, 1, CV_32F, centroid_childCluster0_6.data());
    cv::Mat r0_7 = cv::Mat(numVarANN, 1, CV_32F, centroid_childCluster0_7.data());
    cv::Mat r0_8 = cv::Mat(numVarANN, 1, CV_32F, centroid_childCluster0_8.data());
    cv::Mat r0_9 = cv::Mat(numVarANN, 1, CV_32F, centroid_childCluster0_9.data());
    cv::Mat r0_10 = cv::Mat(numVarANN, 1, CV_32F, centroid_childCluster0_10.data());
    cv::Mat r0_11 = cv::Mat(numVarANN, 1, CV_32F, centroid_childCluster0_11.data());
    cv::Mat r0_12 = cv::Mat(numVarANN, 1, CV_32F, centroid_childCluster0_12.data());
    cv::Mat r0_13 = cv::Mat(numVarANN, 1, CV_32F, centroid_childCluster0_13.data());
    cv::Mat r0_14 = cv::Mat(numVarANN, 1, CV_32F, centroid_childCluster0_14.data());
    cv::Mat r0_15 = cv::Mat(numVarANN, 1, CV_32F, centroid_childCluster0_15.data());
    cv::Mat r0_16 = cv::Mat(numVarANN, 1, CV_32F, centroid_childCluster0_16.data());
    cv::Mat r0_17 = cv::Mat(numVarANN, 1, CV_32F, centroid_childCluster0_17.data());
    cv::Mat r0_18 = cv::Mat(numVarANN, 1, CV_32F, centroid_childCluster0_18.data());
    cv::Mat r0_19 = cv::Mat(numVarANN, 1, CV_32F, centroid_childCluster0_19.data());
    
    // Child of parentCluster1
    std::string centroid_childCluster1_0PATH = casePath + "0.treatedFile/parentCluster1/centroids_parentCluster1_child0.txt";
    vector<float> centroid_childCluster1_0;
    readFromCommaDelimitedFile_Float(centroid_childCluster1_0PATH, centroid_childCluster1_0);
    
    std::string centroid_childCluster1_1PATH = casePath + "0.treatedFile/parentCluster1/centroids_parentCluster1_child1.txt";
    vector<float> centroid_childCluster1_1;
    readFromCommaDelimitedFile_Float(centroid_childCluster1_1PATH, centroid_childCluster1_1);
    
    std::string centroid_childCluster1_2PATH = casePath + "0.treatedFile/parentCluster1/centroids_parentCluster1_child2.txt";
    vector<float> centroid_childCluster1_2;
    readFromCommaDelimitedFile_Float(centroid_childCluster1_2PATH, centroid_childCluster1_2);
    // Transform to openCV MAT
    cv::Mat r1_0 = cv::Mat(numVarANN, 1, CV_32F, centroid_childCluster1_0.data());
    cv::Mat r1_1 = cv::Mat(numVarANN, 1, CV_32F, centroid_childCluster1_1.data());
    cv::Mat r1_2 = cv::Mat(numVarANN, 1, CV_32F, centroid_childCluster1_2.data());

    
    // parentCluster1 - childCluster0 - grandChild
    std::string centroid_childCluster1_0_0PATH = casePath + "0.treatedFile/parentCluster1/childCluster0/centroids_parentCluster1_child0_grandChild0.txt";
    vector<float> centroid_childCluster1_0_0;
    readFromCommaDelimitedFile_Float(centroid_childCluster1_0_0PATH, centroid_childCluster1_0_0);
    
    std::string centroid_childCluster1_0_1PATH = casePath + "0.treatedFile/parentCluster1/childCluster0/centroids_parentCluster1_child0_grandChild1.txt";
    vector<float> centroid_childCluster1_0_1;
    readFromCommaDelimitedFile_Float(centroid_childCluster1_0_1PATH, centroid_childCluster1_0_1);
    
    std::string centroid_childCluster1_0_2PATH = casePath + "0.treatedFile/parentCluster1/childCluster0/centroids_parentCluster1_child0_grandChild2.txt";
    vector<float> centroid_childCluster1_0_2;
    readFromCommaDelimitedFile_Float(centroid_childCluster1_0_2PATH, centroid_childCluster1_0_2);
    
    std::string centroid_childCluster1_0_3PATH = casePath + "0.treatedFile/parentCluster1/childCluster0/centroids_parentCluster1_child0_grandChild3.txt";
    vector<float> centroid_childCluster1_0_3;
    readFromCommaDelimitedFile_Float(centroid_childCluster1_0_3PATH, centroid_childCluster1_0_3);
    
    std::string centroid_childCluster1_0_4PATH = casePath + "0.treatedFile/parentCluster1/childCluster0/centroids_parentCluster1_child0_grandChild4.txt";
    vector<float> centroid_childCluster1_0_4;
    readFromCommaDelimitedFile_Float(centroid_childCluster1_0_4PATH, centroid_childCluster1_0_4);
    
    std::string centroid_childCluster1_0_5PATH = casePath + "0.treatedFile/parentCluster1/childCluster0/centroids_parentCluster1_child0_grandChild5.txt";
    vector<float> centroid_childCluster1_0_5;
    readFromCommaDelimitedFile_Float(centroid_childCluster1_0_5PATH, centroid_childCluster1_0_5);
    
    std::string centroid_childCluster1_0_6PATH = casePath + "0.treatedFile/parentCluster1/childCluster0/centroids_parentCluster1_child0_grandChild6.txt";
    vector<float> centroid_childCluster1_0_6;
    readFromCommaDelimitedFile_Float(centroid_childCluster1_0_6PATH, centroid_childCluster1_0_6);
    
    std::string centroid_childCluster1_0_7PATH = casePath + "0.treatedFile/parentCluster1/childCluster0/centroids_parentCluster1_child0_grandChild7.txt";
    vector<float> centroid_childCluster1_0_7;
    readFromCommaDelimitedFile_Float(centroid_childCluster1_0_7PATH, centroid_childCluster1_0_7);
    
    std::string centroid_childCluster1_0_8PATH = casePath + "0.treatedFile/parentCluster1/childCluster0/centroids_parentCluster1_child0_grandChild8.txt";
    vector<float> centroid_childCluster1_0_8;
    readFromCommaDelimitedFile_Float(centroid_childCluster1_0_8PATH, centroid_childCluster1_0_8);
    
    std::string centroid_childCluster1_0_9PATH = casePath + "0.treatedFile/parentCluster1/childCluster0/centroids_parentCluster1_child0_grandChild9.txt";
    vector<float> centroid_childCluster1_0_9;
    readFromCommaDelimitedFile_Float(centroid_childCluster1_0_9PATH, centroid_childCluster1_0_9);
    
    std::string centroid_childCluster1_0_10PATH = casePath + "0.treatedFile/parentCluster1/childCluster0/centroids_parentCluster1_child0_grandChild10.txt";
    vector<float> centroid_childCluster1_0_10;
    readFromCommaDelimitedFile_Float(centroid_childCluster1_0_10PATH, centroid_childCluster1_0_10);
    
    std::string centroid_childCluster1_0_11PATH = casePath + "0.treatedFile/parentCluster1/childCluster0/centroids_parentCluster1_child0_grandChild11.txt";
    vector<float> centroid_childCluster1_0_11;
    readFromCommaDelimitedFile_Float(centroid_childCluster1_0_11PATH, centroid_childCluster1_0_11);
    
    std::string centroid_childCluster1_0_12PATH = casePath + "0.treatedFile/parentCluster1/childCluster0/centroids_parentCluster1_child0_grandChild12.txt";
    vector<float> centroid_childCluster1_0_12;
    readFromCommaDelimitedFile_Float(centroid_childCluster1_0_12PATH, centroid_childCluster1_0_12);
    
    std::string centroid_childCluster1_0_13PATH = casePath + "0.treatedFile/parentCluster1/childCluster0/centroids_parentCluster1_child0_grandChild13.txt";
    vector<float> centroid_childCluster1_0_13;
    readFromCommaDelimitedFile_Float(centroid_childCluster1_0_13PATH, centroid_childCluster1_0_13);
    
    std::string centroid_childCluster1_0_14PATH = casePath + "0.treatedFile/parentCluster1/childCluster0/centroids_parentCluster1_child0_grandChild14.txt";
    vector<float> centroid_childCluster1_0_14;
    readFromCommaDelimitedFile_Float(centroid_childCluster1_0_14PATH, centroid_childCluster1_0_14);
    
    std::string centroid_childCluster1_0_15PATH = casePath + "0.treatedFile/parentCluster1/childCluster0/centroids_parentCluster1_child0_grandChild15.txt";
    vector<float> centroid_childCluster1_0_15;
    readFromCommaDelimitedFile_Float(centroid_childCluster1_0_15PATH, centroid_childCluster1_0_15);
    
    std::string centroid_childCluster1_0_16PATH = casePath + "0.treatedFile/parentCluster1/childCluster0/centroids_parentCluster1_child0_grandChild16.txt";
    vector<float> centroid_childCluster1_0_16;
    readFromCommaDelimitedFile_Float(centroid_childCluster1_0_16PATH, centroid_childCluster1_0_16);
    
    std::string centroid_childCluster1_0_17PATH = casePath + "0.treatedFile/parentCluster1/childCluster0/centroids_parentCluster1_child0_grandChild17.txt";
    vector<float> centroid_childCluster1_0_17;
    readFromCommaDelimitedFile_Float(centroid_childCluster1_0_17PATH, centroid_childCluster1_0_17);
    
    std::string centroid_childCluster1_0_18PATH = casePath + "0.treatedFile/parentCluster1/childCluster0/centroids_parentCluster1_child0_grandChild18.txt";
    vector<float> centroid_childCluster1_0_18;
    readFromCommaDelimitedFile_Float(centroid_childCluster1_0_18PATH, centroid_childCluster1_0_18);
    
    std::string centroid_childCluster1_0_19PATH = casePath + "0.treatedFile/parentCluster1/childCluster0/centroids_parentCluster1_child0_grandChild19.txt";
    vector<float> centroid_childCluster1_0_19;
    readFromCommaDelimitedFile_Float(centroid_childCluster1_0_19PATH, centroid_childCluster1_0_19);
    
    
    cv::Mat r1_0_0 = cv::Mat(numVarANN, 1, CV_32F, centroid_childCluster1_0_0.data());
    cv::Mat r1_0_1 = cv::Mat(numVarANN, 1, CV_32F, centroid_childCluster1_0_1.data());
    cv::Mat r1_0_2 = cv::Mat(numVarANN, 1, CV_32F, centroid_childCluster1_0_2.data());
    cv::Mat r1_0_3 = cv::Mat(numVarANN, 1, CV_32F, centroid_childCluster1_0_3.data());
    cv::Mat r1_0_4 = cv::Mat(numVarANN, 1, CV_32F, centroid_childCluster1_0_4.data());
    cv::Mat r1_0_5 = cv::Mat(numVarANN, 1, CV_32F, centroid_childCluster1_0_5.data());
    cv::Mat r1_0_6 = cv::Mat(numVarANN, 1, CV_32F, centroid_childCluster1_0_6.data());
    cv::Mat r1_0_7 = cv::Mat(numVarANN, 1, CV_32F, centroid_childCluster1_0_7.data());
    cv::Mat r1_0_8 = cv::Mat(numVarANN, 1, CV_32F, centroid_childCluster1_0_8.data());
    cv::Mat r1_0_9 = cv::Mat(numVarANN, 1, CV_32F, centroid_childCluster1_0_9.data());
    cv::Mat r1_0_10 = cv::Mat(numVarANN, 1, CV_32F, centroid_childCluster1_0_10.data());
    cv::Mat r1_0_11 = cv::Mat(numVarANN, 1, CV_32F, centroid_childCluster1_0_11.data());
    cv::Mat r1_0_12 = cv::Mat(numVarANN, 1, CV_32F, centroid_childCluster1_0_12.data());
    cv::Mat r1_0_13 = cv::Mat(numVarANN, 1, CV_32F, centroid_childCluster1_0_13.data());
    cv::Mat r1_0_14 = cv::Mat(numVarANN, 1, CV_32F, centroid_childCluster1_0_14.data());
    cv::Mat r1_0_15 = cv::Mat(numVarANN, 1, CV_32F, centroid_childCluster1_0_15.data());
    cv::Mat r1_0_16 = cv::Mat(numVarANN, 1, CV_32F, centroid_childCluster1_0_16.data());
    cv::Mat r1_0_17 = cv::Mat(numVarANN, 1, CV_32F, centroid_childCluster1_0_17.data());
    cv::Mat r1_0_18 = cv::Mat(numVarANN, 1, CV_32F, centroid_childCluster1_0_18.data());
    cv::Mat r1_0_19 = cv::Mat(numVarANN, 1, CV_32F, centroid_childCluster1_0_19.data());
    
    
    // parentCluster1 - childCluster1 - grandChild
    std::string centroid_childCluster1_1_0PATH = casePath + "0.treatedFile/parentCluster1/childCluster1/centroids_parentCluster1_child1_grandChild0.txt";
    vector<float> centroid_childCluster1_1_0;
    readFromCommaDelimitedFile_Float(centroid_childCluster1_1_0PATH, centroid_childCluster1_1_0);
    
    std::string centroid_childCluster1_1_1PATH = casePath + "0.treatedFile/parentCluster1/childCluster1/centroids_parentCluster1_child1_grandChild1.txt";
    vector<float> centroid_childCluster1_1_1;
    readFromCommaDelimitedFile_Float(centroid_childCluster1_1_1PATH, centroid_childCluster1_1_1);
    
    std::string centroid_childCluster1_1_2PATH = casePath + "0.treatedFile/parentCluster1/childCluster1/centroids_parentCluster1_child1_grandChild2.txt";
    vector<float> centroid_childCluster1_1_2;
    readFromCommaDelimitedFile_Float(centroid_childCluster1_1_2PATH, centroid_childCluster1_1_2);
    
    std::string centroid_childCluster1_1_3PATH = casePath + "0.treatedFile/parentCluster1/childCluster1/centroids_parentCluster1_child1_grandChild3.txt";
    vector<float> centroid_childCluster1_1_3;
    readFromCommaDelimitedFile_Float(centroid_childCluster1_1_3PATH, centroid_childCluster1_1_3);
    
    std::string centroid_childCluster1_1_4PATH = casePath + "0.treatedFile/parentCluster1/childCluster1/centroids_parentCluster1_child1_grandChild4.txt";
    vector<float> centroid_childCluster1_1_4;
    readFromCommaDelimitedFile_Float(centroid_childCluster1_1_4PATH, centroid_childCluster1_1_4);
    
    std::string centroid_childCluster1_1_5PATH = casePath + "0.treatedFile/parentCluster1/childCluster1/centroids_parentCluster1_child1_grandChild5.txt";
    vector<float> centroid_childCluster1_1_5;
    readFromCommaDelimitedFile_Float(centroid_childCluster1_1_5PATH, centroid_childCluster1_1_5);
    
    std::string centroid_childCluster1_1_6PATH = casePath + "0.treatedFile/parentCluster1/childCluster1/centroids_parentCluster1_child1_grandChild6.txt";
    vector<float> centroid_childCluster1_1_6;
    readFromCommaDelimitedFile_Float(centroid_childCluster1_1_6PATH, centroid_childCluster1_1_6);
    
    std::string centroid_childCluster1_1_7PATH = casePath + "0.treatedFile/parentCluster1/childCluster1/centroids_parentCluster1_child1_grandChild7.txt";
    vector<float> centroid_childCluster1_1_7;
    readFromCommaDelimitedFile_Float(centroid_childCluster1_1_7PATH, centroid_childCluster1_1_7);
    
    std::string centroid_childCluster1_1_8PATH = casePath + "0.treatedFile/parentCluster1/childCluster1/centroids_parentCluster1_child1_grandChild8.txt";
    vector<float> centroid_childCluster1_1_8;
    readFromCommaDelimitedFile_Float(centroid_childCluster1_1_8PATH, centroid_childCluster1_1_8);
    
    std::string centroid_childCluster1_1_9PATH = casePath + "0.treatedFile/parentCluster1/childCluster1/centroids_parentCluster1_child1_grandChild9.txt";
    vector<float> centroid_childCluster1_1_9;
    readFromCommaDelimitedFile_Float(centroid_childCluster1_1_9PATH, centroid_childCluster1_1_9);
    
    cv::Mat r1_1_0 = cv::Mat(numVarANN, 1, CV_32F, centroid_childCluster1_1_0.data());
    cv::Mat r1_1_1 = cv::Mat(numVarANN, 1, CV_32F, centroid_childCluster1_1_1.data());
    cv::Mat r1_1_2 = cv::Mat(numVarANN, 1, CV_32F, centroid_childCluster1_1_2.data());
    cv::Mat r1_1_3 = cv::Mat(numVarANN, 1, CV_32F, centroid_childCluster1_1_3.data());
    cv::Mat r1_1_4 = cv::Mat(numVarANN, 1, CV_32F, centroid_childCluster1_1_4.data());
    cv::Mat r1_1_5 = cv::Mat(numVarANN, 1, CV_32F, centroid_childCluster1_1_5.data());
    cv::Mat r1_1_6 = cv::Mat(numVarANN, 1, CV_32F, centroid_childCluster1_1_6.data());
    cv::Mat r1_1_7 = cv::Mat(numVarANN, 1, CV_32F, centroid_childCluster1_1_7.data());
    cv::Mat r1_1_8 = cv::Mat(numVarANN, 1, CV_32F, centroid_childCluster1_1_8.data());
    cv::Mat r1_1_9 = cv::Mat(numVarANN, 1, CV_32F, centroid_childCluster1_1_9.data());

    
    // parentCluster1 - childCluster2 - grandChild
    std::string centroid_childCluster1_2_0PATH = casePath + "0.treatedFile/parentCluster1/childCluster2/centroids_parentCluster1_child2_grandChild0.txt";
    vector<float> centroid_childCluster1_2_0;
    readFromCommaDelimitedFile_Float(centroid_childCluster1_2_0PATH, centroid_childCluster1_2_0);
    
    std::string centroid_childCluster1_2_1PATH = casePath + "0.treatedFile/parentCluster1/childCluster2/centroids_parentCluster1_child2_grandChild1.txt";
    vector<float> centroid_childCluster1_2_1;
    readFromCommaDelimitedFile_Float(centroid_childCluster1_2_1PATH, centroid_childCluster1_2_1);
    
    std::string centroid_childCluster1_2_2PATH = casePath + "0.treatedFile/parentCluster1/childCluster2/centroids_parentCluster1_child2_grandChild2.txt";
    vector<float> centroid_childCluster1_2_2;
    readFromCommaDelimitedFile_Float(centroid_childCluster1_2_2PATH, centroid_childCluster1_2_2);
    
    std::string centroid_childCluster1_2_3PATH = casePath + "0.treatedFile/parentCluster1/childCluster2/centroids_parentCluster1_child2_grandChild3.txt";
    vector<float> centroid_childCluster1_2_3;
    readFromCommaDelimitedFile_Float(centroid_childCluster1_2_3PATH, centroid_childCluster1_2_3);
    
    std::string centroid_childCluster1_2_4PATH = casePath + "0.treatedFile/parentCluster1/childCluster2/centroids_parentCluster1_child2_grandChild4.txt";
    vector<float> centroid_childCluster1_2_4;
    readFromCommaDelimitedFile_Float(centroid_childCluster1_2_4PATH, centroid_childCluster1_2_4);
    
    std::string centroid_childCluster1_2_5PATH = casePath + "0.treatedFile/parentCluster1/childCluster2/centroids_parentCluster1_child2_grandChild5.txt";
    vector<float> centroid_childCluster1_2_5;
    readFromCommaDelimitedFile_Float(centroid_childCluster1_2_5PATH, centroid_childCluster1_2_5);
    
    std::string centroid_childCluster1_2_6PATH = casePath + "0.treatedFile/parentCluster1/childCluster2/centroids_parentCluster1_child2_grandChild6.txt";
    vector<float> centroid_childCluster1_2_6;
    readFromCommaDelimitedFile_Float(centroid_childCluster1_2_6PATH, centroid_childCluster1_2_6);
    
    std::string centroid_childCluster1_2_7PATH = casePath + "0.treatedFile/parentCluster1/childCluster2/centroids_parentCluster1_child2_grandChild7.txt";
    vector<float> centroid_childCluster1_2_7;
    readFromCommaDelimitedFile_Float(centroid_childCluster1_2_7PATH, centroid_childCluster1_2_7);
    
    std::string centroid_childCluster1_2_8PATH = casePath + "0.treatedFile/parentCluster1/childCluster2/centroids_parentCluster1_child2_grandChild8.txt";
    vector<float> centroid_childCluster1_2_8;
    readFromCommaDelimitedFile_Float(centroid_childCluster1_2_8PATH, centroid_childCluster1_2_8);
    
    std::string centroid_childCluster1_2_9PATH = casePath + "0.treatedFile/parentCluster1/childCluster2/centroids_parentCluster1_child2_grandChild9.txt";
    vector<float> centroid_childCluster1_2_9;
    readFromCommaDelimitedFile_Float(centroid_childCluster1_2_9PATH, centroid_childCluster1_2_9);

    
    cv::Mat r1_2_0 = cv::Mat(numVarANN, 1, CV_32F, centroid_childCluster1_2_0.data());
    cv::Mat r1_2_1 = cv::Mat(numVarANN, 1, CV_32F, centroid_childCluster1_2_1.data());
    cv::Mat r1_2_2 = cv::Mat(numVarANN, 1, CV_32F, centroid_childCluster1_2_2.data());
    cv::Mat r1_2_3 = cv::Mat(numVarANN, 1, CV_32F, centroid_childCluster1_2_3.data());
    cv::Mat r1_2_4 = cv::Mat(numVarANN, 1, CV_32F, centroid_childCluster1_2_4.data());
    cv::Mat r1_2_5 = cv::Mat(numVarANN, 1, CV_32F, centroid_childCluster1_2_5.data());
    cv::Mat r1_2_6 = cv::Mat(numVarANN, 1, CV_32F, centroid_childCluster1_2_6.data());
    cv::Mat r1_2_7 = cv::Mat(numVarANN, 1, CV_32F, centroid_childCluster1_2_7.data());
    cv::Mat r1_2_8 = cv::Mat(numVarANN, 1, CV_32F, centroid_childCluster1_2_8.data());
    cv::Mat r1_2_9 = cv::Mat(numVarANN, 1, CV_32F, centroid_childCluster1_2_9.data());

    
    //Declare CHILD LOCAL INPUT & OUTPUT Standardize for CHILD clusters: standardizedX = (X - meanScale) / standardDeviationScale
    // Child of parentCluster0
    std::string cluster0_0_Folder = "ANNReg_parentCluster0/ANNReg_childCluster0_0/";
    std::string cluster0_0_meanScaleINPUTPATH = casePath + cluster0_0_Folder + "scalerLocal_Input_onlyMean.txt";
    std::string cluster0_0_stdScaleINPUTPATH = casePath + cluster0_0_Folder + "scalerLocal_Input_onlyStandardDeviation.txt";
    std::string cluster0_0_meanScaleOUTPUTPATH = casePath + cluster0_0_Folder + "scalerLocal_Label_onlyMean.txt";
    std::string cluster0_0_stdScaleOUTPUTPATH = casePath + cluster0_0_Folder + "scalerLocal_Label_onlyStandardDeviation.txt";
    vector<float> cluster0_0_meanScaleINPUT;
    vector<float> cluster0_0_stdScaleINPUT;
    vector<float> cluster0_0_meanScaleOUTPUT;
    vector<float> cluster0_0_stdScaleOUTPUT;
    readFromCommaDelimitedFile_Float(cluster0_0_meanScaleINPUTPATH, cluster0_0_meanScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster0_0_stdScaleINPUTPATH, cluster0_0_stdScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster0_0_meanScaleOUTPUTPATH, cluster0_0_meanScaleOUTPUT);
    readFromCommaDelimitedFile_Float(cluster0_0_stdScaleOUTPUTPATH, cluster0_0_stdScaleOUTPUT);
    
    std::string cluster0_1_Folder = "ANNReg_parentCluster0/ANNReg_childCluster0_1/";
    std::string cluster0_1_meanScaleINPUTPATH = casePath + cluster0_1_Folder + "scalerLocal_Input_onlyMean.txt";
    std::string cluster0_1_stdScaleINPUTPATH = casePath + cluster0_1_Folder + "scalerLocal_Input_onlyStandardDeviation.txt";
    std::string cluster0_1_meanScaleOUTPUTPATH = casePath + cluster0_1_Folder + "scalerLocal_Label_onlyMean.txt";
    std::string cluster0_1_stdScaleOUTPUTPATH = casePath + cluster0_1_Folder + "scalerLocal_Label_onlyStandardDeviation.txt";
    vector<float> cluster0_1_meanScaleINPUT;
    vector<float> cluster0_1_stdScaleINPUT;
    vector<float> cluster0_1_meanScaleOUTPUT;
    vector<float> cluster0_1_stdScaleOUTPUT;
    readFromCommaDelimitedFile_Float(cluster0_1_meanScaleINPUTPATH, cluster0_1_meanScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster0_1_stdScaleINPUTPATH, cluster0_1_stdScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster0_1_meanScaleOUTPUTPATH, cluster0_1_meanScaleOUTPUT);
    readFromCommaDelimitedFile_Float(cluster0_1_stdScaleOUTPUTPATH, cluster0_1_stdScaleOUTPUT);
    
    std::string cluster0_2_Folder = "ANNReg_parentCluster0/ANNReg_childCluster0_2/";
    std::string cluster0_2_meanScaleINPUTPATH = casePath + cluster0_2_Folder + "scalerLocal_Input_onlyMean.txt";
    std::string cluster0_2_stdScaleINPUTPATH = casePath + cluster0_2_Folder + "scalerLocal_Input_onlyStandardDeviation.txt";
    std::string cluster0_2_meanScaleOUTPUTPATH = casePath + cluster0_2_Folder + "scalerLocal_Label_onlyMean.txt";
    std::string cluster0_2_stdScaleOUTPUTPATH = casePath + cluster0_2_Folder + "scalerLocal_Label_onlyStandardDeviation.txt";
    vector<float> cluster0_2_meanScaleINPUT;
    vector<float> cluster0_2_stdScaleINPUT;
    vector<float> cluster0_2_meanScaleOUTPUT;
    vector<float> cluster0_2_stdScaleOUTPUT;
    readFromCommaDelimitedFile_Float(cluster0_2_meanScaleINPUTPATH, cluster0_2_meanScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster0_2_stdScaleINPUTPATH, cluster0_2_stdScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster0_2_meanScaleOUTPUTPATH, cluster0_2_meanScaleOUTPUT);
    readFromCommaDelimitedFile_Float(cluster0_2_stdScaleOUTPUTPATH, cluster0_2_stdScaleOUTPUT);
    
    std::string cluster0_3_Folder = "ANNReg_parentCluster0/ANNReg_childCluster0_3/";
    std::string cluster0_3_meanScaleINPUTPATH = casePath + cluster0_3_Folder + "scalerLocal_Input_onlyMean.txt";
    std::string cluster0_3_stdScaleINPUTPATH = casePath + cluster0_3_Folder + "scalerLocal_Input_onlyStandardDeviation.txt";
    std::string cluster0_3_meanScaleOUTPUTPATH = casePath + cluster0_3_Folder + "scalerLocal_Label_onlyMean.txt";
    std::string cluster0_3_stdScaleOUTPUTPATH = casePath + cluster0_3_Folder + "scalerLocal_Label_onlyStandardDeviation.txt";
    vector<float> cluster0_3_meanScaleINPUT;
    vector<float> cluster0_3_stdScaleINPUT;
    vector<float> cluster0_3_meanScaleOUTPUT;
    vector<float> cluster0_3_stdScaleOUTPUT;
    readFromCommaDelimitedFile_Float(cluster0_3_meanScaleINPUTPATH, cluster0_3_meanScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster0_3_stdScaleINPUTPATH, cluster0_3_stdScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster0_3_meanScaleOUTPUTPATH, cluster0_3_meanScaleOUTPUT);
    readFromCommaDelimitedFile_Float(cluster0_3_stdScaleOUTPUTPATH, cluster0_3_stdScaleOUTPUT);
    
    std::string cluster0_4_Folder = "ANNReg_parentCluster0/ANNReg_childCluster0_4/";
    std::string cluster0_4_meanScaleINPUTPATH = casePath + cluster0_4_Folder + "scalerLocal_Input_onlyMean.txt";
    std::string cluster0_4_stdScaleINPUTPATH = casePath + cluster0_4_Folder + "scalerLocal_Input_onlyStandardDeviation.txt";
    std::string cluster0_4_meanScaleOUTPUTPATH = casePath + cluster0_4_Folder + "scalerLocal_Label_onlyMean.txt";
    std::string cluster0_4_stdScaleOUTPUTPATH = casePath + cluster0_4_Folder + "scalerLocal_Label_onlyStandardDeviation.txt";
    vector<float> cluster0_4_meanScaleINPUT;
    vector<float> cluster0_4_stdScaleINPUT;
    vector<float> cluster0_4_meanScaleOUTPUT;
    vector<float> cluster0_4_stdScaleOUTPUT;
    readFromCommaDelimitedFile_Float(cluster0_4_meanScaleINPUTPATH, cluster0_4_meanScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster0_4_stdScaleINPUTPATH, cluster0_4_stdScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster0_4_meanScaleOUTPUTPATH, cluster0_4_meanScaleOUTPUT);
    readFromCommaDelimitedFile_Float(cluster0_4_stdScaleOUTPUTPATH, cluster0_4_stdScaleOUTPUT);
    
    std::string cluster0_5_Folder = "ANNReg_parentCluster0/ANNReg_childCluster0_5/";
    std::string cluster0_5_meanScaleINPUTPATH = casePath + cluster0_5_Folder + "scalerLocal_Input_onlyMean.txt";
    std::string cluster0_5_stdScaleINPUTPATH = casePath + cluster0_5_Folder + "scalerLocal_Input_onlyStandardDeviation.txt";
    std::string cluster0_5_meanScaleOUTPUTPATH = casePath + cluster0_5_Folder + "scalerLocal_Label_onlyMean.txt";
    std::string cluster0_5_stdScaleOUTPUTPATH = casePath + cluster0_5_Folder + "scalerLocal_Label_onlyStandardDeviation.txt";
    vector<float> cluster0_5_meanScaleINPUT;
    vector<float> cluster0_5_stdScaleINPUT;
    vector<float> cluster0_5_meanScaleOUTPUT;
    vector<float> cluster0_5_stdScaleOUTPUT;
    readFromCommaDelimitedFile_Float(cluster0_5_meanScaleINPUTPATH, cluster0_5_meanScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster0_5_stdScaleINPUTPATH, cluster0_5_stdScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster0_5_meanScaleOUTPUTPATH, cluster0_5_meanScaleOUTPUT);
    readFromCommaDelimitedFile_Float(cluster0_5_stdScaleOUTPUTPATH, cluster0_5_stdScaleOUTPUT);
    
    std::string cluster0_6_Folder = "ANNReg_parentCluster0/ANNReg_childCluster0_6/";
    std::string cluster0_6_meanScaleINPUTPATH = casePath + cluster0_6_Folder + "scalerLocal_Input_onlyMean.txt";
    std::string cluster0_6_stdScaleINPUTPATH = casePath + cluster0_6_Folder + "scalerLocal_Input_onlyStandardDeviation.txt";
    std::string cluster0_6_meanScaleOUTPUTPATH = casePath + cluster0_6_Folder + "scalerLocal_Label_onlyMean.txt";
    std::string cluster0_6_stdScaleOUTPUTPATH = casePath + cluster0_6_Folder + "scalerLocal_Label_onlyStandardDeviation.txt";
    vector<float> cluster0_6_meanScaleINPUT;
    vector<float> cluster0_6_stdScaleINPUT;
    vector<float> cluster0_6_meanScaleOUTPUT;
    vector<float> cluster0_6_stdScaleOUTPUT;
    readFromCommaDelimitedFile_Float(cluster0_6_meanScaleINPUTPATH, cluster0_6_meanScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster0_6_stdScaleINPUTPATH, cluster0_6_stdScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster0_6_meanScaleOUTPUTPATH, cluster0_6_meanScaleOUTPUT);
    readFromCommaDelimitedFile_Float(cluster0_6_stdScaleOUTPUTPATH, cluster0_6_stdScaleOUTPUT);
    
    std::string cluster0_7_Folder = "ANNReg_parentCluster0/ANNReg_childCluster0_7/";
    std::string cluster0_7_meanScaleINPUTPATH = casePath + cluster0_7_Folder + "scalerLocal_Input_onlyMean.txt";
    std::string cluster0_7_stdScaleINPUTPATH = casePath + cluster0_7_Folder + "scalerLocal_Input_onlyStandardDeviation.txt";
    std::string cluster0_7_meanScaleOUTPUTPATH = casePath + cluster0_7_Folder + "scalerLocal_Label_onlyMean.txt";
    std::string cluster0_7_stdScaleOUTPUTPATH = casePath + cluster0_7_Folder + "scalerLocal_Label_onlyStandardDeviation.txt";
    vector<float> cluster0_7_meanScaleINPUT;
    vector<float> cluster0_7_stdScaleINPUT;
    vector<float> cluster0_7_meanScaleOUTPUT;
    vector<float> cluster0_7_stdScaleOUTPUT;
    readFromCommaDelimitedFile_Float(cluster0_7_meanScaleINPUTPATH, cluster0_7_meanScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster0_7_stdScaleINPUTPATH, cluster0_7_stdScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster0_7_meanScaleOUTPUTPATH, cluster0_7_meanScaleOUTPUT);
    readFromCommaDelimitedFile_Float(cluster0_7_stdScaleOUTPUTPATH, cluster0_7_stdScaleOUTPUT);
    
    std::string cluster0_8_Folder = "ANNReg_parentCluster0/ANNReg_childCluster0_8/";
    std::string cluster0_8_meanScaleINPUTPATH = casePath + cluster0_8_Folder + "scalerLocal_Input_onlyMean.txt";
    std::string cluster0_8_stdScaleINPUTPATH = casePath + cluster0_8_Folder + "scalerLocal_Input_onlyStandardDeviation.txt";
    std::string cluster0_8_meanScaleOUTPUTPATH = casePath + cluster0_8_Folder + "scalerLocal_Label_onlyMean.txt";
    std::string cluster0_8_stdScaleOUTPUTPATH = casePath + cluster0_8_Folder + "scalerLocal_Label_onlyStandardDeviation.txt";
    vector<float> cluster0_8_meanScaleINPUT;
    vector<float> cluster0_8_stdScaleINPUT;
    vector<float> cluster0_8_meanScaleOUTPUT;
    vector<float> cluster0_8_stdScaleOUTPUT;
    readFromCommaDelimitedFile_Float(cluster0_8_meanScaleINPUTPATH, cluster0_8_meanScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster0_8_stdScaleINPUTPATH, cluster0_8_stdScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster0_8_meanScaleOUTPUTPATH, cluster0_8_meanScaleOUTPUT);
    readFromCommaDelimitedFile_Float(cluster0_8_stdScaleOUTPUTPATH, cluster0_8_stdScaleOUTPUT);
    
    std::string cluster0_9_Folder = "ANNReg_parentCluster0/ANNReg_childCluster0_9/";
    std::string cluster0_9_meanScaleINPUTPATH = casePath + cluster0_9_Folder + "scalerLocal_Input_onlyMean.txt";
    std::string cluster0_9_stdScaleINPUTPATH = casePath + cluster0_9_Folder + "scalerLocal_Input_onlyStandardDeviation.txt";
    std::string cluster0_9_meanScaleOUTPUTPATH = casePath + cluster0_9_Folder + "scalerLocal_Label_onlyMean.txt";
    std::string cluster0_9_stdScaleOUTPUTPATH = casePath + cluster0_9_Folder + "scalerLocal_Label_onlyStandardDeviation.txt";
    vector<float> cluster0_9_meanScaleINPUT;
    vector<float> cluster0_9_stdScaleINPUT;
    vector<float> cluster0_9_meanScaleOUTPUT;
    vector<float> cluster0_9_stdScaleOUTPUT;
    readFromCommaDelimitedFile_Float(cluster0_9_meanScaleINPUTPATH, cluster0_9_meanScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster0_9_stdScaleINPUTPATH, cluster0_9_stdScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster0_9_meanScaleOUTPUTPATH, cluster0_9_meanScaleOUTPUT);
    readFromCommaDelimitedFile_Float(cluster0_9_stdScaleOUTPUTPATH, cluster0_9_stdScaleOUTPUT);
    
    std::string cluster0_10_Folder = "ANNReg_parentCluster0/ANNReg_childCluster0_10/";
    std::string cluster0_10_meanScaleINPUTPATH = casePath + cluster0_10_Folder + "scalerLocal_Input_onlyMean.txt";
    std::string cluster0_10_stdScaleINPUTPATH = casePath + cluster0_10_Folder + "scalerLocal_Input_onlyStandardDeviation.txt";
    std::string cluster0_10_meanScaleOUTPUTPATH = casePath + cluster0_10_Folder + "scalerLocal_Label_onlyMean.txt";
    std::string cluster0_10_stdScaleOUTPUTPATH = casePath + cluster0_10_Folder + "scalerLocal_Label_onlyStandardDeviation.txt";
    vector<float> cluster0_10_meanScaleINPUT;
    vector<float> cluster0_10_stdScaleINPUT;
    vector<float> cluster0_10_meanScaleOUTPUT;
    vector<float> cluster0_10_stdScaleOUTPUT;
    readFromCommaDelimitedFile_Float(cluster0_10_meanScaleINPUTPATH, cluster0_10_meanScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster0_10_stdScaleINPUTPATH, cluster0_10_stdScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster0_10_meanScaleOUTPUTPATH, cluster0_10_meanScaleOUTPUT);
    readFromCommaDelimitedFile_Float(cluster0_10_stdScaleOUTPUTPATH, cluster0_10_stdScaleOUTPUT);
    
    std::string cluster0_11_Folder = "ANNReg_parentCluster0/ANNReg_childCluster0_11/";
    std::string cluster0_11_meanScaleINPUTPATH = casePath + cluster0_11_Folder + "scalerLocal_Input_onlyMean.txt";
    std::string cluster0_11_stdScaleINPUTPATH = casePath + cluster0_11_Folder + "scalerLocal_Input_onlyStandardDeviation.txt";
    std::string cluster0_11_meanScaleOUTPUTPATH = casePath + cluster0_11_Folder + "scalerLocal_Label_onlyMean.txt";
    std::string cluster0_11_stdScaleOUTPUTPATH = casePath + cluster0_11_Folder + "scalerLocal_Label_onlyStandardDeviation.txt";
    vector<float> cluster0_11_meanScaleINPUT;
    vector<float> cluster0_11_stdScaleINPUT;
    vector<float> cluster0_11_meanScaleOUTPUT;
    vector<float> cluster0_11_stdScaleOUTPUT;
    readFromCommaDelimitedFile_Float(cluster0_11_meanScaleINPUTPATH, cluster0_11_meanScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster0_11_stdScaleINPUTPATH, cluster0_11_stdScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster0_11_meanScaleOUTPUTPATH, cluster0_11_meanScaleOUTPUT);
    readFromCommaDelimitedFile_Float(cluster0_11_stdScaleOUTPUTPATH, cluster0_11_stdScaleOUTPUT);
    
    std::string cluster0_12_Folder = "ANNReg_parentCluster0/ANNReg_childCluster0_12/";
    std::string cluster0_12_meanScaleINPUTPATH = casePath + cluster0_12_Folder + "scalerLocal_Input_onlyMean.txt";
    std::string cluster0_12_stdScaleINPUTPATH = casePath + cluster0_12_Folder + "scalerLocal_Input_onlyStandardDeviation.txt";
    std::string cluster0_12_meanScaleOUTPUTPATH = casePath + cluster0_12_Folder + "scalerLocal_Label_onlyMean.txt";
    std::string cluster0_12_stdScaleOUTPUTPATH = casePath + cluster0_12_Folder + "scalerLocal_Label_onlyStandardDeviation.txt";
    vector<float> cluster0_12_meanScaleINPUT;
    vector<float> cluster0_12_stdScaleINPUT;
    vector<float> cluster0_12_meanScaleOUTPUT;
    vector<float> cluster0_12_stdScaleOUTPUT;
    readFromCommaDelimitedFile_Float(cluster0_12_meanScaleINPUTPATH, cluster0_12_meanScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster0_12_stdScaleINPUTPATH, cluster0_12_stdScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster0_12_meanScaleOUTPUTPATH, cluster0_12_meanScaleOUTPUT);
    readFromCommaDelimitedFile_Float(cluster0_12_stdScaleOUTPUTPATH, cluster0_12_stdScaleOUTPUT);
    
    std::string cluster0_13_Folder = "ANNReg_parentCluster0/ANNReg_childCluster0_13/";
    std::string cluster0_13_meanScaleINPUTPATH = casePath + cluster0_13_Folder + "scalerLocal_Input_onlyMean.txt";
    std::string cluster0_13_stdScaleINPUTPATH = casePath + cluster0_13_Folder + "scalerLocal_Input_onlyStandardDeviation.txt";
    std::string cluster0_13_meanScaleOUTPUTPATH = casePath + cluster0_13_Folder + "scalerLocal_Label_onlyMean.txt";
    std::string cluster0_13_stdScaleOUTPUTPATH = casePath + cluster0_13_Folder + "scalerLocal_Label_onlyStandardDeviation.txt";
    vector<float> cluster0_13_meanScaleINPUT;
    vector<float> cluster0_13_stdScaleINPUT;
    vector<float> cluster0_13_meanScaleOUTPUT;
    vector<float> cluster0_13_stdScaleOUTPUT;
    readFromCommaDelimitedFile_Float(cluster0_13_meanScaleINPUTPATH, cluster0_13_meanScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster0_13_stdScaleINPUTPATH, cluster0_13_stdScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster0_13_meanScaleOUTPUTPATH, cluster0_13_meanScaleOUTPUT);
    readFromCommaDelimitedFile_Float(cluster0_13_stdScaleOUTPUTPATH, cluster0_13_stdScaleOUTPUT);
    
    std::string cluster0_14_Folder = "ANNReg_parentCluster0/ANNReg_childCluster0_14/";
    std::string cluster0_14_meanScaleINPUTPATH = casePath + cluster0_14_Folder + "scalerLocal_Input_onlyMean.txt";
    std::string cluster0_14_stdScaleINPUTPATH = casePath + cluster0_14_Folder + "scalerLocal_Input_onlyStandardDeviation.txt";
    std::string cluster0_14_meanScaleOUTPUTPATH = casePath + cluster0_14_Folder + "scalerLocal_Label_onlyMean.txt";
    std::string cluster0_14_stdScaleOUTPUTPATH = casePath + cluster0_14_Folder + "scalerLocal_Label_onlyStandardDeviation.txt";
    vector<float> cluster0_14_meanScaleINPUT;
    vector<float> cluster0_14_stdScaleINPUT;
    vector<float> cluster0_14_meanScaleOUTPUT;
    vector<float> cluster0_14_stdScaleOUTPUT;
    readFromCommaDelimitedFile_Float(cluster0_14_meanScaleINPUTPATH, cluster0_14_meanScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster0_14_stdScaleINPUTPATH, cluster0_14_stdScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster0_14_meanScaleOUTPUTPATH, cluster0_14_meanScaleOUTPUT);
    readFromCommaDelimitedFile_Float(cluster0_14_stdScaleOUTPUTPATH, cluster0_14_stdScaleOUTPUT);
    
    std::string cluster0_15_Folder = "ANNReg_parentCluster0/ANNReg_childCluster0_15/";
    std::string cluster0_15_meanScaleINPUTPATH = casePath + cluster0_15_Folder + "scalerLocal_Input_onlyMean.txt";
    std::string cluster0_15_stdScaleINPUTPATH = casePath + cluster0_15_Folder + "scalerLocal_Input_onlyStandardDeviation.txt";
    std::string cluster0_15_meanScaleOUTPUTPATH = casePath + cluster0_15_Folder + "scalerLocal_Label_onlyMean.txt";
    std::string cluster0_15_stdScaleOUTPUTPATH = casePath + cluster0_15_Folder + "scalerLocal_Label_onlyStandardDeviation.txt";
    vector<float> cluster0_15_meanScaleINPUT;
    vector<float> cluster0_15_stdScaleINPUT;
    vector<float> cluster0_15_meanScaleOUTPUT;
    vector<float> cluster0_15_stdScaleOUTPUT;
    readFromCommaDelimitedFile_Float(cluster0_15_meanScaleINPUTPATH, cluster0_15_meanScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster0_15_stdScaleINPUTPATH, cluster0_15_stdScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster0_15_meanScaleOUTPUTPATH, cluster0_15_meanScaleOUTPUT);
    readFromCommaDelimitedFile_Float(cluster0_15_stdScaleOUTPUTPATH, cluster0_15_stdScaleOUTPUT);
    
    std::string cluster0_16_Folder = "ANNReg_parentCluster0/ANNReg_childCluster0_16/";
    std::string cluster0_16_meanScaleINPUTPATH = casePath + cluster0_16_Folder + "scalerLocal_Input_onlyMean.txt";
    std::string cluster0_16_stdScaleINPUTPATH = casePath + cluster0_16_Folder + "scalerLocal_Input_onlyStandardDeviation.txt";
    std::string cluster0_16_meanScaleOUTPUTPATH = casePath + cluster0_16_Folder + "scalerLocal_Label_onlyMean.txt";
    std::string cluster0_16_stdScaleOUTPUTPATH = casePath + cluster0_16_Folder + "scalerLocal_Label_onlyStandardDeviation.txt";
    vector<float> cluster0_16_meanScaleINPUT;
    vector<float> cluster0_16_stdScaleINPUT;
    vector<float> cluster0_16_meanScaleOUTPUT;
    vector<float> cluster0_16_stdScaleOUTPUT;
    readFromCommaDelimitedFile_Float(cluster0_16_meanScaleINPUTPATH, cluster0_16_meanScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster0_16_stdScaleINPUTPATH, cluster0_16_stdScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster0_16_meanScaleOUTPUTPATH, cluster0_16_meanScaleOUTPUT);
    readFromCommaDelimitedFile_Float(cluster0_16_stdScaleOUTPUTPATH, cluster0_16_stdScaleOUTPUT);
    
    std::string cluster0_17_Folder = "ANNReg_parentCluster0/ANNReg_childCluster0_17/";
    std::string cluster0_17_meanScaleINPUTPATH = casePath + cluster0_17_Folder + "scalerLocal_Input_onlyMean.txt";
    std::string cluster0_17_stdScaleINPUTPATH = casePath + cluster0_17_Folder + "scalerLocal_Input_onlyStandardDeviation.txt";
    std::string cluster0_17_meanScaleOUTPUTPATH = casePath + cluster0_17_Folder + "scalerLocal_Label_onlyMean.txt";
    std::string cluster0_17_stdScaleOUTPUTPATH = casePath + cluster0_17_Folder + "scalerLocal_Label_onlyStandardDeviation.txt";
    vector<float> cluster0_17_meanScaleINPUT;
    vector<float> cluster0_17_stdScaleINPUT;
    vector<float> cluster0_17_meanScaleOUTPUT;
    vector<float> cluster0_17_stdScaleOUTPUT;
    readFromCommaDelimitedFile_Float(cluster0_17_meanScaleINPUTPATH, cluster0_17_meanScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster0_17_stdScaleINPUTPATH, cluster0_17_stdScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster0_17_meanScaleOUTPUTPATH, cluster0_17_meanScaleOUTPUT);
    readFromCommaDelimitedFile_Float(cluster0_17_stdScaleOUTPUTPATH, cluster0_17_stdScaleOUTPUT);
    
    std::string cluster0_18_Folder = "ANNReg_parentCluster0/ANNReg_childCluster0_18/";
    std::string cluster0_18_meanScaleINPUTPATH = casePath + cluster0_18_Folder + "scalerLocal_Input_onlyMean.txt";
    std::string cluster0_18_stdScaleINPUTPATH = casePath + cluster0_18_Folder + "scalerLocal_Input_onlyStandardDeviation.txt";
    std::string cluster0_18_meanScaleOUTPUTPATH = casePath + cluster0_18_Folder + "scalerLocal_Label_onlyMean.txt";
    std::string cluster0_18_stdScaleOUTPUTPATH = casePath + cluster0_18_Folder + "scalerLocal_Label_onlyStandardDeviation.txt";
    vector<float> cluster0_18_meanScaleINPUT;
    vector<float> cluster0_18_stdScaleINPUT;
    vector<float> cluster0_18_meanScaleOUTPUT;
    vector<float> cluster0_18_stdScaleOUTPUT;
    readFromCommaDelimitedFile_Float(cluster0_18_meanScaleINPUTPATH, cluster0_18_meanScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster0_18_stdScaleINPUTPATH, cluster0_18_stdScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster0_18_meanScaleOUTPUTPATH, cluster0_18_meanScaleOUTPUT);
    readFromCommaDelimitedFile_Float(cluster0_18_stdScaleOUTPUTPATH, cluster0_18_stdScaleOUTPUT);
    
    std::string cluster0_19_Folder = "ANNReg_parentCluster0/ANNReg_childCluster0_19/";
    std::string cluster0_19_meanScaleINPUTPATH = casePath + cluster0_19_Folder + "scalerLocal_Input_onlyMean.txt";
    std::string cluster0_19_stdScaleINPUTPATH = casePath + cluster0_19_Folder + "scalerLocal_Input_onlyStandardDeviation.txt";
    std::string cluster0_19_meanScaleOUTPUTPATH = casePath + cluster0_19_Folder + "scalerLocal_Label_onlyMean.txt";
    std::string cluster0_19_stdScaleOUTPUTPATH = casePath + cluster0_19_Folder + "scalerLocal_Label_onlyStandardDeviation.txt";
    vector<float> cluster0_19_meanScaleINPUT;
    vector<float> cluster0_19_stdScaleINPUT;
    vector<float> cluster0_19_meanScaleOUTPUT;
    vector<float> cluster0_19_stdScaleOUTPUT;
    readFromCommaDelimitedFile_Float(cluster0_19_meanScaleINPUTPATH, cluster0_19_meanScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster0_19_stdScaleINPUTPATH, cluster0_19_stdScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster0_19_meanScaleOUTPUTPATH, cluster0_19_meanScaleOUTPUT);
    readFromCommaDelimitedFile_Float(cluster0_19_stdScaleOUTPUTPATH, cluster0_19_stdScaleOUTPUT);
    
    
    // Child of parentCluster1: child1_0, child1_1, child1_2
    std::string cluster1_0_meanScaleINPUTPATH = casePath + "0.treatedFile/parentCluster1/childCluster0/scalerParentCluster1_childCluster0_OnlyMean.txt";
    vector<float> cluster1_0_meanScaleINPUT;
    readFromCommaDelimitedFile_Float(cluster1_0_meanScaleINPUTPATH, cluster1_0_meanScaleINPUT);
    std::string cluster1_0_stdScaleINPUTPATH = casePath + "0.treatedFile/parentCluster1/childCluster0/scalerParentCluster1_childCluster0_OnlyStandardDeviation.txt";
    vector<float> cluster1_0_stdScaleINPUT;
    readFromCommaDelimitedFile_Float(cluster1_0_stdScaleINPUTPATH, cluster1_0_stdScaleINPUT);
    
    std::string cluster1_1_meanScaleINPUTPATH = casePath + "0.treatedFile/parentCluster1/childCluster1/scalerParentCluster1_childCluster1_OnlyMean.txt";
    vector<float> cluster1_1_meanScaleINPUT;
    readFromCommaDelimitedFile_Float(cluster1_1_meanScaleINPUTPATH, cluster1_1_meanScaleINPUT);
    std::string cluster1_1_stdScaleINPUTPATH = casePath + "0.treatedFile/parentCluster1/childCluster1/scalerParentCluster1_childCluster1_OnlyStandardDeviation.txt";
    vector<float> cluster1_1_stdScaleINPUT;
    readFromCommaDelimitedFile_Float(cluster1_1_stdScaleINPUTPATH, cluster1_1_stdScaleINPUT);
    
    std::string cluster1_2_meanScaleINPUTPATH = casePath + "0.treatedFile/parentCluster1/childCluster2/scalerParentCluster1_childCluster2_OnlyMean.txt";
    vector<float> cluster1_2_meanScaleINPUT;
    readFromCommaDelimitedFile_Float(cluster1_2_meanScaleINPUTPATH, cluster1_2_meanScaleINPUT);
    std::string cluster1_2_stdScaleINPUTPATH = casePath + "0.treatedFile/parentCluster1/childCluster2/scalerParentCluster1_childCluster2_OnlyStandardDeviation.txt";
    vector<float> cluster1_2_stdScaleINPUT;
    readFromCommaDelimitedFile_Float(cluster1_2_stdScaleINPUTPATH, cluster1_2_stdScaleINPUT);
    

    // Grand child of parentCluster1
        // parentCluster1 - childCluster0
    std::string cluster1_0_0_Folder = "ANNReg_parentCluster1/ANNReg_childCluster1_0_40grandChild/grandChild1_0_0/";
    std::string cluster1_0_0_meanScaleINPUTPATH = casePath + cluster1_0_0_Folder + "scalerLocal_Input_onlyMean.txt";
    std::string cluster1_0_0_stdScaleINPUTPATH = casePath + cluster1_0_0_Folder + "scalerLocal_Input_onlyStandardDeviation.txt";
    std::string cluster1_0_0_meanScaleOUTPUTPATH = casePath + cluster1_0_0_Folder + "scalerLocal_Label_onlyMean.txt";
    std::string cluster1_0_0_stdScaleOUTPUTPATH = casePath + cluster1_0_0_Folder + "scalerLocal_Label_onlyStandardDeviation.txt";
    vector<float> cluster1_0_0_meanScaleINPUT;
    vector<float> cluster1_0_0_stdScaleINPUT;
    vector<float> cluster1_0_0_meanScaleOUTPUT;
    vector<float> cluster1_0_0_stdScaleOUTPUT;
    readFromCommaDelimitedFile_Float(cluster1_0_0_meanScaleINPUTPATH, cluster1_0_0_meanScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster1_0_0_stdScaleINPUTPATH, cluster1_0_0_stdScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster1_0_0_meanScaleOUTPUTPATH, cluster1_0_0_meanScaleOUTPUT);
    readFromCommaDelimitedFile_Float(cluster1_0_0_stdScaleOUTPUTPATH, cluster1_0_0_stdScaleOUTPUT);
    
    std::string cluster1_0_1_Folder = "ANNReg_parentCluster1/ANNReg_childCluster1_0_40grandChild/grandChild1_0_1/";
    std::string cluster1_0_1_meanScaleINPUTPATH = casePath + cluster1_0_1_Folder + "scalerLocal_Input_onlyMean.txt";
    std::string cluster1_0_1_stdScaleINPUTPATH = casePath + cluster1_0_1_Folder + "scalerLocal_Input_onlyStandardDeviation.txt";
    std::string cluster1_0_1_meanScaleOUTPUTPATH = casePath + cluster1_0_1_Folder + "scalerLocal_Label_onlyMean.txt";
    std::string cluster1_0_1_stdScaleOUTPUTPATH = casePath + cluster1_0_1_Folder + "scalerLocal_Label_onlyStandardDeviation.txt";
    vector<float> cluster1_0_1_meanScaleINPUT;
    vector<float> cluster1_0_1_stdScaleINPUT;
    vector<float> cluster1_0_1_meanScaleOUTPUT;
    vector<float> cluster1_0_1_stdScaleOUTPUT;
    readFromCommaDelimitedFile_Float(cluster1_0_1_meanScaleINPUTPATH, cluster1_0_1_meanScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster1_0_1_stdScaleINPUTPATH, cluster1_0_1_stdScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster1_0_1_meanScaleOUTPUTPATH, cluster1_0_1_meanScaleOUTPUT);
    readFromCommaDelimitedFile_Float(cluster1_0_1_stdScaleOUTPUTPATH, cluster1_0_1_stdScaleOUTPUT);
    
    std::string cluster1_0_2_Folder = "ANNReg_parentCluster1/ANNReg_childCluster1_0_40grandChild/grandChild1_0_2/";
    std::string cluster1_0_2_meanScaleINPUTPATH = casePath + cluster1_0_2_Folder + "scalerLocal_Input_onlyMean.txt";
    std::string cluster1_0_2_stdScaleINPUTPATH = casePath + cluster1_0_2_Folder + "scalerLocal_Input_onlyStandardDeviation.txt";
    std::string cluster1_0_2_meanScaleOUTPUTPATH = casePath + cluster1_0_2_Folder + "scalerLocal_Label_onlyMean.txt";
    std::string cluster1_0_2_stdScaleOUTPUTPATH = casePath + cluster1_0_2_Folder + "scalerLocal_Label_onlyStandardDeviation.txt";
    vector<float> cluster1_0_2_meanScaleINPUT;
    vector<float> cluster1_0_2_stdScaleINPUT;
    vector<float> cluster1_0_2_meanScaleOUTPUT;
    vector<float> cluster1_0_2_stdScaleOUTPUT;
    readFromCommaDelimitedFile_Float(cluster1_0_2_meanScaleINPUTPATH, cluster1_0_2_meanScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster1_0_2_stdScaleINPUTPATH, cluster1_0_2_stdScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster1_0_2_meanScaleOUTPUTPATH, cluster1_0_2_meanScaleOUTPUT);
    readFromCommaDelimitedFile_Float(cluster1_0_2_stdScaleOUTPUTPATH, cluster1_0_2_stdScaleOUTPUT);
    
    std::string cluster1_0_3_Folder = "ANNReg_parentCluster1/ANNReg_childCluster1_0_40grandChild/grandChild1_0_3/";
    std::string cluster1_0_3_meanScaleINPUTPATH = casePath + cluster1_0_3_Folder + "scalerLocal_Input_onlyMean.txt";
    std::string cluster1_0_3_stdScaleINPUTPATH = casePath + cluster1_0_3_Folder + "scalerLocal_Input_onlyStandardDeviation.txt";
    std::string cluster1_0_3_meanScaleOUTPUTPATH = casePath + cluster1_0_3_Folder + "scalerLocal_Label_onlyMean.txt";
    std::string cluster1_0_3_stdScaleOUTPUTPATH = casePath + cluster1_0_3_Folder + "scalerLocal_Label_onlyStandardDeviation.txt";
    vector<float> cluster1_0_3_meanScaleINPUT;
    vector<float> cluster1_0_3_stdScaleINPUT;
    vector<float> cluster1_0_3_meanScaleOUTPUT;
    vector<float> cluster1_0_3_stdScaleOUTPUT;
    readFromCommaDelimitedFile_Float(cluster1_0_3_meanScaleINPUTPATH, cluster1_0_3_meanScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster1_0_3_stdScaleINPUTPATH, cluster1_0_3_stdScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster1_0_3_meanScaleOUTPUTPATH, cluster1_0_3_meanScaleOUTPUT);
    readFromCommaDelimitedFile_Float(cluster1_0_3_stdScaleOUTPUTPATH, cluster1_0_3_stdScaleOUTPUT);
    
    std::string cluster1_0_4_Folder = "ANNReg_parentCluster1/ANNReg_childCluster1_0_40grandChild/grandChild1_0_4/";
    std::string cluster1_0_4_meanScaleINPUTPATH = casePath + cluster1_0_4_Folder + "scalerLocal_Input_onlyMean.txt";
    std::string cluster1_0_4_stdScaleINPUTPATH = casePath + cluster1_0_4_Folder + "scalerLocal_Input_onlyStandardDeviation.txt";
    std::string cluster1_0_4_meanScaleOUTPUTPATH = casePath + cluster1_0_4_Folder + "scalerLocal_Label_onlyMean.txt";
    std::string cluster1_0_4_stdScaleOUTPUTPATH = casePath + cluster1_0_4_Folder + "scalerLocal_Label_onlyStandardDeviation.txt";
    vector<float> cluster1_0_4_meanScaleINPUT;
    vector<float> cluster1_0_4_stdScaleINPUT;
    vector<float> cluster1_0_4_meanScaleOUTPUT;
    vector<float> cluster1_0_4_stdScaleOUTPUT;
    readFromCommaDelimitedFile_Float(cluster1_0_4_meanScaleINPUTPATH, cluster1_0_4_meanScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster1_0_4_stdScaleINPUTPATH, cluster1_0_4_stdScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster1_0_4_meanScaleOUTPUTPATH, cluster1_0_4_meanScaleOUTPUT);
    readFromCommaDelimitedFile_Float(cluster1_0_4_stdScaleOUTPUTPATH, cluster1_0_4_stdScaleOUTPUT);
    
    std::string cluster1_0_5_Folder = "ANNReg_parentCluster1/ANNReg_childCluster1_0_40grandChild/grandChild1_0_5/";
    std::string cluster1_0_5_meanScaleINPUTPATH = casePath + cluster1_0_5_Folder + "scalerLocal_Input_onlyMean.txt";
    std::string cluster1_0_5_stdScaleINPUTPATH = casePath + cluster1_0_5_Folder + "scalerLocal_Input_onlyStandardDeviation.txt";
    std::string cluster1_0_5_meanScaleOUTPUTPATH = casePath + cluster1_0_5_Folder + "scalerLocal_Label_onlyMean.txt";
    std::string cluster1_0_5_stdScaleOUTPUTPATH = casePath + cluster1_0_5_Folder + "scalerLocal_Label_onlyStandardDeviation.txt";
    vector<float> cluster1_0_5_meanScaleINPUT;
    vector<float> cluster1_0_5_stdScaleINPUT;
    vector<float> cluster1_0_5_meanScaleOUTPUT;
    vector<float> cluster1_0_5_stdScaleOUTPUT;
    readFromCommaDelimitedFile_Float(cluster1_0_5_meanScaleINPUTPATH, cluster1_0_5_meanScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster1_0_5_stdScaleINPUTPATH, cluster1_0_5_stdScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster1_0_5_meanScaleOUTPUTPATH, cluster1_0_5_meanScaleOUTPUT);
    readFromCommaDelimitedFile_Float(cluster1_0_5_stdScaleOUTPUTPATH, cluster1_0_5_stdScaleOUTPUT);
    
    std::string cluster1_0_6_Folder = "ANNReg_parentCluster1/ANNReg_childCluster1_0_40grandChild/grandChild1_0_6/";
    std::string cluster1_0_6_meanScaleINPUTPATH = casePath + cluster1_0_6_Folder + "scalerLocal_Input_onlyMean.txt";
    std::string cluster1_0_6_stdScaleINPUTPATH = casePath + cluster1_0_6_Folder + "scalerLocal_Input_onlyStandardDeviation.txt";
    std::string cluster1_0_6_meanScaleOUTPUTPATH = casePath + cluster1_0_6_Folder + "scalerLocal_Label_onlyMean.txt";
    std::string cluster1_0_6_stdScaleOUTPUTPATH = casePath + cluster1_0_6_Folder + "scalerLocal_Label_onlyStandardDeviation.txt";
    vector<float> cluster1_0_6_meanScaleINPUT;
    vector<float> cluster1_0_6_stdScaleINPUT;
    vector<float> cluster1_0_6_meanScaleOUTPUT;
    vector<float> cluster1_0_6_stdScaleOUTPUT;
    readFromCommaDelimitedFile_Float(cluster1_0_6_meanScaleINPUTPATH, cluster1_0_6_meanScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster1_0_6_stdScaleINPUTPATH, cluster1_0_6_stdScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster1_0_6_meanScaleOUTPUTPATH, cluster1_0_6_meanScaleOUTPUT);
    readFromCommaDelimitedFile_Float(cluster1_0_6_stdScaleOUTPUTPATH, cluster1_0_6_stdScaleOUTPUT);
    
    std::string cluster1_0_7_Folder = "ANNReg_parentCluster1/ANNReg_childCluster1_0_40grandChild/grandChild1_0_7/";
    std::string cluster1_0_7_meanScaleINPUTPATH = casePath + cluster1_0_7_Folder + "scalerLocal_Input_onlyMean.txt";
    std::string cluster1_0_7_stdScaleINPUTPATH = casePath + cluster1_0_7_Folder + "scalerLocal_Input_onlyStandardDeviation.txt";
    std::string cluster1_0_7_meanScaleOUTPUTPATH = casePath + cluster1_0_7_Folder + "scalerLocal_Label_onlyMean.txt";
    std::string cluster1_0_7_stdScaleOUTPUTPATH = casePath + cluster1_0_7_Folder + "scalerLocal_Label_onlyStandardDeviation.txt";
    vector<float> cluster1_0_7_meanScaleINPUT;
    vector<float> cluster1_0_7_stdScaleINPUT;
    vector<float> cluster1_0_7_meanScaleOUTPUT;
    vector<float> cluster1_0_7_stdScaleOUTPUT;
    readFromCommaDelimitedFile_Float(cluster1_0_7_meanScaleINPUTPATH, cluster1_0_7_meanScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster1_0_7_stdScaleINPUTPATH, cluster1_0_7_stdScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster1_0_7_meanScaleOUTPUTPATH, cluster1_0_7_meanScaleOUTPUT);
    readFromCommaDelimitedFile_Float(cluster1_0_7_stdScaleOUTPUTPATH, cluster1_0_7_stdScaleOUTPUT);
    
    std::string cluster1_0_8_Folder = "ANNReg_parentCluster1/ANNReg_childCluster1_0_40grandChild/grandChild1_0_8/";
    std::string cluster1_0_8_meanScaleINPUTPATH = casePath + cluster1_0_8_Folder + "scalerLocal_Input_onlyMean.txt";
    std::string cluster1_0_8_stdScaleINPUTPATH = casePath + cluster1_0_8_Folder + "scalerLocal_Input_onlyStandardDeviation.txt";
    std::string cluster1_0_8_meanScaleOUTPUTPATH = casePath + cluster1_0_8_Folder + "scalerLocal_Label_onlyMean.txt";
    std::string cluster1_0_8_stdScaleOUTPUTPATH = casePath + cluster1_0_8_Folder + "scalerLocal_Label_onlyStandardDeviation.txt";
    vector<float> cluster1_0_8_meanScaleINPUT;
    vector<float> cluster1_0_8_stdScaleINPUT;
    vector<float> cluster1_0_8_meanScaleOUTPUT;
    vector<float> cluster1_0_8_stdScaleOUTPUT;
    readFromCommaDelimitedFile_Float(cluster1_0_8_meanScaleINPUTPATH, cluster1_0_8_meanScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster1_0_8_stdScaleINPUTPATH, cluster1_0_8_stdScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster1_0_8_meanScaleOUTPUTPATH, cluster1_0_8_meanScaleOUTPUT);
    readFromCommaDelimitedFile_Float(cluster1_0_8_stdScaleOUTPUTPATH, cluster1_0_8_stdScaleOUTPUT);
    
    std::string cluster1_0_9_Folder = "ANNReg_parentCluster1/ANNReg_childCluster1_0_40grandChild/grandChild1_0_9/";
    std::string cluster1_0_9_meanScaleINPUTPATH = casePath + cluster1_0_9_Folder + "scalerLocal_Input_onlyMean.txt";
    std::string cluster1_0_9_stdScaleINPUTPATH = casePath + cluster1_0_9_Folder + "scalerLocal_Input_onlyStandardDeviation.txt";
    std::string cluster1_0_9_meanScaleOUTPUTPATH = casePath + cluster1_0_9_Folder + "scalerLocal_Label_onlyMean.txt";
    std::string cluster1_0_9_stdScaleOUTPUTPATH = casePath + cluster1_0_9_Folder + "scalerLocal_Label_onlyStandardDeviation.txt";
    vector<float> cluster1_0_9_meanScaleINPUT;
    vector<float> cluster1_0_9_stdScaleINPUT;
    vector<float> cluster1_0_9_meanScaleOUTPUT;
    vector<float> cluster1_0_9_stdScaleOUTPUT;
    readFromCommaDelimitedFile_Float(cluster1_0_9_meanScaleINPUTPATH, cluster1_0_9_meanScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster1_0_9_stdScaleINPUTPATH, cluster1_0_9_stdScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster1_0_9_meanScaleOUTPUTPATH, cluster1_0_9_meanScaleOUTPUT);
    readFromCommaDelimitedFile_Float(cluster1_0_9_stdScaleOUTPUTPATH, cluster1_0_9_stdScaleOUTPUT);
    
    std::string cluster1_0_10_Folder = "ANNReg_parentCluster1/ANNReg_childCluster1_0_40grandChild/grandChild1_0_10/";
    std::string cluster1_0_10_meanScaleINPUTPATH = casePath + cluster1_0_10_Folder + "scalerLocal_Input_onlyMean.txt";
    std::string cluster1_0_10_stdScaleINPUTPATH = casePath + cluster1_0_10_Folder + "scalerLocal_Input_onlyStandardDeviation.txt";
    std::string cluster1_0_10_meanScaleOUTPUTPATH = casePath + cluster1_0_10_Folder + "scalerLocal_Label_onlyMean.txt";
    std::string cluster1_0_10_stdScaleOUTPUTPATH = casePath + cluster1_0_10_Folder + "scalerLocal_Label_onlyStandardDeviation.txt";
    vector<float> cluster1_0_10_meanScaleINPUT;
    vector<float> cluster1_0_10_stdScaleINPUT;
    vector<float> cluster1_0_10_meanScaleOUTPUT;
    vector<float> cluster1_0_10_stdScaleOUTPUT;
    readFromCommaDelimitedFile_Float(cluster1_0_10_meanScaleINPUTPATH, cluster1_0_10_meanScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster1_0_10_stdScaleINPUTPATH, cluster1_0_10_stdScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster1_0_10_meanScaleOUTPUTPATH, cluster1_0_10_meanScaleOUTPUT);
    readFromCommaDelimitedFile_Float(cluster1_0_10_stdScaleOUTPUTPATH, cluster1_0_10_stdScaleOUTPUT);
    
    std::string cluster1_0_11_Folder = "ANNReg_parentCluster1/ANNReg_childCluster1_0_40grandChild/grandChild1_0_11/";
    std::string cluster1_0_11_meanScaleINPUTPATH = casePath + cluster1_0_11_Folder + "scalerLocal_Input_onlyMean.txt";
    std::string cluster1_0_11_stdScaleINPUTPATH = casePath + cluster1_0_11_Folder + "scalerLocal_Input_onlyStandardDeviation.txt";
    std::string cluster1_0_11_meanScaleOUTPUTPATH = casePath + cluster1_0_11_Folder + "scalerLocal_Label_onlyMean.txt";
    std::string cluster1_0_11_stdScaleOUTPUTPATH = casePath + cluster1_0_11_Folder + "scalerLocal_Label_onlyStandardDeviation.txt";
    vector<float> cluster1_0_11_meanScaleINPUT;
    vector<float> cluster1_0_11_stdScaleINPUT;
    vector<float> cluster1_0_11_meanScaleOUTPUT;
    vector<float> cluster1_0_11_stdScaleOUTPUT;
    readFromCommaDelimitedFile_Float(cluster1_0_11_meanScaleINPUTPATH, cluster1_0_11_meanScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster1_0_11_stdScaleINPUTPATH, cluster1_0_11_stdScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster1_0_11_meanScaleOUTPUTPATH, cluster1_0_11_meanScaleOUTPUT);
    readFromCommaDelimitedFile_Float(cluster1_0_11_stdScaleOUTPUTPATH, cluster1_0_11_stdScaleOUTPUT);
    
    std::string cluster1_0_12_Folder = "ANNReg_parentCluster1/ANNReg_childCluster1_0_40grandChild/grandChild1_0_12/";
    std::string cluster1_0_12_meanScaleINPUTPATH = casePath + cluster1_0_12_Folder + "scalerLocal_Input_onlyMean.txt";
    std::string cluster1_0_12_stdScaleINPUTPATH = casePath + cluster1_0_12_Folder + "scalerLocal_Input_onlyStandardDeviation.txt";
    std::string cluster1_0_12_meanScaleOUTPUTPATH = casePath + cluster1_0_12_Folder + "scalerLocal_Label_onlyMean.txt";
    std::string cluster1_0_12_stdScaleOUTPUTPATH = casePath + cluster1_0_12_Folder + "scalerLocal_Label_onlyStandardDeviation.txt";
    vector<float> cluster1_0_12_meanScaleINPUT;
    vector<float> cluster1_0_12_stdScaleINPUT;
    vector<float> cluster1_0_12_meanScaleOUTPUT;
    vector<float> cluster1_0_12_stdScaleOUTPUT;
    readFromCommaDelimitedFile_Float(cluster1_0_12_meanScaleINPUTPATH, cluster1_0_12_meanScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster1_0_12_stdScaleINPUTPATH, cluster1_0_12_stdScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster1_0_12_meanScaleOUTPUTPATH, cluster1_0_12_meanScaleOUTPUT);
    readFromCommaDelimitedFile_Float(cluster1_0_12_stdScaleOUTPUTPATH, cluster1_0_12_stdScaleOUTPUT);
    
    std::string cluster1_0_13_Folder = "ANNReg_parentCluster1/ANNReg_childCluster1_0_40grandChild/grandChild1_0_13/";
    std::string cluster1_0_13_meanScaleINPUTPATH = casePath + cluster1_0_13_Folder + "scalerLocal_Input_onlyMean.txt";
    std::string cluster1_0_13_stdScaleINPUTPATH = casePath + cluster1_0_13_Folder + "scalerLocal_Input_onlyStandardDeviation.txt";
    std::string cluster1_0_13_meanScaleOUTPUTPATH = casePath + cluster1_0_13_Folder + "scalerLocal_Label_onlyMean.txt";
    std::string cluster1_0_13_stdScaleOUTPUTPATH = casePath + cluster1_0_13_Folder + "scalerLocal_Label_onlyStandardDeviation.txt";
    vector<float> cluster1_0_13_meanScaleINPUT;
    vector<float> cluster1_0_13_stdScaleINPUT;
    vector<float> cluster1_0_13_meanScaleOUTPUT;
    vector<float> cluster1_0_13_stdScaleOUTPUT;
    readFromCommaDelimitedFile_Float(cluster1_0_13_meanScaleINPUTPATH, cluster1_0_13_meanScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster1_0_13_stdScaleINPUTPATH, cluster1_0_13_stdScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster1_0_13_meanScaleOUTPUTPATH, cluster1_0_13_meanScaleOUTPUT);
    readFromCommaDelimitedFile_Float(cluster1_0_13_stdScaleOUTPUTPATH, cluster1_0_13_stdScaleOUTPUT);
    
    std::string cluster1_0_14_Folder = "ANNReg_parentCluster1/ANNReg_childCluster1_0_40grandChild/grandChild1_0_14/";
    std::string cluster1_0_14_meanScaleINPUTPATH = casePath + cluster1_0_14_Folder + "scalerLocal_Input_onlyMean.txt";
    std::string cluster1_0_14_stdScaleINPUTPATH = casePath + cluster1_0_14_Folder + "scalerLocal_Input_onlyStandardDeviation.txt";
    std::string cluster1_0_14_meanScaleOUTPUTPATH = casePath + cluster1_0_14_Folder + "scalerLocal_Label_onlyMean.txt";
    std::string cluster1_0_14_stdScaleOUTPUTPATH = casePath + cluster1_0_14_Folder + "scalerLocal_Label_onlyStandardDeviation.txt";
    vector<float> cluster1_0_14_meanScaleINPUT;
    vector<float> cluster1_0_14_stdScaleINPUT;
    vector<float> cluster1_0_14_meanScaleOUTPUT;
    vector<float> cluster1_0_14_stdScaleOUTPUT;
    readFromCommaDelimitedFile_Float(cluster1_0_14_meanScaleINPUTPATH, cluster1_0_14_meanScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster1_0_14_stdScaleINPUTPATH, cluster1_0_14_stdScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster1_0_14_meanScaleOUTPUTPATH, cluster1_0_14_meanScaleOUTPUT);
    readFromCommaDelimitedFile_Float(cluster1_0_14_stdScaleOUTPUTPATH, cluster1_0_14_stdScaleOUTPUT);
    
    std::string cluster1_0_15_Folder = "ANNReg_parentCluster1/ANNReg_childCluster1_0_40grandChild/grandChild1_0_15/";
    std::string cluster1_0_15_meanScaleINPUTPATH = casePath + cluster1_0_15_Folder + "scalerLocal_Input_onlyMean.txt";
    std::string cluster1_0_15_stdScaleINPUTPATH = casePath + cluster1_0_15_Folder + "scalerLocal_Input_onlyStandardDeviation.txt";
    std::string cluster1_0_15_meanScaleOUTPUTPATH = casePath + cluster1_0_15_Folder + "scalerLocal_Label_onlyMean.txt";
    std::string cluster1_0_15_stdScaleOUTPUTPATH = casePath + cluster1_0_15_Folder + "scalerLocal_Label_onlyStandardDeviation.txt";
    vector<float> cluster1_0_15_meanScaleINPUT;
    vector<float> cluster1_0_15_stdScaleINPUT;
    vector<float> cluster1_0_15_meanScaleOUTPUT;
    vector<float> cluster1_0_15_stdScaleOUTPUT;
    readFromCommaDelimitedFile_Float(cluster1_0_15_meanScaleINPUTPATH, cluster1_0_15_meanScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster1_0_15_stdScaleINPUTPATH, cluster1_0_15_stdScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster1_0_15_meanScaleOUTPUTPATH, cluster1_0_15_meanScaleOUTPUT);
    readFromCommaDelimitedFile_Float(cluster1_0_15_stdScaleOUTPUTPATH, cluster1_0_15_stdScaleOUTPUT);
    
    std::string cluster1_0_16_Folder = "ANNReg_parentCluster1/ANNReg_childCluster1_0_40grandChild/grandChild1_0_16/";
    std::string cluster1_0_16_meanScaleINPUTPATH = casePath + cluster1_0_16_Folder + "scalerLocal_Input_onlyMean.txt";
    std::string cluster1_0_16_stdScaleINPUTPATH = casePath + cluster1_0_16_Folder + "scalerLocal_Input_onlyStandardDeviation.txt";
    std::string cluster1_0_16_meanScaleOUTPUTPATH = casePath + cluster1_0_16_Folder + "scalerLocal_Label_onlyMean.txt";
    std::string cluster1_0_16_stdScaleOUTPUTPATH = casePath + cluster1_0_16_Folder + "scalerLocal_Label_onlyStandardDeviation.txt";
    vector<float> cluster1_0_16_meanScaleINPUT;
    vector<float> cluster1_0_16_stdScaleINPUT;
    vector<float> cluster1_0_16_meanScaleOUTPUT;
    vector<float> cluster1_0_16_stdScaleOUTPUT;
    readFromCommaDelimitedFile_Float(cluster1_0_16_meanScaleINPUTPATH, cluster1_0_16_meanScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster1_0_16_stdScaleINPUTPATH, cluster1_0_16_stdScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster1_0_16_meanScaleOUTPUTPATH, cluster1_0_16_meanScaleOUTPUT);
    readFromCommaDelimitedFile_Float(cluster1_0_16_stdScaleOUTPUTPATH, cluster1_0_16_stdScaleOUTPUT);
    
    std::string cluster1_0_17_Folder = "ANNReg_parentCluster1/ANNReg_childCluster1_0_40grandChild/grandChild1_0_17/";
    std::string cluster1_0_17_meanScaleINPUTPATH = casePath + cluster1_0_17_Folder + "scalerLocal_Input_onlyMean.txt";
    std::string cluster1_0_17_stdScaleINPUTPATH = casePath + cluster1_0_17_Folder + "scalerLocal_Input_onlyStandardDeviation.txt";
    std::string cluster1_0_17_meanScaleOUTPUTPATH = casePath + cluster1_0_17_Folder + "scalerLocal_Label_onlyMean.txt";
    std::string cluster1_0_17_stdScaleOUTPUTPATH = casePath + cluster1_0_17_Folder + "scalerLocal_Label_onlyStandardDeviation.txt";
    vector<float> cluster1_0_17_meanScaleINPUT;
    vector<float> cluster1_0_17_stdScaleINPUT;
    vector<float> cluster1_0_17_meanScaleOUTPUT;
    vector<float> cluster1_0_17_stdScaleOUTPUT;
    readFromCommaDelimitedFile_Float(cluster1_0_17_meanScaleINPUTPATH, cluster1_0_17_meanScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster1_0_17_stdScaleINPUTPATH, cluster1_0_17_stdScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster1_0_17_meanScaleOUTPUTPATH, cluster1_0_17_meanScaleOUTPUT);
    readFromCommaDelimitedFile_Float(cluster1_0_17_stdScaleOUTPUTPATH, cluster1_0_17_stdScaleOUTPUT);
    
    std::string cluster1_0_18_Folder = "ANNReg_parentCluster1/ANNReg_childCluster1_0_40grandChild/grandChild1_0_18/";
    std::string cluster1_0_18_meanScaleINPUTPATH = casePath + cluster1_0_18_Folder + "scalerLocal_Input_onlyMean.txt";
    std::string cluster1_0_18_stdScaleINPUTPATH = casePath + cluster1_0_18_Folder + "scalerLocal_Input_onlyStandardDeviation.txt";
    std::string cluster1_0_18_meanScaleOUTPUTPATH = casePath + cluster1_0_18_Folder + "scalerLocal_Label_onlyMean.txt";
    std::string cluster1_0_18_stdScaleOUTPUTPATH = casePath + cluster1_0_18_Folder + "scalerLocal_Label_onlyStandardDeviation.txt";
    vector<float> cluster1_0_18_meanScaleINPUT;
    vector<float> cluster1_0_18_stdScaleINPUT;
    vector<float> cluster1_0_18_meanScaleOUTPUT;
    vector<float> cluster1_0_18_stdScaleOUTPUT;
    readFromCommaDelimitedFile_Float(cluster1_0_18_meanScaleINPUTPATH, cluster1_0_18_meanScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster1_0_18_stdScaleINPUTPATH, cluster1_0_18_stdScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster1_0_18_meanScaleOUTPUTPATH, cluster1_0_18_meanScaleOUTPUT);
    readFromCommaDelimitedFile_Float(cluster1_0_18_stdScaleOUTPUTPATH, cluster1_0_18_stdScaleOUTPUT);
    
    std::string cluster1_0_19_Folder = "ANNReg_parentCluster1/ANNReg_childCluster1_0_40grandChild/grandChild1_0_19/";
    std::string cluster1_0_19_meanScaleINPUTPATH = casePath + cluster1_0_19_Folder + "scalerLocal_Input_onlyMean.txt";
    std::string cluster1_0_19_stdScaleINPUTPATH = casePath + cluster1_0_19_Folder + "scalerLocal_Input_onlyStandardDeviation.txt";
    std::string cluster1_0_19_meanScaleOUTPUTPATH = casePath + cluster1_0_19_Folder + "scalerLocal_Label_onlyMean.txt";
    std::string cluster1_0_19_stdScaleOUTPUTPATH = casePath + cluster1_0_19_Folder + "scalerLocal_Label_onlyStandardDeviation.txt";
    vector<float> cluster1_0_19_meanScaleINPUT;
    vector<float> cluster1_0_19_stdScaleINPUT;
    vector<float> cluster1_0_19_meanScaleOUTPUT;
    vector<float> cluster1_0_19_stdScaleOUTPUT;
    readFromCommaDelimitedFile_Float(cluster1_0_19_meanScaleINPUTPATH, cluster1_0_19_meanScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster1_0_19_stdScaleINPUTPATH, cluster1_0_19_stdScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster1_0_19_meanScaleOUTPUTPATH, cluster1_0_19_meanScaleOUTPUT);
    readFromCommaDelimitedFile_Float(cluster1_0_19_stdScaleOUTPUTPATH, cluster1_0_19_stdScaleOUTPUT);
    

        // parentCluster1 - childCluster1
    std::string cluster1_1_0_Folder = "ANNReg_parentCluster1/ANNReg_childCluster1_1_15grandChild/grandChild1_1_0/";
    std::string cluster1_1_0_meanScaleINPUTPATH = casePath + cluster1_1_0_Folder + "scalerLocal_Input_onlyMean.txt";
    std::string cluster1_1_0_stdScaleINPUTPATH = casePath + cluster1_1_0_Folder + "scalerLocal_Input_onlyStandardDeviation.txt";
    std::string cluster1_1_0_meanScaleOUTPUTPATH = casePath + cluster1_1_0_Folder + "scalerLocal_Label_onlyMean.txt";
    std::string cluster1_1_0_stdScaleOUTPUTPATH = casePath + cluster1_1_0_Folder + "scalerLocal_Label_onlyStandardDeviation.txt";
    vector<float> cluster1_1_0_meanScaleINPUT;
    vector<float> cluster1_1_0_stdScaleINPUT;
    vector<float> cluster1_1_0_meanScaleOUTPUT;
    vector<float> cluster1_1_0_stdScaleOUTPUT;
    readFromCommaDelimitedFile_Float(cluster1_1_0_meanScaleINPUTPATH, cluster1_1_0_meanScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster1_1_0_stdScaleINPUTPATH, cluster1_1_0_stdScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster1_1_0_meanScaleOUTPUTPATH, cluster1_1_0_meanScaleOUTPUT);
    readFromCommaDelimitedFile_Float(cluster1_1_0_stdScaleOUTPUTPATH, cluster1_1_0_stdScaleOUTPUT);
    
    std::string cluster1_1_1_Folder = "ANNReg_parentCluster1/ANNReg_childCluster1_1_15grandChild/grandChild1_1_1/";
    std::string cluster1_1_1_meanScaleINPUTPATH = casePath + cluster1_1_1_Folder + "scalerLocal_Input_onlyMean.txt";
    std::string cluster1_1_1_stdScaleINPUTPATH = casePath + cluster1_1_1_Folder + "scalerLocal_Input_onlyStandardDeviation.txt";
    std::string cluster1_1_1_meanScaleOUTPUTPATH = casePath + cluster1_1_1_Folder + "scalerLocal_Label_onlyMean.txt";
    std::string cluster1_1_1_stdScaleOUTPUTPATH = casePath + cluster1_1_1_Folder + "scalerLocal_Label_onlyStandardDeviation.txt";
    vector<float> cluster1_1_1_meanScaleINPUT;
    vector<float> cluster1_1_1_stdScaleINPUT;
    vector<float> cluster1_1_1_meanScaleOUTPUT;
    vector<float> cluster1_1_1_stdScaleOUTPUT;
    readFromCommaDelimitedFile_Float(cluster1_1_1_meanScaleINPUTPATH, cluster1_1_1_meanScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster1_1_1_stdScaleINPUTPATH, cluster1_1_1_stdScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster1_1_1_meanScaleOUTPUTPATH, cluster1_1_1_meanScaleOUTPUT);
    readFromCommaDelimitedFile_Float(cluster1_1_1_stdScaleOUTPUTPATH, cluster1_1_1_stdScaleOUTPUT);
    
    std::string cluster1_1_2_Folder = "ANNReg_parentCluster1/ANNReg_childCluster1_1_15grandChild/grandChild1_1_2/";
    std::string cluster1_1_2_meanScaleINPUTPATH = casePath + cluster1_1_2_Folder + "scalerLocal_Input_onlyMean.txt";
    std::string cluster1_1_2_stdScaleINPUTPATH = casePath + cluster1_1_2_Folder + "scalerLocal_Input_onlyStandardDeviation.txt";
    std::string cluster1_1_2_meanScaleOUTPUTPATH = casePath + cluster1_1_2_Folder + "scalerLocal_Label_onlyMean.txt";
    std::string cluster1_1_2_stdScaleOUTPUTPATH = casePath + cluster1_1_2_Folder + "scalerLocal_Label_onlyStandardDeviation.txt";
    vector<float> cluster1_1_2_meanScaleINPUT;
    vector<float> cluster1_1_2_stdScaleINPUT;
    vector<float> cluster1_1_2_meanScaleOUTPUT;
    vector<float> cluster1_1_2_stdScaleOUTPUT;
    readFromCommaDelimitedFile_Float(cluster1_1_2_meanScaleINPUTPATH, cluster1_1_2_meanScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster1_1_2_stdScaleINPUTPATH, cluster1_1_2_stdScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster1_1_2_meanScaleOUTPUTPATH, cluster1_1_2_meanScaleOUTPUT);
    readFromCommaDelimitedFile_Float(cluster1_1_2_stdScaleOUTPUTPATH, cluster1_1_2_stdScaleOUTPUT);
    
    std::string cluster1_1_3_Folder = "ANNReg_parentCluster1/ANNReg_childCluster1_1_15grandChild/grandChild1_1_3/";
    std::string cluster1_1_3_meanScaleINPUTPATH = casePath + cluster1_1_3_Folder + "scalerLocal_Input_onlyMean.txt";
    std::string cluster1_1_3_stdScaleINPUTPATH = casePath + cluster1_1_3_Folder + "scalerLocal_Input_onlyStandardDeviation.txt";
    std::string cluster1_1_3_meanScaleOUTPUTPATH = casePath + cluster1_1_3_Folder + "scalerLocal_Label_onlyMean.txt";
    std::string cluster1_1_3_stdScaleOUTPUTPATH = casePath + cluster1_1_3_Folder + "scalerLocal_Label_onlyStandardDeviation.txt";
    vector<float> cluster1_1_3_meanScaleINPUT;
    vector<float> cluster1_1_3_stdScaleINPUT;
    vector<float> cluster1_1_3_meanScaleOUTPUT;
    vector<float> cluster1_1_3_stdScaleOUTPUT;
    readFromCommaDelimitedFile_Float(cluster1_1_3_meanScaleINPUTPATH, cluster1_1_3_meanScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster1_1_3_stdScaleINPUTPATH, cluster1_1_3_stdScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster1_1_3_meanScaleOUTPUTPATH, cluster1_1_3_meanScaleOUTPUT);
    readFromCommaDelimitedFile_Float(cluster1_1_3_stdScaleOUTPUTPATH, cluster1_1_3_stdScaleOUTPUT);
    
    std::string cluster1_1_4_Folder = "ANNReg_parentCluster1/ANNReg_childCluster1_1_15grandChild/grandChild1_1_4/";
    std::string cluster1_1_4_meanScaleINPUTPATH = casePath + cluster1_1_4_Folder + "scalerLocal_Input_onlyMean.txt";
    std::string cluster1_1_4_stdScaleINPUTPATH = casePath + cluster1_1_4_Folder + "scalerLocal_Input_onlyStandardDeviation.txt";
    std::string cluster1_1_4_meanScaleOUTPUTPATH = casePath + cluster1_1_4_Folder + "scalerLocal_Label_onlyMean.txt";
    std::string cluster1_1_4_stdScaleOUTPUTPATH = casePath + cluster1_1_4_Folder + "scalerLocal_Label_onlyStandardDeviation.txt";
    vector<float> cluster1_1_4_meanScaleINPUT;
    vector<float> cluster1_1_4_stdScaleINPUT;
    vector<float> cluster1_1_4_meanScaleOUTPUT;
    vector<float> cluster1_1_4_stdScaleOUTPUT;
    readFromCommaDelimitedFile_Float(cluster1_1_4_meanScaleINPUTPATH, cluster1_1_4_meanScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster1_1_4_stdScaleINPUTPATH, cluster1_1_4_stdScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster1_1_4_meanScaleOUTPUTPATH, cluster1_1_4_meanScaleOUTPUT);
    readFromCommaDelimitedFile_Float(cluster1_1_4_stdScaleOUTPUTPATH, cluster1_1_4_stdScaleOUTPUT);
    
    std::string cluster1_1_5_Folder = "ANNReg_parentCluster1/ANNReg_childCluster1_1_15grandChild/grandChild1_1_5/";
    std::string cluster1_1_5_meanScaleINPUTPATH = casePath + cluster1_1_5_Folder + "scalerLocal_Input_onlyMean.txt";
    std::string cluster1_1_5_stdScaleINPUTPATH = casePath + cluster1_1_5_Folder + "scalerLocal_Input_onlyStandardDeviation.txt";
    std::string cluster1_1_5_meanScaleOUTPUTPATH = casePath + cluster1_1_5_Folder + "scalerLocal_Label_onlyMean.txt";
    std::string cluster1_1_5_stdScaleOUTPUTPATH = casePath + cluster1_1_5_Folder + "scalerLocal_Label_onlyStandardDeviation.txt";
    vector<float> cluster1_1_5_meanScaleINPUT;
    vector<float> cluster1_1_5_stdScaleINPUT;
    vector<float> cluster1_1_5_meanScaleOUTPUT;
    vector<float> cluster1_1_5_stdScaleOUTPUT;
    readFromCommaDelimitedFile_Float(cluster1_1_5_meanScaleINPUTPATH, cluster1_1_5_meanScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster1_1_5_stdScaleINPUTPATH, cluster1_1_5_stdScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster1_1_5_meanScaleOUTPUTPATH, cluster1_1_5_meanScaleOUTPUT);
    readFromCommaDelimitedFile_Float(cluster1_1_5_stdScaleOUTPUTPATH, cluster1_1_5_stdScaleOUTPUT);
    
    std::string cluster1_1_6_Folder = "ANNReg_parentCluster1/ANNReg_childCluster1_1_15grandChild/grandChild1_1_6/";
    std::string cluster1_1_6_meanScaleINPUTPATH = casePath + cluster1_1_6_Folder + "scalerLocal_Input_onlyMean.txt";
    std::string cluster1_1_6_stdScaleINPUTPATH = casePath + cluster1_1_6_Folder + "scalerLocal_Input_onlyStandardDeviation.txt";
    std::string cluster1_1_6_meanScaleOUTPUTPATH = casePath + cluster1_1_6_Folder + "scalerLocal_Label_onlyMean.txt";
    std::string cluster1_1_6_stdScaleOUTPUTPATH = casePath + cluster1_1_6_Folder + "scalerLocal_Label_onlyStandardDeviation.txt";
    vector<float> cluster1_1_6_meanScaleINPUT;
    vector<float> cluster1_1_6_stdScaleINPUT;
    vector<float> cluster1_1_6_meanScaleOUTPUT;
    vector<float> cluster1_1_6_stdScaleOUTPUT;
    readFromCommaDelimitedFile_Float(cluster1_1_6_meanScaleINPUTPATH, cluster1_1_6_meanScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster1_1_6_stdScaleINPUTPATH, cluster1_1_6_stdScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster1_1_6_meanScaleOUTPUTPATH, cluster1_1_6_meanScaleOUTPUT);
    readFromCommaDelimitedFile_Float(cluster1_1_6_stdScaleOUTPUTPATH, cluster1_1_6_stdScaleOUTPUT);
    
    std::string cluster1_1_7_Folder = "ANNReg_parentCluster1/ANNReg_childCluster1_1_15grandChild/grandChild1_1_7/";
    std::string cluster1_1_7_meanScaleINPUTPATH = casePath + cluster1_1_7_Folder + "scalerLocal_Input_onlyMean.txt";
    std::string cluster1_1_7_stdScaleINPUTPATH = casePath + cluster1_1_7_Folder + "scalerLocal_Input_onlyStandardDeviation.txt";
    std::string cluster1_1_7_meanScaleOUTPUTPATH = casePath + cluster1_1_7_Folder + "scalerLocal_Label_onlyMean.txt";
    std::string cluster1_1_7_stdScaleOUTPUTPATH = casePath + cluster1_1_7_Folder + "scalerLocal_Label_onlyStandardDeviation.txt";
    vector<float> cluster1_1_7_meanScaleINPUT;
    vector<float> cluster1_1_7_stdScaleINPUT;
    vector<float> cluster1_1_7_meanScaleOUTPUT;
    vector<float> cluster1_1_7_stdScaleOUTPUT;
    readFromCommaDelimitedFile_Float(cluster1_1_7_meanScaleINPUTPATH, cluster1_1_7_meanScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster1_1_7_stdScaleINPUTPATH, cluster1_1_7_stdScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster1_1_7_meanScaleOUTPUTPATH, cluster1_1_7_meanScaleOUTPUT);
    readFromCommaDelimitedFile_Float(cluster1_1_7_stdScaleOUTPUTPATH, cluster1_1_7_stdScaleOUTPUT);
    
    std::string cluster1_1_8_Folder = "ANNReg_parentCluster1/ANNReg_childCluster1_1_15grandChild/grandChild1_1_8/";
    std::string cluster1_1_8_meanScaleINPUTPATH = casePath + cluster1_1_8_Folder + "scalerLocal_Input_onlyMean.txt";
    std::string cluster1_1_8_stdScaleINPUTPATH = casePath + cluster1_1_8_Folder + "scalerLocal_Input_onlyStandardDeviation.txt";
    std::string cluster1_1_8_meanScaleOUTPUTPATH = casePath + cluster1_1_8_Folder + "scalerLocal_Label_onlyMean.txt";
    std::string cluster1_1_8_stdScaleOUTPUTPATH = casePath + cluster1_1_8_Folder + "scalerLocal_Label_onlyStandardDeviation.txt";
    vector<float> cluster1_1_8_meanScaleINPUT;
    vector<float> cluster1_1_8_stdScaleINPUT;
    vector<float> cluster1_1_8_meanScaleOUTPUT;
    vector<float> cluster1_1_8_stdScaleOUTPUT;
    readFromCommaDelimitedFile_Float(cluster1_1_8_meanScaleINPUTPATH, cluster1_1_8_meanScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster1_1_8_stdScaleINPUTPATH, cluster1_1_8_stdScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster1_1_8_meanScaleOUTPUTPATH, cluster1_1_8_meanScaleOUTPUT);
    readFromCommaDelimitedFile_Float(cluster1_1_8_stdScaleOUTPUTPATH, cluster1_1_8_stdScaleOUTPUT);
    
    std::string cluster1_1_9_Folder = "ANNReg_parentCluster1/ANNReg_childCluster1_1_15grandChild/grandChild1_1_9/";
    std::string cluster1_1_9_meanScaleINPUTPATH = casePath + cluster1_1_9_Folder + "scalerLocal_Input_onlyMean.txt";
    std::string cluster1_1_9_stdScaleINPUTPATH = casePath + cluster1_1_9_Folder + "scalerLocal_Input_onlyStandardDeviation.txt";
    std::string cluster1_1_9_meanScaleOUTPUTPATH = casePath + cluster1_1_9_Folder + "scalerLocal_Label_onlyMean.txt";
    std::string cluster1_1_9_stdScaleOUTPUTPATH = casePath + cluster1_1_9_Folder + "scalerLocal_Label_onlyStandardDeviation.txt";
    vector<float> cluster1_1_9_meanScaleINPUT;
    vector<float> cluster1_1_9_stdScaleINPUT;
    vector<float> cluster1_1_9_meanScaleOUTPUT;
    vector<float> cluster1_1_9_stdScaleOUTPUT;
    readFromCommaDelimitedFile_Float(cluster1_1_9_meanScaleINPUTPATH, cluster1_1_9_meanScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster1_1_9_stdScaleINPUTPATH, cluster1_1_9_stdScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster1_1_9_meanScaleOUTPUTPATH, cluster1_1_9_meanScaleOUTPUT);
    readFromCommaDelimitedFile_Float(cluster1_1_9_stdScaleOUTPUTPATH, cluster1_1_9_stdScaleOUTPUT);
    

        // parentCluster1 - childCluster2
    std::string cluster1_2_0_Folder = "ANNReg_parentCluster1/ANNReg_childCluster1_2_15grandChild/grandChild1_2_0/";
    std::string cluster1_2_0_meanScaleINPUTPATH = casePath + cluster1_2_0_Folder + "scalerLocal_Input_onlyMean.txt";
    std::string cluster1_2_0_stdScaleINPUTPATH = casePath + cluster1_2_0_Folder + "scalerLocal_Input_onlyStandardDeviation.txt";
    std::string cluster1_2_0_meanScaleOUTPUTPATH = casePath + cluster1_2_0_Folder + "scalerLocal_Label_onlyMean.txt";
    std::string cluster1_2_0_stdScaleOUTPUTPATH = casePath + cluster1_2_0_Folder + "scalerLocal_Label_onlyStandardDeviation.txt";
    vector<float> cluster1_2_0_meanScaleINPUT;
    vector<float> cluster1_2_0_stdScaleINPUT;
    vector<float> cluster1_2_0_meanScaleOUTPUT;
    vector<float> cluster1_2_0_stdScaleOUTPUT;
    readFromCommaDelimitedFile_Float(cluster1_2_0_meanScaleINPUTPATH, cluster1_2_0_meanScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster1_2_0_stdScaleINPUTPATH, cluster1_2_0_stdScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster1_2_0_meanScaleOUTPUTPATH, cluster1_2_0_meanScaleOUTPUT);
    readFromCommaDelimitedFile_Float(cluster1_2_0_stdScaleOUTPUTPATH, cluster1_2_0_stdScaleOUTPUT);
    
    std::string cluster1_2_1_Folder = "ANNReg_parentCluster1/ANNReg_childCluster1_2_15grandChild/grandChild1_2_1/";
    std::string cluster1_2_1_meanScaleINPUTPATH = casePath + cluster1_2_1_Folder + "scalerLocal_Input_onlyMean.txt";
    std::string cluster1_2_1_stdScaleINPUTPATH = casePath + cluster1_2_1_Folder + "scalerLocal_Input_onlyStandardDeviation.txt";
    std::string cluster1_2_1_meanScaleOUTPUTPATH = casePath + cluster1_2_1_Folder + "scalerLocal_Label_onlyMean.txt";
    std::string cluster1_2_1_stdScaleOUTPUTPATH = casePath + cluster1_2_1_Folder + "scalerLocal_Label_onlyStandardDeviation.txt";
    vector<float> cluster1_2_1_meanScaleINPUT;
    vector<float> cluster1_2_1_stdScaleINPUT;
    vector<float> cluster1_2_1_meanScaleOUTPUT;
    vector<float> cluster1_2_1_stdScaleOUTPUT;
    readFromCommaDelimitedFile_Float(cluster1_2_1_meanScaleINPUTPATH, cluster1_2_1_meanScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster1_2_1_stdScaleINPUTPATH, cluster1_2_1_stdScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster1_2_1_meanScaleOUTPUTPATH, cluster1_2_1_meanScaleOUTPUT);
    readFromCommaDelimitedFile_Float(cluster1_2_1_stdScaleOUTPUTPATH, cluster1_2_1_stdScaleOUTPUT);
    
    std::string cluster1_2_2_Folder = "ANNReg_parentCluster1/ANNReg_childCluster1_2_15grandChild/grandChild1_2_2/";
    std::string cluster1_2_2_meanScaleINPUTPATH = casePath + cluster1_2_2_Folder + "scalerLocal_Input_onlyMean.txt";
    std::string cluster1_2_2_stdScaleINPUTPATH = casePath + cluster1_2_2_Folder + "scalerLocal_Input_onlyStandardDeviation.txt";
    std::string cluster1_2_2_meanScaleOUTPUTPATH = casePath + cluster1_2_2_Folder + "scalerLocal_Label_onlyMean.txt";
    std::string cluster1_2_2_stdScaleOUTPUTPATH = casePath + cluster1_2_2_Folder + "scalerLocal_Label_onlyStandardDeviation.txt";
    vector<float> cluster1_2_2_meanScaleINPUT;
    vector<float> cluster1_2_2_stdScaleINPUT;
    vector<float> cluster1_2_2_meanScaleOUTPUT;
    vector<float> cluster1_2_2_stdScaleOUTPUT;
    readFromCommaDelimitedFile_Float(cluster1_2_2_meanScaleINPUTPATH, cluster1_2_2_meanScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster1_2_2_stdScaleINPUTPATH, cluster1_2_2_stdScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster1_2_2_meanScaleOUTPUTPATH, cluster1_2_2_meanScaleOUTPUT);
    readFromCommaDelimitedFile_Float(cluster1_2_2_stdScaleOUTPUTPATH, cluster1_2_2_stdScaleOUTPUT);
    
    std::string cluster1_2_3_Folder = "ANNReg_parentCluster1/ANNReg_childCluster1_2_15grandChild/grandChild1_2_3/";
    std::string cluster1_2_3_meanScaleINPUTPATH = casePath + cluster1_2_3_Folder + "scalerLocal_Input_onlyMean.txt";
    std::string cluster1_2_3_stdScaleINPUTPATH = casePath + cluster1_2_3_Folder + "scalerLocal_Input_onlyStandardDeviation.txt";
    std::string cluster1_2_3_meanScaleOUTPUTPATH = casePath + cluster1_2_3_Folder + "scalerLocal_Label_onlyMean.txt";
    std::string cluster1_2_3_stdScaleOUTPUTPATH = casePath + cluster1_2_3_Folder + "scalerLocal_Label_onlyStandardDeviation.txt";
    vector<float> cluster1_2_3_meanScaleINPUT;
    vector<float> cluster1_2_3_stdScaleINPUT;
    vector<float> cluster1_2_3_meanScaleOUTPUT;
    vector<float> cluster1_2_3_stdScaleOUTPUT;
    readFromCommaDelimitedFile_Float(cluster1_2_3_meanScaleINPUTPATH, cluster1_2_3_meanScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster1_2_3_stdScaleINPUTPATH, cluster1_2_3_stdScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster1_2_3_meanScaleOUTPUTPATH, cluster1_2_3_meanScaleOUTPUT);
    readFromCommaDelimitedFile_Float(cluster1_2_3_stdScaleOUTPUTPATH, cluster1_2_3_stdScaleOUTPUT);
    
    std::string cluster1_2_4_Folder = "ANNReg_parentCluster1/ANNReg_childCluster1_2_15grandChild/grandChild1_2_4/";
    std::string cluster1_2_4_meanScaleINPUTPATH = casePath + cluster1_2_4_Folder + "scalerLocal_Input_onlyMean.txt";
    std::string cluster1_2_4_stdScaleINPUTPATH = casePath + cluster1_2_4_Folder + "scalerLocal_Input_onlyStandardDeviation.txt";
    std::string cluster1_2_4_meanScaleOUTPUTPATH = casePath + cluster1_2_4_Folder + "scalerLocal_Label_onlyMean.txt";
    std::string cluster1_2_4_stdScaleOUTPUTPATH = casePath + cluster1_2_4_Folder + "scalerLocal_Label_onlyStandardDeviation.txt";
    vector<float> cluster1_2_4_meanScaleINPUT;
    vector<float> cluster1_2_4_stdScaleINPUT;
    vector<float> cluster1_2_4_meanScaleOUTPUT;
    vector<float> cluster1_2_4_stdScaleOUTPUT;
    readFromCommaDelimitedFile_Float(cluster1_2_4_meanScaleINPUTPATH, cluster1_2_4_meanScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster1_2_4_stdScaleINPUTPATH, cluster1_2_4_stdScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster1_2_4_meanScaleOUTPUTPATH, cluster1_2_4_meanScaleOUTPUT);
    readFromCommaDelimitedFile_Float(cluster1_2_4_stdScaleOUTPUTPATH, cluster1_2_4_stdScaleOUTPUT);
    
    std::string cluster1_2_5_Folder = "ANNReg_parentCluster1/ANNReg_childCluster1_2_15grandChild/grandChild1_2_5/";
    std::string cluster1_2_5_meanScaleINPUTPATH = casePath + cluster1_2_5_Folder + "scalerLocal_Input_onlyMean.txt";
    std::string cluster1_2_5_stdScaleINPUTPATH = casePath + cluster1_2_5_Folder + "scalerLocal_Input_onlyStandardDeviation.txt";
    std::string cluster1_2_5_meanScaleOUTPUTPATH = casePath + cluster1_2_5_Folder + "scalerLocal_Label_onlyMean.txt";
    std::string cluster1_2_5_stdScaleOUTPUTPATH = casePath + cluster1_2_5_Folder + "scalerLocal_Label_onlyStandardDeviation.txt";
    vector<float> cluster1_2_5_meanScaleINPUT;
    vector<float> cluster1_2_5_stdScaleINPUT;
    vector<float> cluster1_2_5_meanScaleOUTPUT;
    vector<float> cluster1_2_5_stdScaleOUTPUT;
    readFromCommaDelimitedFile_Float(cluster1_2_5_meanScaleINPUTPATH, cluster1_2_5_meanScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster1_2_5_stdScaleINPUTPATH, cluster1_2_5_stdScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster1_2_5_meanScaleOUTPUTPATH, cluster1_2_5_meanScaleOUTPUT);
    readFromCommaDelimitedFile_Float(cluster1_2_5_stdScaleOUTPUTPATH, cluster1_2_5_stdScaleOUTPUT);
    
    std::string cluster1_2_6_Folder = "ANNReg_parentCluster1/ANNReg_childCluster1_2_15grandChild/grandChild1_2_6/";
    std::string cluster1_2_6_meanScaleINPUTPATH = casePath + cluster1_2_6_Folder + "scalerLocal_Input_onlyMean.txt";
    std::string cluster1_2_6_stdScaleINPUTPATH = casePath + cluster1_2_6_Folder + "scalerLocal_Input_onlyStandardDeviation.txt";
    std::string cluster1_2_6_meanScaleOUTPUTPATH = casePath + cluster1_2_6_Folder + "scalerLocal_Label_onlyMean.txt";
    std::string cluster1_2_6_stdScaleOUTPUTPATH = casePath + cluster1_2_6_Folder + "scalerLocal_Label_onlyStandardDeviation.txt";
    vector<float> cluster1_2_6_meanScaleINPUT;
    vector<float> cluster1_2_6_stdScaleINPUT;
    vector<float> cluster1_2_6_meanScaleOUTPUT;
    vector<float> cluster1_2_6_stdScaleOUTPUT;
    readFromCommaDelimitedFile_Float(cluster1_2_6_meanScaleINPUTPATH, cluster1_2_6_meanScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster1_2_6_stdScaleINPUTPATH, cluster1_2_6_stdScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster1_2_6_meanScaleOUTPUTPATH, cluster1_2_6_meanScaleOUTPUT);
    readFromCommaDelimitedFile_Float(cluster1_2_6_stdScaleOUTPUTPATH, cluster1_2_6_stdScaleOUTPUT);
    
    std::string cluster1_2_7_Folder = "ANNReg_parentCluster1/ANNReg_childCluster1_2_15grandChild/grandChild1_2_7/";
    std::string cluster1_2_7_meanScaleINPUTPATH = casePath + cluster1_2_7_Folder + "scalerLocal_Input_onlyMean.txt";
    std::string cluster1_2_7_stdScaleINPUTPATH = casePath + cluster1_2_7_Folder + "scalerLocal_Input_onlyStandardDeviation.txt";
    std::string cluster1_2_7_meanScaleOUTPUTPATH = casePath + cluster1_2_7_Folder + "scalerLocal_Label_onlyMean.txt";
    std::string cluster1_2_7_stdScaleOUTPUTPATH = casePath + cluster1_2_7_Folder + "scalerLocal_Label_onlyStandardDeviation.txt";
    vector<float> cluster1_2_7_meanScaleINPUT;
    vector<float> cluster1_2_7_stdScaleINPUT;
    vector<float> cluster1_2_7_meanScaleOUTPUT;
    vector<float> cluster1_2_7_stdScaleOUTPUT;
    readFromCommaDelimitedFile_Float(cluster1_2_7_meanScaleINPUTPATH, cluster1_2_7_meanScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster1_2_7_stdScaleINPUTPATH, cluster1_2_7_stdScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster1_2_7_meanScaleOUTPUTPATH, cluster1_2_7_meanScaleOUTPUT);
    readFromCommaDelimitedFile_Float(cluster1_2_7_stdScaleOUTPUTPATH, cluster1_2_7_stdScaleOUTPUT);
    
    std::string cluster1_2_8_Folder = "ANNReg_parentCluster1/ANNReg_childCluster1_2_15grandChild/grandChild1_2_8/";
    std::string cluster1_2_8_meanScaleINPUTPATH = casePath + cluster1_2_8_Folder + "scalerLocal_Input_onlyMean.txt";
    std::string cluster1_2_8_stdScaleINPUTPATH = casePath + cluster1_2_8_Folder + "scalerLocal_Input_onlyStandardDeviation.txt";
    std::string cluster1_2_8_meanScaleOUTPUTPATH = casePath + cluster1_2_8_Folder + "scalerLocal_Label_onlyMean.txt";
    std::string cluster1_2_8_stdScaleOUTPUTPATH = casePath + cluster1_2_8_Folder + "scalerLocal_Label_onlyStandardDeviation.txt";
    vector<float> cluster1_2_8_meanScaleINPUT;
    vector<float> cluster1_2_8_stdScaleINPUT;
    vector<float> cluster1_2_8_meanScaleOUTPUT;
    vector<float> cluster1_2_8_stdScaleOUTPUT;
    readFromCommaDelimitedFile_Float(cluster1_2_8_meanScaleINPUTPATH, cluster1_2_8_meanScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster1_2_8_stdScaleINPUTPATH, cluster1_2_8_stdScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster1_2_8_meanScaleOUTPUTPATH, cluster1_2_8_meanScaleOUTPUT);
    readFromCommaDelimitedFile_Float(cluster1_2_8_stdScaleOUTPUTPATH, cluster1_2_8_stdScaleOUTPUT);
    
    std::string cluster1_2_9_Folder = "ANNReg_parentCluster1/ANNReg_childCluster1_2_15grandChild/grandChild1_2_9/";
    std::string cluster1_2_9_meanScaleINPUTPATH = casePath + cluster1_2_9_Folder + "scalerLocal_Input_onlyMean.txt";
    std::string cluster1_2_9_stdScaleINPUTPATH = casePath + cluster1_2_9_Folder + "scalerLocal_Input_onlyStandardDeviation.txt";
    std::string cluster1_2_9_meanScaleOUTPUTPATH = casePath + cluster1_2_9_Folder + "scalerLocal_Label_onlyMean.txt";
    std::string cluster1_2_9_stdScaleOUTPUTPATH = casePath + cluster1_2_9_Folder + "scalerLocal_Label_onlyStandardDeviation.txt";
    vector<float> cluster1_2_9_meanScaleINPUT;
    vector<float> cluster1_2_9_stdScaleINPUT;
    vector<float> cluster1_2_9_meanScaleOUTPUT;
    vector<float> cluster1_2_9_stdScaleOUTPUT;
    readFromCommaDelimitedFile_Float(cluster1_2_9_meanScaleINPUTPATH, cluster1_2_9_meanScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster1_2_9_stdScaleINPUTPATH, cluster1_2_9_stdScaleINPUT);
    readFromCommaDelimitedFile_Float(cluster1_2_9_meanScaleOUTPUTPATH, cluster1_2_9_meanScaleOUTPUT);
    readFromCommaDelimitedFile_Float(cluster1_2_9_stdScaleOUTPUTPATH, cluster1_2_9_stdScaleOUTPUT);
    

    // Declare Max,Min of OUTPUT (MaxMinVar - deltaY deltaT) standardized to drive the prediction
    // Child of parentCluster0
    std::string cluster0_0_maxOUTPUTPATH = casePath + cluster0_0_Folder + "onlyMax_Var.txt";
    std::string cluster0_0_minOUTPUTPATH = casePath + cluster0_0_Folder + "onlyMin_Var.txt";
    vector<float> cluster0_0_maxOUTPUT;
    vector<float> cluster0_0_minOUTPUT;
    readFromCommaDelimitedFile_Float(cluster0_0_maxOUTPUTPATH, cluster0_0_maxOUTPUT);
    readFromCommaDelimitedFile_Float(cluster0_0_minOUTPUTPATH, cluster0_0_minOUTPUT);
    
    std::string cluster0_1_maxOUTPUTPATH = casePath + cluster0_1_Folder + "onlyMax_Var.txt";
    std::string cluster0_1_minOUTPUTPATH = casePath + cluster0_1_Folder + "onlyMin_Var.txt";
    vector<float> cluster0_1_maxOUTPUT;
    vector<float> cluster0_1_minOUTPUT;
    readFromCommaDelimitedFile_Float(cluster0_1_maxOUTPUTPATH, cluster0_1_maxOUTPUT);
    readFromCommaDelimitedFile_Float(cluster0_1_minOUTPUTPATH, cluster0_1_minOUTPUT);
    
    std::string cluster0_2_maxOUTPUTPATH = casePath + cluster0_2_Folder + "onlyMax_Var.txt";
    std::string cluster0_2_minOUTPUTPATH = casePath + cluster0_2_Folder + "onlyMin_Var.txt";
    vector<float> cluster0_2_maxOUTPUT;
    vector<float> cluster0_2_minOUTPUT;
    readFromCommaDelimitedFile_Float(cluster0_2_maxOUTPUTPATH, cluster0_2_maxOUTPUT);
    readFromCommaDelimitedFile_Float(cluster0_2_minOUTPUTPATH, cluster0_2_minOUTPUT);
 
    std::string cluster0_3_maxOUTPUTPATH = casePath + cluster0_3_Folder + "onlyMax_Var.txt";
    std::string cluster0_3_minOUTPUTPATH = casePath + cluster0_3_Folder + "onlyMin_Var.txt";
    vector<float> cluster0_3_maxOUTPUT;
    vector<float> cluster0_3_minOUTPUT;
    readFromCommaDelimitedFile_Float(cluster0_3_maxOUTPUTPATH, cluster0_3_maxOUTPUT);
    readFromCommaDelimitedFile_Float(cluster0_3_minOUTPUTPATH, cluster0_3_minOUTPUT);
    
    std::string cluster0_4_maxOUTPUTPATH = casePath + cluster0_4_Folder + "onlyMax_Var.txt";
    std::string cluster0_4_minOUTPUTPATH = casePath + cluster0_4_Folder + "onlyMin_Var.txt";
    vector<float> cluster0_4_maxOUTPUT;
    vector<float> cluster0_4_minOUTPUT;
    readFromCommaDelimitedFile_Float(cluster0_4_maxOUTPUTPATH, cluster0_4_maxOUTPUT);
    readFromCommaDelimitedFile_Float(cluster0_4_minOUTPUTPATH, cluster0_4_minOUTPUT);
    
    std::string cluster0_5_maxOUTPUTPATH = casePath + cluster0_5_Folder + "onlyMax_Var.txt";
    std::string cluster0_5_minOUTPUTPATH = casePath + cluster0_5_Folder + "onlyMin_Var.txt";
    vector<float> cluster0_5_maxOUTPUT;
    vector<float> cluster0_5_minOUTPUT;
    readFromCommaDelimitedFile_Float(cluster0_5_maxOUTPUTPATH, cluster0_5_maxOUTPUT);
    readFromCommaDelimitedFile_Float(cluster0_5_minOUTPUTPATH, cluster0_5_minOUTPUT);
    
    std::string cluster0_6_maxOUTPUTPATH = casePath + cluster0_6_Folder + "onlyMax_Var.txt";
    std::string cluster0_6_minOUTPUTPATH = casePath + cluster0_6_Folder + "onlyMin_Var.txt";
    vector<float> cluster0_6_maxOUTPUT;
    vector<float> cluster0_6_minOUTPUT;
    readFromCommaDelimitedFile_Float(cluster0_6_maxOUTPUTPATH, cluster0_6_maxOUTPUT);
    readFromCommaDelimitedFile_Float(cluster0_6_minOUTPUTPATH, cluster0_6_minOUTPUT);
    
    std::string cluster0_7_maxOUTPUTPATH = casePath + cluster0_7_Folder + "onlyMax_Var.txt";
    std::string cluster0_7_minOUTPUTPATH = casePath + cluster0_7_Folder + "onlyMin_Var.txt";
    vector<float> cluster0_7_maxOUTPUT;
    vector<float> cluster0_7_minOUTPUT;
    readFromCommaDelimitedFile_Float(cluster0_7_maxOUTPUTPATH, cluster0_7_maxOUTPUT);
    readFromCommaDelimitedFile_Float(cluster0_7_minOUTPUTPATH, cluster0_7_minOUTPUT);
    
    std::string cluster0_8_maxOUTPUTPATH = casePath + cluster0_8_Folder + "onlyMax_Var.txt";
    std::string cluster0_8_minOUTPUTPATH = casePath + cluster0_8_Folder + "onlyMin_Var.txt";
    vector<float> cluster0_8_maxOUTPUT;
    vector<float> cluster0_8_minOUTPUT;
    readFromCommaDelimitedFile_Float(cluster0_8_maxOUTPUTPATH, cluster0_8_maxOUTPUT);
    readFromCommaDelimitedFile_Float(cluster0_8_minOUTPUTPATH, cluster0_8_minOUTPUT);
    
    std::string cluster0_9_maxOUTPUTPATH = casePath + cluster0_9_Folder + "onlyMax_Var.txt";
    std::string cluster0_9_minOUTPUTPATH = casePath + cluster0_9_Folder + "onlyMin_Var.txt";
    vector<float> cluster0_9_maxOUTPUT;
    vector<float> cluster0_9_minOUTPUT;
    readFromCommaDelimitedFile_Float(cluster0_9_maxOUTPUTPATH, cluster0_9_maxOUTPUT);
    readFromCommaDelimitedFile_Float(cluster0_9_minOUTPUTPATH, cluster0_9_minOUTPUT);
    
    std::string cluster0_10_maxOUTPUTPATH = casePath + cluster0_10_Folder + "onlyMax_Var.txt";
    std::string cluster0_10_minOUTPUTPATH = casePath + cluster0_10_Folder + "onlyMin_Var.txt";
    vector<float> cluster0_10_maxOUTPUT;
    vector<float> cluster0_10_minOUTPUT;
    readFromCommaDelimitedFile_Float(cluster0_10_maxOUTPUTPATH, cluster0_10_maxOUTPUT);
    readFromCommaDelimitedFile_Float(cluster0_10_minOUTPUTPATH, cluster0_10_minOUTPUT);
    
    std::string cluster0_11_maxOUTPUTPATH = casePath + cluster0_11_Folder + "onlyMax_Var.txt";
    std::string cluster0_11_minOUTPUTPATH = casePath + cluster0_11_Folder + "onlyMin_Var.txt";
    vector<float> cluster0_11_maxOUTPUT;
    vector<float> cluster0_11_minOUTPUT;
    readFromCommaDelimitedFile_Float(cluster0_11_maxOUTPUTPATH, cluster0_11_maxOUTPUT);
    readFromCommaDelimitedFile_Float(cluster0_11_minOUTPUTPATH, cluster0_11_minOUTPUT);
    
    std::string cluster0_12_maxOUTPUTPATH = casePath + cluster0_12_Folder + "onlyMax_Var.txt";
    std::string cluster0_12_minOUTPUTPATH = casePath + cluster0_12_Folder + "onlyMin_Var.txt";
    vector<float> cluster0_12_maxOUTPUT;
    vector<float> cluster0_12_minOUTPUT;
    readFromCommaDelimitedFile_Float(cluster0_12_maxOUTPUTPATH, cluster0_12_maxOUTPUT);
    readFromCommaDelimitedFile_Float(cluster0_12_minOUTPUTPATH, cluster0_12_minOUTPUT);
    
    std::string cluster0_13_maxOUTPUTPATH = casePath + cluster0_13_Folder + "onlyMax_Var.txt";
    std::string cluster0_13_minOUTPUTPATH = casePath + cluster0_13_Folder + "onlyMin_Var.txt";
    vector<float> cluster0_13_maxOUTPUT;
    vector<float> cluster0_13_minOUTPUT;
    readFromCommaDelimitedFile_Float(cluster0_13_maxOUTPUTPATH, cluster0_13_maxOUTPUT);
    readFromCommaDelimitedFile_Float(cluster0_13_minOUTPUTPATH, cluster0_13_minOUTPUT);
    
    std::string cluster0_14_maxOUTPUTPATH = casePath + cluster0_14_Folder + "onlyMax_Var.txt";
    std::string cluster0_14_minOUTPUTPATH = casePath + cluster0_14_Folder + "onlyMin_Var.txt";
    vector<float> cluster0_14_maxOUTPUT;
    vector<float> cluster0_14_minOUTPUT;
    readFromCommaDelimitedFile_Float(cluster0_14_maxOUTPUTPATH, cluster0_14_maxOUTPUT);
    readFromCommaDelimitedFile_Float(cluster0_14_minOUTPUTPATH, cluster0_14_minOUTPUT);
    
    std::string cluster0_15_maxOUTPUTPATH = casePath + cluster0_15_Folder + "onlyMax_Var.txt";
    std::string cluster0_15_minOUTPUTPATH = casePath + cluster0_15_Folder + "onlyMin_Var.txt";
    vector<float> cluster0_15_maxOUTPUT;
    vector<float> cluster0_15_minOUTPUT;
    readFromCommaDelimitedFile_Float(cluster0_15_maxOUTPUTPATH, cluster0_15_maxOUTPUT);
    readFromCommaDelimitedFile_Float(cluster0_15_minOUTPUTPATH, cluster0_15_minOUTPUT);
    
    std::string cluster0_16_maxOUTPUTPATH = casePath + cluster0_16_Folder + "onlyMax_Var.txt";
    std::string cluster0_16_minOUTPUTPATH = casePath + cluster0_16_Folder + "onlyMin_Var.txt";
    vector<float> cluster0_16_maxOUTPUT;
    vector<float> cluster0_16_minOUTPUT;
    readFromCommaDelimitedFile_Float(cluster0_16_maxOUTPUTPATH, cluster0_16_maxOUTPUT);
    readFromCommaDelimitedFile_Float(cluster0_16_minOUTPUTPATH, cluster0_16_minOUTPUT);
    
    std::string cluster0_17_maxOUTPUTPATH = casePath + cluster0_17_Folder + "onlyMax_Var.txt";
    std::string cluster0_17_minOUTPUTPATH = casePath + cluster0_17_Folder + "onlyMin_Var.txt";
    vector<float> cluster0_17_maxOUTPUT;
    vector<float> cluster0_17_minOUTPUT;
    readFromCommaDelimitedFile_Float(cluster0_17_maxOUTPUTPATH, cluster0_17_maxOUTPUT);
    readFromCommaDelimitedFile_Float(cluster0_17_minOUTPUTPATH, cluster0_17_minOUTPUT);
    
    std::string cluster0_18_maxOUTPUTPATH = casePath + cluster0_18_Folder + "onlyMax_Var.txt";
    std::string cluster0_18_minOUTPUTPATH = casePath + cluster0_18_Folder + "onlyMin_Var.txt";
    vector<float> cluster0_18_maxOUTPUT;
    vector<float> cluster0_18_minOUTPUT;
    readFromCommaDelimitedFile_Float(cluster0_18_maxOUTPUTPATH, cluster0_18_maxOUTPUT);
    readFromCommaDelimitedFile_Float(cluster0_18_minOUTPUTPATH, cluster0_18_minOUTPUT);
    
    std::string cluster0_19_maxOUTPUTPATH = casePath + cluster0_19_Folder + "onlyMax_Var.txt";
    std::string cluster0_19_minOUTPUTPATH = casePath + cluster0_19_Folder + "onlyMin_Var.txt";
    vector<float> cluster0_19_maxOUTPUT;
    vector<float> cluster0_19_minOUTPUT;
    readFromCommaDelimitedFile_Float(cluster0_19_maxOUTPUTPATH, cluster0_19_maxOUTPUT);
    readFromCommaDelimitedFile_Float(cluster0_19_minOUTPUTPATH, cluster0_19_minOUTPUT);
    
    
    // Child of parentCluster1
    // parentCluster1 - chilCluster0 - grandChild //It should be maxVar, minVar
    std::string cluster1_0_0_maxOUTPUTPATH = casePath + cluster1_0_0_Folder + "onlyMax_Var.txt";
    std::string cluster1_0_0_minOUTPUTPATH = casePath + cluster1_0_0_Folder + "onlyMin_Var.txt";
    vector<float> cluster1_0_0_maxOUTPUT;
    vector<float> cluster1_0_0_minOUTPUT;
    readFromCommaDelimitedFile_Float(cluster1_0_0_maxOUTPUTPATH, cluster1_0_0_maxOUTPUT);
    readFromCommaDelimitedFile_Float(cluster1_0_0_minOUTPUTPATH, cluster1_0_0_minOUTPUT);
    
    std::string cluster1_0_1_maxOUTPUTPATH = casePath + cluster1_0_1_Folder + "onlyMax_Var.txt";
    std::string cluster1_0_1_minOUTPUTPATH = casePath + cluster1_0_1_Folder + "onlyMin_Var.txt";
    vector<float> cluster1_0_1_maxOUTPUT;
    vector<float> cluster1_0_1_minOUTPUT;
    readFromCommaDelimitedFile_Float(cluster1_0_1_maxOUTPUTPATH, cluster1_0_1_maxOUTPUT);
    readFromCommaDelimitedFile_Float(cluster1_0_1_minOUTPUTPATH, cluster1_0_1_minOUTPUT);
    
    std::string cluster1_0_2_maxOUTPUTPATH = casePath + cluster1_0_2_Folder + "onlyMax_Var.txt";
    std::string cluster1_0_2_minOUTPUTPATH = casePath + cluster1_0_2_Folder + "onlyMin_Var.txt";
    vector<float> cluster1_0_2_maxOUTPUT;
    vector<float> cluster1_0_2_minOUTPUT;
    readFromCommaDelimitedFile_Float(cluster1_0_2_maxOUTPUTPATH, cluster1_0_2_maxOUTPUT);
    readFromCommaDelimitedFile_Float(cluster1_0_2_minOUTPUTPATH, cluster1_0_2_minOUTPUT);
    
    std::string cluster1_0_3_maxOUTPUTPATH = casePath + cluster1_0_3_Folder + "onlyMax_Var.txt";
    std::string cluster1_0_3_minOUTPUTPATH = casePath + cluster1_0_3_Folder + "onlyMin_Var.txt";
    vector<float> cluster1_0_3_maxOUTPUT;
    vector<float> cluster1_0_3_minOUTPUT;
    readFromCommaDelimitedFile_Float(cluster1_0_3_maxOUTPUTPATH, cluster1_0_3_maxOUTPUT);
    readFromCommaDelimitedFile_Float(cluster1_0_3_minOUTPUTPATH, cluster1_0_3_minOUTPUT);
    
    std::string cluster1_0_4_maxOUTPUTPATH = casePath + cluster1_0_4_Folder + "onlyMax_Var.txt";
    std::string cluster1_0_4_minOUTPUTPATH = casePath + cluster1_0_4_Folder + "onlyMin_Var.txt";
    vector<float> cluster1_0_4_maxOUTPUT;
    vector<float> cluster1_0_4_minOUTPUT;
    readFromCommaDelimitedFile_Float(cluster1_0_4_maxOUTPUTPATH, cluster1_0_4_maxOUTPUT);
    readFromCommaDelimitedFile_Float(cluster1_0_4_minOUTPUTPATH, cluster1_0_4_minOUTPUT);
    
    std::string cluster1_0_5_maxOUTPUTPATH = casePath + cluster1_0_5_Folder + "onlyMax_Var.txt";
    std::string cluster1_0_5_minOUTPUTPATH = casePath + cluster1_0_5_Folder + "onlyMin_Var.txt";
    vector<float> cluster1_0_5_maxOUTPUT;
    vector<float> cluster1_0_5_minOUTPUT;
    readFromCommaDelimitedFile_Float(cluster1_0_5_maxOUTPUTPATH, cluster1_0_5_maxOUTPUT);
    readFromCommaDelimitedFile_Float(cluster1_0_5_minOUTPUTPATH, cluster1_0_5_minOUTPUT);
    
    std::string cluster1_0_6_maxOUTPUTPATH = casePath + cluster1_0_6_Folder + "onlyMax_Var.txt";
    std::string cluster1_0_6_minOUTPUTPATH = casePath + cluster1_0_6_Folder + "onlyMin_Var.txt";
    vector<float> cluster1_0_6_maxOUTPUT;
    vector<float> cluster1_0_6_minOUTPUT;
    readFromCommaDelimitedFile_Float(cluster1_0_6_maxOUTPUTPATH, cluster1_0_6_maxOUTPUT);
    readFromCommaDelimitedFile_Float(cluster1_0_6_minOUTPUTPATH, cluster1_0_6_minOUTPUT);
    
    std::string cluster1_0_7_maxOUTPUTPATH = casePath + cluster1_0_7_Folder + "onlyMax_Var.txt";
    std::string cluster1_0_7_minOUTPUTPATH = casePath + cluster1_0_7_Folder + "onlyMin_Var.txt";
    vector<float> cluster1_0_7_maxOUTPUT;
    vector<float> cluster1_0_7_minOUTPUT;
    readFromCommaDelimitedFile_Float(cluster1_0_7_maxOUTPUTPATH, cluster1_0_7_maxOUTPUT);
    readFromCommaDelimitedFile_Float(cluster1_0_7_minOUTPUTPATH, cluster1_0_7_minOUTPUT);
    
    std::string cluster1_0_8_maxOUTPUTPATH = casePath + cluster1_0_8_Folder + "onlyMax_Var.txt";
    std::string cluster1_0_8_minOUTPUTPATH = casePath + cluster1_0_8_Folder + "onlyMin_Var.txt";
    vector<float> cluster1_0_8_maxOUTPUT;
    vector<float> cluster1_0_8_minOUTPUT;
    readFromCommaDelimitedFile_Float(cluster1_0_8_maxOUTPUTPATH, cluster1_0_8_maxOUTPUT);
    readFromCommaDelimitedFile_Float(cluster1_0_8_minOUTPUTPATH, cluster1_0_8_minOUTPUT);
    
    std::string cluster1_0_9_maxOUTPUTPATH = casePath + cluster1_0_9_Folder + "onlyMax_Var.txt";
    std::string cluster1_0_9_minOUTPUTPATH = casePath + cluster1_0_9_Folder + "onlyMin_Var.txt";
    vector<float> cluster1_0_9_maxOUTPUT;
    vector<float> cluster1_0_9_minOUTPUT;
    readFromCommaDelimitedFile_Float(cluster1_0_9_maxOUTPUTPATH, cluster1_0_9_maxOUTPUT);
    readFromCommaDelimitedFile_Float(cluster1_0_9_minOUTPUTPATH, cluster1_0_9_minOUTPUT);
    
    std::string cluster1_0_10_maxOUTPUTPATH = casePath + cluster1_0_10_Folder + "onlyMax_Var.txt";
    std::string cluster1_0_10_minOUTPUTPATH = casePath + cluster1_0_10_Folder + "onlyMin_Var.txt";
    vector<float> cluster1_0_10_maxOUTPUT;
    vector<float> cluster1_0_10_minOUTPUT;
    readFromCommaDelimitedFile_Float(cluster1_0_10_maxOUTPUTPATH, cluster1_0_10_maxOUTPUT);
    readFromCommaDelimitedFile_Float(cluster1_0_10_minOUTPUTPATH, cluster1_0_10_minOUTPUT);
    
    std::string cluster1_0_11_maxOUTPUTPATH = casePath + cluster1_0_11_Folder + "onlyMax_Var.txt";
    std::string cluster1_0_11_minOUTPUTPATH = casePath + cluster1_0_11_Folder + "onlyMin_Var.txt";
    vector<float> cluster1_0_11_maxOUTPUT;
    vector<float> cluster1_0_11_minOUTPUT;
    readFromCommaDelimitedFile_Float(cluster1_0_11_maxOUTPUTPATH, cluster1_0_11_maxOUTPUT);
    readFromCommaDelimitedFile_Float(cluster1_0_11_minOUTPUTPATH, cluster1_0_11_minOUTPUT);
    
    std::string cluster1_0_12_maxOUTPUTPATH = casePath + cluster1_0_12_Folder + "onlyMax_Var.txt";
    std::string cluster1_0_12_minOUTPUTPATH = casePath + cluster1_0_12_Folder + "onlyMin_Var.txt";
    vector<float> cluster1_0_12_maxOUTPUT;
    vector<float> cluster1_0_12_minOUTPUT;
    readFromCommaDelimitedFile_Float(cluster1_0_12_maxOUTPUTPATH, cluster1_0_12_maxOUTPUT);
    readFromCommaDelimitedFile_Float(cluster1_0_12_minOUTPUTPATH, cluster1_0_12_minOUTPUT);
    
    std::string cluster1_0_13_maxOUTPUTPATH = casePath + cluster1_0_13_Folder + "onlyMax_Var.txt";
    std::string cluster1_0_13_minOUTPUTPATH = casePath + cluster1_0_13_Folder + "onlyMin_Var.txt";
    vector<float> cluster1_0_13_maxOUTPUT;
    vector<float> cluster1_0_13_minOUTPUT;
    readFromCommaDelimitedFile_Float(cluster1_0_13_maxOUTPUTPATH, cluster1_0_13_maxOUTPUT);
    readFromCommaDelimitedFile_Float(cluster1_0_13_minOUTPUTPATH, cluster1_0_13_minOUTPUT);
    
    std::string cluster1_0_14_maxOUTPUTPATH = casePath + cluster1_0_14_Folder + "onlyMax_Var.txt";
    std::string cluster1_0_14_minOUTPUTPATH = casePath + cluster1_0_14_Folder + "onlyMin_Var.txt";
    vector<float> cluster1_0_14_maxOUTPUT;
    vector<float> cluster1_0_14_minOUTPUT;
    readFromCommaDelimitedFile_Float(cluster1_0_14_maxOUTPUTPATH, cluster1_0_14_maxOUTPUT);
    readFromCommaDelimitedFile_Float(cluster1_0_14_minOUTPUTPATH, cluster1_0_14_minOUTPUT);
    
    std::string cluster1_0_15_maxOUTPUTPATH = casePath + cluster1_0_15_Folder + "onlyMax_Var.txt";
    std::string cluster1_0_15_minOUTPUTPATH = casePath + cluster1_0_15_Folder + "onlyMin_Var.txt";
    vector<float> cluster1_0_15_maxOUTPUT;
    vector<float> cluster1_0_15_minOUTPUT;
    readFromCommaDelimitedFile_Float(cluster1_0_15_maxOUTPUTPATH, cluster1_0_15_maxOUTPUT);
    readFromCommaDelimitedFile_Float(cluster1_0_15_minOUTPUTPATH, cluster1_0_15_minOUTPUT);
    
    std::string cluster1_0_16_maxOUTPUTPATH = casePath + cluster1_0_16_Folder + "onlyMax_Var.txt";
    std::string cluster1_0_16_minOUTPUTPATH = casePath + cluster1_0_16_Folder + "onlyMin_Var.txt";
    vector<float> cluster1_0_16_maxOUTPUT;
    vector<float> cluster1_0_16_minOUTPUT;
    readFromCommaDelimitedFile_Float(cluster1_0_16_maxOUTPUTPATH, cluster1_0_16_maxOUTPUT);
    readFromCommaDelimitedFile_Float(cluster1_0_16_minOUTPUTPATH, cluster1_0_16_minOUTPUT);
    
    std::string cluster1_0_17_maxOUTPUTPATH = casePath + cluster1_0_17_Folder + "onlyMax_Var.txt";
    std::string cluster1_0_17_minOUTPUTPATH = casePath + cluster1_0_17_Folder + "onlyMin_Var.txt";
    vector<float> cluster1_0_17_maxOUTPUT;
    vector<float> cluster1_0_17_minOUTPUT;
    readFromCommaDelimitedFile_Float(cluster1_0_17_maxOUTPUTPATH, cluster1_0_17_maxOUTPUT);
    readFromCommaDelimitedFile_Float(cluster1_0_17_minOUTPUTPATH, cluster1_0_17_minOUTPUT);
    
    std::string cluster1_0_18_maxOUTPUTPATH = casePath + cluster1_0_18_Folder + "onlyMax_Var.txt";
    std::string cluster1_0_18_minOUTPUTPATH = casePath + cluster1_0_18_Folder + "onlyMin_Var.txt";
    vector<float> cluster1_0_18_maxOUTPUT;
    vector<float> cluster1_0_18_minOUTPUT;
    readFromCommaDelimitedFile_Float(cluster1_0_18_maxOUTPUTPATH, cluster1_0_18_maxOUTPUT);
    readFromCommaDelimitedFile_Float(cluster1_0_18_minOUTPUTPATH, cluster1_0_18_minOUTPUT);
    
    std::string cluster1_0_19_maxOUTPUTPATH = casePath + cluster1_0_19_Folder + "onlyMax_Var.txt";
    std::string cluster1_0_19_minOUTPUTPATH = casePath + cluster1_0_19_Folder + "onlyMin_Var.txt";
    vector<float> cluster1_0_19_maxOUTPUT;
    vector<float> cluster1_0_19_minOUTPUT;
    readFromCommaDelimitedFile_Float(cluster1_0_19_maxOUTPUTPATH, cluster1_0_19_maxOUTPUT);
    readFromCommaDelimitedFile_Float(cluster1_0_19_minOUTPUTPATH, cluster1_0_19_minOUTPUT);
    
    
    // parentCluster1 - chilCluster1 - grandChild
    std::string cluster1_1_0_maxOUTPUTPATH = casePath + cluster1_1_0_Folder + "onlyMax_Var.txt";
    std::string cluster1_1_0_minOUTPUTPATH = casePath + cluster1_1_0_Folder + "onlyMin_Var.txt";
    vector<float> cluster1_1_0_maxOUTPUT;
    vector<float> cluster1_1_0_minOUTPUT;
    readFromCommaDelimitedFile_Float(cluster1_1_0_maxOUTPUTPATH, cluster1_1_0_maxOUTPUT);
    readFromCommaDelimitedFile_Float(cluster1_1_0_minOUTPUTPATH, cluster1_1_0_minOUTPUT);
    
    std::string cluster1_1_1_maxOUTPUTPATH = casePath + cluster1_1_1_Folder + "onlyMax_Var.txt";
    std::string cluster1_1_1_minOUTPUTPATH = casePath + cluster1_1_1_Folder + "onlyMin_Var.txt";
    vector<float> cluster1_1_1_maxOUTPUT;
    vector<float> cluster1_1_1_minOUTPUT;
    readFromCommaDelimitedFile_Float(cluster1_1_1_maxOUTPUTPATH, cluster1_1_1_maxOUTPUT);
    readFromCommaDelimitedFile_Float(cluster1_1_1_minOUTPUTPATH, cluster1_1_1_minOUTPUT);
    
    std::string cluster1_1_2_maxOUTPUTPATH = casePath + cluster1_1_2_Folder + "onlyMax_Var.txt";
    std::string cluster1_1_2_minOUTPUTPATH = casePath + cluster1_1_2_Folder + "onlyMin_Var.txt";
    vector<float> cluster1_1_2_maxOUTPUT;
    vector<float> cluster1_1_2_minOUTPUT;
    readFromCommaDelimitedFile_Float(cluster1_1_2_maxOUTPUTPATH, cluster1_1_2_maxOUTPUT);
    readFromCommaDelimitedFile_Float(cluster1_1_2_minOUTPUTPATH, cluster1_1_2_minOUTPUT);
    
    std::string cluster1_1_3_maxOUTPUTPATH = casePath + cluster1_1_3_Folder + "onlyMax_Var.txt";
    std::string cluster1_1_3_minOUTPUTPATH = casePath + cluster1_1_3_Folder + "onlyMin_Var.txt";
    vector<float> cluster1_1_3_maxOUTPUT;
    vector<float> cluster1_1_3_minOUTPUT;
    readFromCommaDelimitedFile_Float(cluster1_1_3_maxOUTPUTPATH, cluster1_1_3_maxOUTPUT);
    readFromCommaDelimitedFile_Float(cluster1_1_3_minOUTPUTPATH, cluster1_1_3_minOUTPUT);
    
    std::string cluster1_1_4_maxOUTPUTPATH = casePath + cluster1_1_4_Folder + "onlyMax_Var.txt";
    std::string cluster1_1_4_minOUTPUTPATH = casePath + cluster1_1_4_Folder + "onlyMin_Var.txt";
    vector<float> cluster1_1_4_maxOUTPUT;
    vector<float> cluster1_1_4_minOUTPUT;
    readFromCommaDelimitedFile_Float(cluster1_1_4_maxOUTPUTPATH, cluster1_1_4_maxOUTPUT);
    readFromCommaDelimitedFile_Float(cluster1_1_4_minOUTPUTPATH, cluster1_1_4_minOUTPUT);
    
    std::string cluster1_1_5_maxOUTPUTPATH = casePath + cluster1_1_5_Folder + "onlyMax_Var.txt";
    std::string cluster1_1_5_minOUTPUTPATH = casePath + cluster1_1_5_Folder + "onlyMin_Var.txt";
    vector<float> cluster1_1_5_maxOUTPUT;
    vector<float> cluster1_1_5_minOUTPUT;
    readFromCommaDelimitedFile_Float(cluster1_1_5_maxOUTPUTPATH, cluster1_1_5_maxOUTPUT);
    readFromCommaDelimitedFile_Float(cluster1_1_5_minOUTPUTPATH, cluster1_1_5_minOUTPUT);
    
    std::string cluster1_1_6_maxOUTPUTPATH = casePath + cluster1_1_6_Folder + "onlyMax_Var.txt";
    std::string cluster1_1_6_minOUTPUTPATH = casePath + cluster1_1_6_Folder + "onlyMin_Var.txt";
    vector<float> cluster1_1_6_maxOUTPUT;
    vector<float> cluster1_1_6_minOUTPUT;
    readFromCommaDelimitedFile_Float(cluster1_1_6_maxOUTPUTPATH, cluster1_1_6_maxOUTPUT);
    readFromCommaDelimitedFile_Float(cluster1_1_6_minOUTPUTPATH, cluster1_1_6_minOUTPUT);
    
    std::string cluster1_1_7_maxOUTPUTPATH = casePath + cluster1_1_7_Folder + "onlyMax_Var.txt";
    std::string cluster1_1_7_minOUTPUTPATH = casePath + cluster1_1_7_Folder + "onlyMin_Var.txt";
    vector<float> cluster1_1_7_maxOUTPUT;
    vector<float> cluster1_1_7_minOUTPUT;
    readFromCommaDelimitedFile_Float(cluster1_1_7_maxOUTPUTPATH, cluster1_1_7_maxOUTPUT);
    readFromCommaDelimitedFile_Float(cluster1_1_7_minOUTPUTPATH, cluster1_1_7_minOUTPUT);
    
    std::string cluster1_1_8_maxOUTPUTPATH = casePath + cluster1_1_8_Folder + "onlyMax_Var.txt";
    std::string cluster1_1_8_minOUTPUTPATH = casePath + cluster1_1_8_Folder + "onlyMin_Var.txt";
    vector<float> cluster1_1_8_maxOUTPUT;
    vector<float> cluster1_1_8_minOUTPUT;
    readFromCommaDelimitedFile_Float(cluster1_1_8_maxOUTPUTPATH, cluster1_1_8_maxOUTPUT);
    readFromCommaDelimitedFile_Float(cluster1_1_8_minOUTPUTPATH, cluster1_1_8_minOUTPUT);
    
    std::string cluster1_1_9_maxOUTPUTPATH = casePath + cluster1_1_9_Folder + "onlyMax_Var.txt";
    std::string cluster1_1_9_minOUTPUTPATH = casePath + cluster1_1_9_Folder + "onlyMin_Var.txt";
    vector<float> cluster1_1_9_maxOUTPUT;
    vector<float> cluster1_1_9_minOUTPUT;
    readFromCommaDelimitedFile_Float(cluster1_1_9_maxOUTPUTPATH, cluster1_1_9_maxOUTPUT);
    readFromCommaDelimitedFile_Float(cluster1_1_9_minOUTPUTPATH, cluster1_1_9_minOUTPUT);
    
    
    // parentCluster1 - chilCluster2 - grandChild
    std::string cluster1_2_0_maxOUTPUTPATH = casePath + cluster1_2_0_Folder + "onlyMax_Var.txt";
    std::string cluster1_2_0_minOUTPUTPATH = casePath + cluster1_2_0_Folder + "onlyMin_Var.txt";
    vector<float> cluster1_2_0_maxOUTPUT;
    vector<float> cluster1_2_0_minOUTPUT;
    readFromCommaDelimitedFile_Float(cluster1_2_0_maxOUTPUTPATH, cluster1_2_0_maxOUTPUT);
    readFromCommaDelimitedFile_Float(cluster1_2_0_minOUTPUTPATH, cluster1_2_0_minOUTPUT);
    
    std::string cluster1_2_1_maxOUTPUTPATH = casePath + cluster1_2_1_Folder + "onlyMax_Var.txt";
    std::string cluster1_2_1_minOUTPUTPATH = casePath + cluster1_2_1_Folder + "onlyMin_Var.txt";
    vector<float> cluster1_2_1_maxOUTPUT;
    vector<float> cluster1_2_1_minOUTPUT;
    readFromCommaDelimitedFile_Float(cluster1_2_1_maxOUTPUTPATH, cluster1_2_1_maxOUTPUT);
    readFromCommaDelimitedFile_Float(cluster1_2_1_minOUTPUTPATH, cluster1_2_1_minOUTPUT);
    
    std::string cluster1_2_2_maxOUTPUTPATH = casePath + cluster1_2_2_Folder + "onlyMax_Var.txt";
    std::string cluster1_2_2_minOUTPUTPATH = casePath + cluster1_2_2_Folder + "onlyMin_Var.txt";
    vector<float> cluster1_2_2_maxOUTPUT;
    vector<float> cluster1_2_2_minOUTPUT;
    readFromCommaDelimitedFile_Float(cluster1_2_2_maxOUTPUTPATH, cluster1_2_2_maxOUTPUT);
    readFromCommaDelimitedFile_Float(cluster1_2_2_minOUTPUTPATH, cluster1_2_2_minOUTPUT);
    
    std::string cluster1_2_3_maxOUTPUTPATH = casePath + cluster1_2_3_Folder + "onlyMax_Var.txt";
    std::string cluster1_2_3_minOUTPUTPATH = casePath + cluster1_2_3_Folder + "onlyMin_Var.txt";
    vector<float> cluster1_2_3_maxOUTPUT;
    vector<float> cluster1_2_3_minOUTPUT;
    readFromCommaDelimitedFile_Float(cluster1_2_3_maxOUTPUTPATH, cluster1_2_3_maxOUTPUT);
    readFromCommaDelimitedFile_Float(cluster1_2_3_minOUTPUTPATH, cluster1_2_3_minOUTPUT);
    
    std::string cluster1_2_4_maxOUTPUTPATH = casePath + cluster1_2_4_Folder + "onlyMax_Var.txt";
    std::string cluster1_2_4_minOUTPUTPATH = casePath + cluster1_2_4_Folder + "onlyMin_Var.txt";
    vector<float> cluster1_2_4_maxOUTPUT;
    vector<float> cluster1_2_4_minOUTPUT;
    readFromCommaDelimitedFile_Float(cluster1_2_4_maxOUTPUTPATH, cluster1_2_4_maxOUTPUT);
    readFromCommaDelimitedFile_Float(cluster1_2_4_minOUTPUTPATH, cluster1_2_4_minOUTPUT);
    
    std::string cluster1_2_5_maxOUTPUTPATH = casePath + cluster1_2_5_Folder + "onlyMax_Var.txt";
    std::string cluster1_2_5_minOUTPUTPATH = casePath + cluster1_2_5_Folder + "onlyMin_Var.txt";
    vector<float> cluster1_2_5_maxOUTPUT;
    vector<float> cluster1_2_5_minOUTPUT;
    readFromCommaDelimitedFile_Float(cluster1_2_5_maxOUTPUTPATH, cluster1_2_5_maxOUTPUT);
    readFromCommaDelimitedFile_Float(cluster1_2_5_minOUTPUTPATH, cluster1_2_5_minOUTPUT);
    
    std::string cluster1_2_6_maxOUTPUTPATH = casePath + cluster1_2_6_Folder + "onlyMax_Var.txt";
    std::string cluster1_2_6_minOUTPUTPATH = casePath + cluster1_2_6_Folder + "onlyMin_Var.txt";
    vector<float> cluster1_2_6_maxOUTPUT;
    vector<float> cluster1_2_6_minOUTPUT;
    readFromCommaDelimitedFile_Float(cluster1_2_6_maxOUTPUTPATH, cluster1_2_6_maxOUTPUT);
    readFromCommaDelimitedFile_Float(cluster1_2_6_minOUTPUTPATH, cluster1_2_6_minOUTPUT);
    
    std::string cluster1_2_7_maxOUTPUTPATH = casePath + cluster1_2_7_Folder + "onlyMax_Var.txt";
    std::string cluster1_2_7_minOUTPUTPATH = casePath + cluster1_2_7_Folder + "onlyMin_Var.txt";
    vector<float> cluster1_2_7_maxOUTPUT;
    vector<float> cluster1_2_7_minOUTPUT;
    readFromCommaDelimitedFile_Float(cluster1_2_7_maxOUTPUTPATH, cluster1_2_7_maxOUTPUT);
    readFromCommaDelimitedFile_Float(cluster1_2_7_minOUTPUTPATH, cluster1_2_7_minOUTPUT);
    
    std::string cluster1_2_8_maxOUTPUTPATH = casePath + cluster1_2_8_Folder + "onlyMax_Var.txt";
    std::string cluster1_2_8_minOUTPUTPATH = casePath + cluster1_2_8_Folder + "onlyMin_Var.txt";
    vector<float> cluster1_2_8_maxOUTPUT;
    vector<float> cluster1_2_8_minOUTPUT;
    readFromCommaDelimitedFile_Float(cluster1_2_8_maxOUTPUTPATH, cluster1_2_8_maxOUTPUT);
    readFromCommaDelimitedFile_Float(cluster1_2_8_minOUTPUTPATH, cluster1_2_8_minOUTPUT);
    
    std::string cluster1_2_9_maxOUTPUTPATH = casePath + cluster1_2_9_Folder + "onlyMax_Var.txt";
    std::string cluster1_2_9_minOUTPUTPATH = casePath + cluster1_2_9_Folder + "onlyMin_Var.txt";
    vector<float> cluster1_2_9_maxOUTPUT;
    vector<float> cluster1_2_9_minOUTPUT;
    readFromCommaDelimitedFile_Float(cluster1_2_9_maxOUTPUTPATH, cluster1_2_9_maxOUTPUT);
    readFromCommaDelimitedFile_Float(cluster1_2_9_minOUTPUTPATH, cluster1_2_9_minOUTPUT);
    

/* HuuTri@20211006: Commented - ANN without PCA
    // Declare PCA parameters for each cluster: mean, eigenVal, eigenVec then Load
    // Child of parentCluster0
    std::string cluster0_0_pcaMeanLoadPATH = casePath + cluster0_0_Folder + "PCA_mean.txt";
    std::string cluster0_0_pcaEigenValLoadPATH = casePath + cluster0_0_Folder + "PCA_eigenValues.txt";
    std::string cluster0_0_pcaEigenVecLoadPATH = casePath + cluster0_0_Folder + "PCA_eigenVectors.txt";
    vector<float> cluster0_0_pcaMeanLoad;
    vector<float> cluster0_0_pcaEigenValLoad;
    vector<float> cluster0_0_pcaEigenVecLoad;
    readFromCommaDelimitedFile_Float(cluster0_0_pcaMeanLoadPATH, cluster0_0_pcaMeanLoad);
    readFromCommaDelimitedFile_Float(cluster0_0_pcaEigenValLoadPATH, cluster0_0_pcaEigenValLoad);
    cv::PCA pcaANN0_0;
    pcaANN0_0.mean = cv::Mat(1, numVarANN, CV_32FC1, cluster0_0_pcaMeanLoad.data());
    pcaANN0_0.eigenvalues = cv::Mat(numComponentPCA,1, CV_32FC1, cluster0_0_pcaEigenValLoad.data());
    pcaANN0_0.eigenvectors = cv::Mat(numComponentPCA,numVarANN, CV_32FC1, cluster0_0_pcaEigenVecLoad.data());
    
    std::string cluster0_1_pcaMeanLoadPATH = casePath + cluster0_1_Folder + "PCA_mean.txt";
    std::string cluster0_1_pcaEigenValLoadPATH = casePath + cluster0_1_Folder + "PCA_eigenValues.txt";
    std::string cluster0_1_pcaEigenVecLoadPATH = casePath + cluster0_1_Folder + "PCA_eigenVectors.txt";
    vector<float> cluster0_1_pcaMeanLoad;
    vector<float> cluster0_1_pcaEigenValLoad;
    vector<float> cluster0_1_pcaEigenVecLoad;
    readFromCommaDelimitedFile_Float(cluster0_1_pcaMeanLoadPATH, cluster0_1_pcaMeanLoad);
    readFromCommaDelimitedFile_Float(cluster0_1_pcaEigenValLoadPATH, cluster0_1_pcaEigenValLoad);
    readFromCommaDelimitedFile_Float(cluster0_1_pcaEigenVecLoadPATH, cluster0_1_pcaEigenVecLoad);
    cv::PCA pcaANN0_1;
    pcaANN0_1.mean = cv::Mat(1, numVarANN, CV_32FC1, cluster0_1_pcaMeanLoad.data());
    pcaANN0_1.eigenvalues = cv::Mat(numComponentPCA,1, CV_32FC1, cluster0_1_pcaEigenValLoad.data());
    pcaANN0_1.eigenvectors = cv::Mat(numComponentPCA,numVarANN, CV_32FC1, cluster0_1_pcaEigenVecLoad.data());
    
    std::string cluster0_2_pcaMeanLoadPATH = casePath + cluster0_2_Folder + "PCA_mean.txt";
    std::string cluster0_2_pcaEigenValLoadPATH = casePath + cluster0_2_Folder + "PCA_eigenValues.txt";
    std::string cluster0_2_pcaEigenVecLoadPATH = casePath + cluster0_2_Folder + "PCA_eigenVectors.txt";
    vector<float> cluster0_2_pcaMeanLoad;
    vector<float> cluster0_2_pcaEigenValLoad;
    vector<float> cluster0_2_pcaEigenVecLoad;
    readFromCommaDelimitedFile_Float(cluster0_2_pcaMeanLoadPATH, cluster0_2_pcaMeanLoad);
    readFromCommaDelimitedFile_Float(cluster0_2_pcaEigenValLoadPATH, cluster0_2_pcaEigenValLoad);
    readFromCommaDelimitedFile_Float(cluster0_2_pcaEigenVecLoadPATH, cluster0_2_pcaEigenVecLoad);
    cv::PCA pcaANN0_2;
    pcaANN0_2.mean = cv::Mat(1, numVarANN, CV_32FC1, cluster0_2_pcaMeanLoad.data());
    pcaANN0_2.eigenvalues = cv::Mat(numComponentPCA,1, CV_32FC1, cluster0_2_pcaEigenValLoad.data());
    pcaANN0_2.eigenvectors = cv::Mat(numComponentPCA,numVarANN, CV_32FC1, cluster0_2_pcaEigenVecLoad.data());
    
    std::string cluster0_3_pcaMeanLoadPATH = casePath + cluster0_3_Folder + "PCA_mean.txt";
    std::string cluster0_3_pcaEigenValLoadPATH = casePath + cluster0_3_Folder + "PCA_eigenValues.txt";
    std::string cluster0_3_pcaEigenVecLoadPATH = casePath + cluster0_3_Folder + "PCA_eigenVectors.txt";
    vector<float> cluster0_3_pcaMeanLoad;
    vector<float> cluster0_3_pcaEigenValLoad;
    vector<float> cluster0_3_pcaEigenVecLoad;
    readFromCommaDelimitedFile_Float(cluster0_3_pcaMeanLoadPATH, cluster0_3_pcaMeanLoad);
    readFromCommaDelimitedFile_Float(cluster0_3_pcaEigenValLoadPATH, cluster0_3_pcaEigenValLoad);
    readFromCommaDelimitedFile_Float(cluster0_3_pcaEigenVecLoadPATH, cluster0_3_pcaEigenVecLoad);
    cv::PCA pcaANN0_3;
    pcaANN0_3.mean = cv::Mat(1, numVarANN, CV_32FC1, cluster0_3_pcaMeanLoad.data());
    pcaANN0_3.eigenvalues = cv::Mat(numComponentPCA,1, CV_32FC1, cluster0_3_pcaEigenValLoad.data());
    pcaANN0_3.eigenvectors = cv::Mat(numComponentPCA,numVarANN, CV_32FC1, cluster0_3_pcaEigenVecLoad.data());
    
    std::string cluster0_4_pcaMeanLoadPATH = casePath + cluster0_4_Folder + "PCA_mean.txt";
    std::string cluster0_4_pcaEigenValLoadPATH = casePath + cluster0_4_Folder + "PCA_eigenValues.txt";
    std::string cluster0_4_pcaEigenVecLoadPATH = casePath + cluster0_4_Folder + "PCA_eigenVectors.txt";
    vector<float> cluster0_4_pcaMeanLoad;
    vector<float> cluster0_4_pcaEigenValLoad;
    vector<float> cluster0_4_pcaEigenVecLoad;
    readFromCommaDelimitedFile_Float(cluster0_4_pcaMeanLoadPATH, cluster0_4_pcaMeanLoad);
    readFromCommaDelimitedFile_Float(cluster0_4_pcaEigenValLoadPATH, cluster0_4_pcaEigenValLoad);
    readFromCommaDelimitedFile_Float(cluster0_4_pcaEigenVecLoadPATH, cluster0_4_pcaEigenVecLoad);
    cv::PCA pcaANN0_4;
    pcaANN0_4.mean = cv::Mat(1, numVarANN, CV_32FC1, cluster0_4_pcaMeanLoad.data());
    pcaANN0_4.eigenvalues = cv::Mat(numComponentPCA,1, CV_32FC1, cluster0_4_pcaEigenValLoad.data());
    pcaANN0_4.eigenvectors = cv::Mat(numComponentPCA,numVarANN, CV_32FC1, cluster0_4_pcaEigenVecLoad.data());
    
    std::string cluster0_5_pcaMeanLoadPATH = casePath + cluster0_5_Folder + "PCA_mean.txt";
    std::string cluster0_5_pcaEigenValLoadPATH = casePath + cluster0_5_Folder + "PCA_eigenValues.txt";
    std::string cluster0_5_pcaEigenVecLoadPATH = casePath + cluster0_5_Folder + "PCA_eigenVectors.txt";
    vector<float> cluster0_5_pcaMeanLoad;
    vector<float> cluster0_5_pcaEigenValLoad;
    vector<float> cluster0_5_pcaEigenVecLoad;
    readFromCommaDelimitedFile_Float(cluster0_5_pcaMeanLoadPATH, cluster0_5_pcaMeanLoad);
    readFromCommaDelimitedFile_Float(cluster0_5_pcaEigenValLoadPATH, cluster0_5_pcaEigenValLoad);
    readFromCommaDelimitedFile_Float(cluster0_5_pcaEigenVecLoadPATH, cluster0_5_pcaEigenVecLoad);
    cv::PCA pcaANN0_5;
    pcaANN0_5.mean = cv::Mat(1, numVarANN, CV_32FC1, cluster0_5_pcaMeanLoad.data());
    pcaANN0_5.eigenvalues = cv::Mat(numComponentPCA,1, CV_32FC1, cluster0_5_pcaEigenValLoad.data());
    pcaANN0_5.eigenvectors = cv::Mat(numComponentPCA,numVarANN, CV_32FC1, cluster0_5_pcaEigenVecLoad.data());
    
    std::string cluster0_6_pcaMeanLoadPATH = casePath + cluster0_6_Folder + "PCA_mean.txt";
    std::string cluster0_6_pcaEigenValLoadPATH = casePath + cluster0_6_Folder + "PCA_eigenValues.txt";
    std::string cluster0_6_pcaEigenVecLoadPATH = casePath + cluster0_6_Folder + "PCA_eigenVectors.txt";
    vector<float> cluster0_6_pcaMeanLoad;
    vector<float> cluster0_6_pcaEigenValLoad;
    vector<float> cluster0_6_pcaEigenVecLoad;
    readFromCommaDelimitedFile_Float(cluster0_6_pcaMeanLoadPATH, cluster0_6_pcaMeanLoad);
    readFromCommaDelimitedFile_Float(cluster0_6_pcaEigenValLoadPATH, cluster0_6_pcaEigenValLoad);
    readFromCommaDelimitedFile_Float(cluster0_6_pcaEigenVecLoadPATH, cluster0_6_pcaEigenVecLoad);
    cv::PCA pcaANN0_6;
    pcaANN0_6.mean = cv::Mat(1, numVarANN, CV_32FC1, cluster0_6_pcaMeanLoad.data());
    pcaANN0_6.eigenvalues = cv::Mat(numComponentPCA,1, CV_32FC1, cluster0_6_pcaEigenValLoad.data());
    pcaANN0_6.eigenvectors = cv::Mat(numComponentPCA,numVarANN, CV_32FC1, cluster0_6_pcaEigenVecLoad.data());
    
    std::string cluster0_7_pcaMeanLoadPATH = casePath + cluster0_7_Folder + "PCA_mean.txt";
    std::string cluster0_7_pcaEigenValLoadPATH = casePath + cluster0_7_Folder + "PCA_eigenValues.txt";
    std::string cluster0_7_pcaEigenVecLoadPATH = casePath + cluster0_7_Folder + "PCA_eigenVectors.txt";
    vector<float> cluster0_7_pcaMeanLoad;
    vector<float> cluster0_7_pcaEigenValLoad;
    vector<float> cluster0_7_pcaEigenVecLoad;
    readFromCommaDelimitedFile_Float(cluster0_7_pcaMeanLoadPATH, cluster0_7_pcaMeanLoad);
    readFromCommaDelimitedFile_Float(cluster0_7_pcaEigenValLoadPATH, cluster0_7_pcaEigenValLoad);
    readFromCommaDelimitedFile_Float(cluster0_7_pcaEigenVecLoadPATH, cluster0_7_pcaEigenVecLoad);
    cv::PCA pcaANN0_7;
    pcaANN0_7.mean = cv::Mat(1, numVarANN, CV_32FC1, cluster0_7_pcaMeanLoad.data());
    pcaANN0_7.eigenvalues = cv::Mat(numComponentPCA,1, CV_32FC1, cluster0_7_pcaEigenValLoad.data());
    pcaANN0_7.eigenvectors = cv::Mat(numComponentPCA,numVarANN, CV_32FC1, cluster0_7_pcaEigenVecLoad.data());
    
    std::string cluster0_8_pcaMeanLoadPATH = casePath + cluster0_8_Folder + "PCA_mean.txt";
    std::string cluster0_8_pcaEigenValLoadPATH = casePath + cluster0_8_Folder + "PCA_eigenValues.txt";
    std::string cluster0_8_pcaEigenVecLoadPATH = casePath + cluster0_8_Folder + "PCA_eigenVectors.txt";
    vector<float> cluster0_8_pcaMeanLoad;
    vector<float> cluster0_8_pcaEigenValLoad;
    vector<float> cluster0_8_pcaEigenVecLoad;
    readFromCommaDelimitedFile_Float(cluster0_8_pcaMeanLoadPATH, cluster0_8_pcaMeanLoad);
    readFromCommaDelimitedFile_Float(cluster0_8_pcaEigenValLoadPATH, cluster0_8_pcaEigenValLoad);
    readFromCommaDelimitedFile_Float(cluster0_8_pcaEigenVecLoadPATH, cluster0_8_pcaEigenVecLoad);
    cv::PCA pcaANN0_8;
    pcaANN0_8.mean = cv::Mat(1, numVarANN, CV_32FC1, cluster0_8_pcaMeanLoad.data());
    pcaANN0_8.eigenvalues = cv::Mat(numComponentPCA,1, CV_32FC1, cluster0_8_pcaEigenValLoad.data());
    pcaANN0_8.eigenvectors = cv::Mat(numComponentPCA,numVarANN, CV_32FC1, cluster0_8_pcaEigenVecLoad.data());
    
    std::string cluster0_9_pcaMeanLoadPATH = casePath + cluster0_9_Folder + "PCA_mean.txt";
    std::string cluster0_9_pcaEigenValLoadPATH = casePath + cluster0_9_Folder + "PCA_eigenValues.txt";
    std::string cluster0_9_pcaEigenVecLoadPATH = casePath + cluster0_9_Folder + "PCA_eigenVectors.txt";
    vector<float> cluster0_9_pcaMeanLoad;
    vector<float> cluster0_9_pcaEigenValLoad;
    vector<float> cluster0_9_pcaEigenVecLoad;
    readFromCommaDelimitedFile_Float(cluster0_9_pcaMeanLoadPATH, cluster0_9_pcaMeanLoad);
    readFromCommaDelimitedFile_Float(cluster0_9_pcaEigenValLoadPATH, cluster0_9_pcaEigenValLoad);
    readFromCommaDelimitedFile_Float(cluster0_9_pcaEigenVecLoadPATH, cluster0_9_pcaEigenVecLoad);
    cv::PCA pcaANN0_9;
    pcaANN0_9.mean = cv::Mat(1, numVarANN, CV_32FC1, cluster0_9_pcaMeanLoad.data());
    pcaANN0_9.eigenvalues = cv::Mat(numComponentPCA,1, CV_32FC1, cluster0_9_pcaEigenValLoad.data());
    pcaANN0_9.eigenvectors = cv::Mat(numComponentPCA,numVarANN, CV_32FC1, cluster0_9_pcaEigenVecLoad.data());
    
    std::string cluster0_10_pcaMeanLoadPATH = casePath + cluster0_10_Folder + "PCA_mean.txt";
    std::string cluster0_10_pcaEigenValLoadPATH = casePath + cluster0_10_Folder + "PCA_eigenValues.txt";
    std::string cluster0_10_pcaEigenVecLoadPATH = casePath + cluster0_10_Folder + "PCA_eigenVectors.txt";
    vector<float> cluster0_10_pcaMeanLoad;
    vector<float> cluster0_10_pcaEigenValLoad;
    vector<float> cluster0_10_pcaEigenVecLoad;
    readFromCommaDelimitedFile_Float(cluster0_10_pcaMeanLoadPATH, cluster0_10_pcaMeanLoad);
    readFromCommaDelimitedFile_Float(cluster0_10_pcaEigenValLoadPATH, cluster0_10_pcaEigenValLoad);
    readFromCommaDelimitedFile_Float(cluster0_10_pcaEigenVecLoadPATH, cluster0_10_pcaEigenVecLoad);
    cv::PCA pcaANN0_10;
    pcaANN0_10.mean = cv::Mat(1, numVarANN, CV_32FC1, cluster0_10_pcaMeanLoad.data());
    pcaANN0_10.eigenvalues = cv::Mat(numComponentPCA,1, CV_32FC1, cluster0_10_pcaEigenValLoad.data());
    pcaANN0_10.eigenvectors = cv::Mat(numComponentPCA,numVarANN, CV_32FC1, cluster0_10_pcaEigenVecLoad.data());
    
    std::string cluster0_11_pcaMeanLoadPATH = casePath + cluster0_11_Folder + "PCA_mean.txt";
    std::string cluster0_11_pcaEigenValLoadPATH = casePath + cluster0_11_Folder + "PCA_eigenValues.txt";
    std::string cluster0_11_pcaEigenVecLoadPATH = casePath + cluster0_11_Folder + "PCA_eigenVectors.txt";
    vector<float> cluster0_11_pcaMeanLoad;
    vector<float> cluster0_11_pcaEigenValLoad;
    vector<float> cluster0_11_pcaEigenVecLoad;
    readFromCommaDelimitedFile_Float(cluster0_11_pcaMeanLoadPATH, cluster0_11_pcaMeanLoad);
    readFromCommaDelimitedFile_Float(cluster0_11_pcaEigenValLoadPATH, cluster0_11_pcaEigenValLoad);
    readFromCommaDelimitedFile_Float(cluster0_11_pcaEigenVecLoadPATH, cluster0_11_pcaEigenVecLoad);
    cv::PCA pcaANN0_11;
    pcaANN0_11.mean = cv::Mat(1, numVarANN, CV_32FC1, cluster0_11_pcaMeanLoad.data());
    pcaANN0_11.eigenvalues = cv::Mat(numComponentPCA,1, CV_32FC1, cluster0_11_pcaEigenValLoad.data());
    pcaANN0_11.eigenvectors = cv::Mat(numComponentPCA,numVarANN, CV_32FC1, cluster0_11_pcaEigenVecLoad.data());
    
    std::string cluster0_12_pcaMeanLoadPATH = casePath + cluster0_12_Folder + "PCA_mean.txt";
    std::string cluster0_12_pcaEigenValLoadPATH = casePath + cluster0_12_Folder + "PCA_eigenValues.txt";
    std::string cluster0_12_pcaEigenVecLoadPATH = casePath + cluster0_12_Folder + "PCA_eigenVectors.txt";
    vector<float> cluster0_12_pcaMeanLoad;
    vector<float> cluster0_12_pcaEigenValLoad;
    vector<float> cluster0_12_pcaEigenVecLoad;
    readFromCommaDelimitedFile_Float(cluster0_12_pcaMeanLoadPATH, cluster0_12_pcaMeanLoad);
    readFromCommaDelimitedFile_Float(cluster0_12_pcaEigenValLoadPATH, cluster0_12_pcaEigenValLoad);
    readFromCommaDelimitedFile_Float(cluster0_12_pcaEigenVecLoadPATH, cluster0_12_pcaEigenVecLoad);
    cv::PCA pcaANN0_12;
    pcaANN0_12.mean = cv::Mat(1, numVarANN, CV_32FC1, cluster0_12_pcaMeanLoad.data());
    pcaANN0_12.eigenvalues = cv::Mat(numComponentPCA,1, CV_32FC1, cluster0_12_pcaEigenValLoad.data());
    pcaANN0_12.eigenvectors = cv::Mat(numComponentPCA,numVarANN, CV_32FC1, cluster0_12_pcaEigenVecLoad.data());
    
    std::string cluster0_13_pcaMeanLoadPATH = casePath + cluster0_13_Folder + "PCA_mean.txt";
    std::string cluster0_13_pcaEigenValLoadPATH = casePath + cluster0_13_Folder + "PCA_eigenValues.txt";
    std::string cluster0_13_pcaEigenVecLoadPATH = casePath + cluster0_13_Folder + "PCA_eigenVectors.txt";
    vector<float> cluster0_13_pcaMeanLoad;
    vector<float> cluster0_13_pcaEigenValLoad;
    vector<float> cluster0_13_pcaEigenVecLoad;
    readFromCommaDelimitedFile_Float(cluster0_13_pcaMeanLoadPATH, cluster0_13_pcaMeanLoad);
    readFromCommaDelimitedFile_Float(cluster0_13_pcaEigenValLoadPATH, cluster0_13_pcaEigenValLoad);
    readFromCommaDelimitedFile_Float(cluster0_13_pcaEigenVecLoadPATH, cluster0_13_pcaEigenVecLoad);
    cv::PCA pcaANN0_13;
    pcaANN0_13.mean = cv::Mat(1, numVarANN, CV_32FC1, cluster0_13_pcaMeanLoad.data());
    pcaANN0_13.eigenvalues = cv::Mat(numComponentPCA,1, CV_32FC1, cluster0_13_pcaEigenValLoad.data());
    pcaANN0_13.eigenvectors = cv::Mat(numComponentPCA,numVarANN, CV_32FC1, cluster0_13_pcaEigenVecLoad.data());
    
    std::string cluster0_14_pcaMeanLoadPATH = casePath + cluster0_14_Folder + "PCA_mean.txt";
    std::string cluster0_14_pcaEigenValLoadPATH = casePath + cluster0_14_Folder + "PCA_eigenValues.txt";
    std::string cluster0_14_pcaEigenVecLoadPATH = casePath + cluster0_14_Folder + "PCA_eigenVectors.txt";
    vector<float> cluster0_14_pcaMeanLoad;
    vector<float> cluster0_14_pcaEigenValLoad;
    vector<float> cluster0_14_pcaEigenVecLoad;
    readFromCommaDelimitedFile_Float(cluster0_14_pcaMeanLoadPATH, cluster0_14_pcaMeanLoad);
    readFromCommaDelimitedFile_Float(cluster0_14_pcaEigenValLoadPATH, cluster0_14_pcaEigenValLoad);
    readFromCommaDelimitedFile_Float(cluster0_14_pcaEigenVecLoadPATH, cluster0_14_pcaEigenVecLoad);
    cv::PCA pcaANN0_14;
    pcaANN0_14.mean = cv::Mat(1, numVarANN, CV_32FC1, cluster0_14_pcaMeanLoad.data());
    pcaANN0_14.eigenvalues = cv::Mat(numComponentPCA,1, CV_32FC1, cluster0_14_pcaEigenValLoad.data());
    pcaANN0_14.eigenvectors = cv::Mat(numComponentPCA,numVarANN, CV_32FC1, cluster0_14_pcaEigenVecLoad.data());
    
    std::string cluster0_15_pcaMeanLoadPATH = casePath + cluster0_15_Folder + "PCA_mean.txt";
    std::string cluster0_15_pcaEigenValLoadPATH = casePath + cluster0_15_Folder + "PCA_eigenValues.txt";
    std::string cluster0_15_pcaEigenVecLoadPATH = casePath + cluster0_15_Folder + "PCA_eigenVectors.txt";
    vector<float> cluster0_15_pcaMeanLoad;
    vector<float> cluster0_15_pcaEigenValLoad;
    vector<float> cluster0_15_pcaEigenVecLoad;
    readFromCommaDelimitedFile_Float(cluster0_15_pcaMeanLoadPATH, cluster0_15_pcaMeanLoad);
    readFromCommaDelimitedFile_Float(cluster0_15_pcaEigenValLoadPATH, cluster0_15_pcaEigenValLoad);
    readFromCommaDelimitedFile_Float(cluster0_15_pcaEigenVecLoadPATH, cluster0_15_pcaEigenVecLoad);
    cv::PCA pcaANN0_15;
    pcaANN0_15.mean = cv::Mat(1, numVarANN, CV_32FC1, cluster0_15_pcaMeanLoad.data());
    pcaANN0_15.eigenvalues = cv::Mat(numComponentPCA,1, CV_32FC1, cluster0_15_pcaEigenValLoad.data());
    pcaANN0_15.eigenvectors = cv::Mat(numComponentPCA,numVarANN, CV_32FC1, cluster0_15_pcaEigenVecLoad.data());
    
    std::string cluster0_16_pcaMeanLoadPATH = casePath + cluster0_16_Folder + "PCA_mean.txt";
    std::string cluster0_16_pcaEigenValLoadPATH = casePath + cluster0_16_Folder + "PCA_eigenValues.txt";
    std::string cluster0_16_pcaEigenVecLoadPATH = casePath + cluster0_16_Folder + "PCA_eigenVectors.txt";
    vector<float> cluster0_16_pcaMeanLoad;
    vector<float> cluster0_16_pcaEigenValLoad;
    vector<float> cluster0_16_pcaEigenVecLoad;
    readFromCommaDelimitedFile_Float(cluster0_16_pcaMeanLoadPATH, cluster0_16_pcaMeanLoad);
    readFromCommaDelimitedFile_Float(cluster0_16_pcaEigenValLoadPATH, cluster0_16_pcaEigenValLoad);
    readFromCommaDelimitedFile_Float(cluster0_16_pcaEigenVecLoadPATH, cluster0_16_pcaEigenVecLoad);
    cv::PCA pcaANN0_16;
    pcaANN0_16.mean = cv::Mat(1, numVarANN, CV_32FC1, cluster0_16_pcaMeanLoad.data());
    pcaANN0_16.eigenvalues = cv::Mat(numComponentPCA,1, CV_32FC1, cluster0_16_pcaEigenValLoad.data());
    pcaANN0_16.eigenvectors = cv::Mat(numComponentPCA,numVarANN, CV_32FC1, cluster0_16_pcaEigenVecLoad.data());
    
    std::string cluster0_17_pcaMeanLoadPATH = casePath + cluster0_17_Folder + "PCA_mean.txt";
    std::string cluster0_17_pcaEigenValLoadPATH = casePath + cluster0_17_Folder + "PCA_eigenValues.txt";
    std::string cluster0_17_pcaEigenVecLoadPATH = casePath + cluster0_17_Folder + "PCA_eigenVectors.txt";
    vector<float> cluster0_17_pcaMeanLoad;
    vector<float> cluster0_17_pcaEigenValLoad;
    vector<float> cluster0_17_pcaEigenVecLoad;
    readFromCommaDelimitedFile_Float(cluster0_17_pcaMeanLoadPATH, cluster0_17_pcaMeanLoad);
    readFromCommaDelimitedFile_Float(cluster0_17_pcaEigenValLoadPATH, cluster0_17_pcaEigenValLoad);
    readFromCommaDelimitedFile_Float(cluster0_17_pcaEigenVecLoadPATH, cluster0_17_pcaEigenVecLoad);
    cv::PCA pcaANN0_17;
    pcaANN0_17.mean = cv::Mat(1, numVarANN, CV_32FC1, cluster0_17_pcaMeanLoad.data());
    pcaANN0_17.eigenvalues = cv::Mat(numComponentPCA,1, CV_32FC1, cluster0_17_pcaEigenValLoad.data());
    pcaANN0_17.eigenvectors = cv::Mat(numComponentPCA,numVarANN, CV_32FC1, cluster0_17_pcaEigenVecLoad.data());
    
    std::string cluster0_18_pcaMeanLoadPATH = casePath + cluster0_18_Folder + "PCA_mean.txt";
    std::string cluster0_18_pcaEigenValLoadPATH = casePath + cluster0_18_Folder + "PCA_eigenValues.txt";
    std::string cluster0_18_pcaEigenVecLoadPATH = casePath + cluster0_18_Folder + "PCA_eigenVectors.txt";
    vector<float> cluster0_18_pcaMeanLoad;
    vector<float> cluster0_18_pcaEigenValLoad;
    vector<float> cluster0_18_pcaEigenVecLoad;
    readFromCommaDelimitedFile_Float(cluster0_18_pcaMeanLoadPATH, cluster0_18_pcaMeanLoad);
    readFromCommaDelimitedFile_Float(cluster0_18_pcaEigenValLoadPATH, cluster0_18_pcaEigenValLoad);
    readFromCommaDelimitedFile_Float(cluster0_18_pcaEigenVecLoadPATH, cluster0_18_pcaEigenVecLoad);
    cv::PCA pcaANN0_18;
    pcaANN0_18.mean = cv::Mat(1, numVarANN, CV_32FC1, cluster0_18_pcaMeanLoad.data());
    pcaANN0_18.eigenvalues = cv::Mat(numComponentPCA,1, CV_32FC1, cluster0_18_pcaEigenValLoad.data());
    pcaANN0_18.eigenvectors = cv::Mat(numComponentPCA,numVarANN, CV_32FC1, cluster0_18_pcaEigenVecLoad.data());
    
    std::string cluster0_19_pcaMeanLoadPATH = casePath + cluster0_19_Folder + "PCA_mean.txt";
    std::string cluster0_19_pcaEigenValLoadPATH = casePath + cluster0_19_Folder + "PCA_eigenValues.txt";
    std::string cluster0_19_pcaEigenVecLoadPATH = casePath + cluster0_19_Folder + "PCA_eigenVectors.txt";
    vector<float> cluster0_19_pcaMeanLoad;
    vector<float> cluster0_19_pcaEigenValLoad;
    vector<float> cluster0_19_pcaEigenVecLoad;
    readFromCommaDelimitedFile_Float(cluster0_19_pcaMeanLoadPATH, cluster0_19_pcaMeanLoad);
    readFromCommaDelimitedFile_Float(cluster0_19_pcaEigenValLoadPATH, cluster0_19_pcaEigenValLoad);
    readFromCommaDelimitedFile_Float(cluster0_19_pcaEigenVecLoadPATH, cluster0_19_pcaEigenVecLoad);
    cv::PCA pcaANN0_19;
    pcaANN0_19.mean = cv::Mat(1, numVarANN, CV_32FC1, cluster0_19_pcaMeanLoad.data());
    pcaANN0_19.eigenvalues = cv::Mat(numComponentPCA,1, CV_32FC1, cluster0_19_pcaEigenValLoad.data());
    pcaANN0_19.eigenvectors = cv::Mat(numComponentPCA,numVarANN, CV_32FC1, cluster0_19_pcaEigenVecLoad.data());
    
    // Child of parentCluster1 : Not necessary
    // GrandChild of parentCluster1 - childCluster0 - grandChild
    std::string cluster1_0_0_pcaMeanLoadPATH = casePath + cluster1_0_0_Folder + "PCA_mean.txt";
    std::string cluster1_0_0_pcaEigenValLoadPATH = casePath + cluster1_0_0_Folder + "PCA_eigenValues.txt";
    std::string cluster1_0_0_pcaEigenVecLoadPATH = casePath + cluster1_0_0_Folder + "PCA_eigenVectors.txt";
    vector<float> cluster1_0_0_pcaMeanLoad;
    vector<float> cluster1_0_0_pcaEigenValLoad;
    vector<float> cluster1_0_0_pcaEigenVecLoad;
    readFromCommaDelimitedFile_Float(cluster1_0_0_pcaMeanLoadPATH, cluster1_0_0_pcaMeanLoad);
    readFromCommaDelimitedFile_Float(cluster1_0_0_pcaEigenValLoadPATH, cluster1_0_0_pcaEigenValLoad);
    readFromCommaDelimitedFile_Float(cluster1_0_0_pcaEigenVecLoadPATH, cluster1_0_0_pcaEigenVecLoad);
    cv::PCA pcaANN1_0_0;
    pcaANN1_0_0.mean = cv::Mat(1, numVarANN, CV_32FC1, cluster1_0_0_pcaMeanLoad.data());
    pcaANN1_0_0.eigenvalues = cv::Mat(numComponentPCA,1, CV_32FC1, cluster1_0_0_pcaEigenValLoad.data());
    pcaANN1_0_0.eigenvectors = cv::Mat(numComponentPCA,numVarANN, CV_32FC1, cluster1_0_0_pcaEigenVecLoad.data());
    
    std::string cluster1_0_1_pcaMeanLoadPATH = casePath + cluster1_0_1_Folder + "PCA_mean.txt";
    std::string cluster1_0_1_pcaEigenValLoadPATH = casePath + cluster1_0_1_Folder + "PCA_eigenValues.txt";
    std::string cluster1_0_1_pcaEigenVecLoadPATH = casePath + cluster1_0_1_Folder + "PCA_eigenVectors.txt";
    vector<float> cluster1_0_1_pcaMeanLoad;
    vector<float> cluster1_0_1_pcaEigenValLoad;
    vector<float> cluster1_0_1_pcaEigenVecLoad;
    readFromCommaDelimitedFile_Float(cluster1_0_1_pcaMeanLoadPATH, cluster1_0_1_pcaMeanLoad);
    readFromCommaDelimitedFile_Float(cluster1_0_1_pcaEigenValLoadPATH, cluster1_0_1_pcaEigenValLoad);
    readFromCommaDelimitedFile_Float(cluster1_0_1_pcaEigenVecLoadPATH, cluster1_0_1_pcaEigenVecLoad);
    cv::PCA pcaANN1_0_1;
    pcaANN1_0_1.mean = cv::Mat(1, numVarANN, CV_32FC1, cluster1_0_1_pcaMeanLoad.data());
    pcaANN1_0_1.eigenvalues = cv::Mat(numComponentPCA,1, CV_32FC1, cluster1_0_1_pcaEigenValLoad.data());
    pcaANN1_0_1.eigenvectors = cv::Mat(numComponentPCA,numVarANN, CV_32FC1, cluster1_0_1_pcaEigenVecLoad.data());
    
    std::string cluster1_0_2_pcaMeanLoadPATH = casePath + cluster1_0_2_Folder + "PCA_mean.txt";
    std::string cluster1_0_2_pcaEigenValLoadPATH = casePath + cluster1_0_2_Folder + "PCA_eigenValues.txt";
    std::string cluster1_0_2_pcaEigenVecLoadPATH = casePath + cluster1_0_2_Folder + "PCA_eigenVectors.txt";
    vector<float> cluster1_0_2_pcaMeanLoad;
    vector<float> cluster1_0_2_pcaEigenValLoad;
    vector<float> cluster1_0_2_pcaEigenVecLoad;
    readFromCommaDelimitedFile_Float(cluster1_0_2_pcaMeanLoadPATH, cluster1_0_2_pcaMeanLoad);
    readFromCommaDelimitedFile_Float(cluster1_0_2_pcaEigenValLoadPATH, cluster1_0_2_pcaEigenValLoad);
    readFromCommaDelimitedFile_Float(cluster1_0_2_pcaEigenVecLoadPATH, cluster1_0_2_pcaEigenVecLoad);
    cv::PCA pcaANN1_0_2;
    pcaANN1_0_2.mean = cv::Mat(1, numVarANN, CV_32FC1, cluster1_0_2_pcaMeanLoad.data());
    pcaANN1_0_2.eigenvalues = cv::Mat(numComponentPCA,1, CV_32FC1, cluster1_0_2_pcaEigenValLoad.data());
    pcaANN1_0_2.eigenvectors = cv::Mat(numComponentPCA,numVarANN, CV_32FC1, cluster1_0_2_pcaEigenVecLoad.data());
    
    std::string cluster1_0_3_pcaMeanLoadPATH = casePath + cluster1_0_3_Folder + "PCA_mean.txt";
    std::string cluster1_0_3_pcaEigenValLoadPATH = casePath + cluster1_0_3_Folder + "PCA_eigenValues.txt";
    std::string cluster1_0_3_pcaEigenVecLoadPATH = casePath + cluster1_0_3_Folder + "PCA_eigenVectors.txt";
    vector<float> cluster1_0_3_pcaMeanLoad;
    vector<float> cluster1_0_3_pcaEigenValLoad;
    vector<float> cluster1_0_3_pcaEigenVecLoad;
    readFromCommaDelimitedFile_Float(cluster1_0_3_pcaMeanLoadPATH, cluster1_0_3_pcaMeanLoad);
    readFromCommaDelimitedFile_Float(cluster1_0_3_pcaEigenValLoadPATH, cluster1_0_3_pcaEigenValLoad);
    readFromCommaDelimitedFile_Float(cluster1_0_3_pcaEigenVecLoadPATH, cluster1_0_3_pcaEigenVecLoad);
    cv::PCA pcaANN1_0_3;
    pcaANN1_0_3.mean = cv::Mat(1, numVarANN, CV_32FC1, cluster1_0_3_pcaMeanLoad.data());
    pcaANN1_0_3.eigenvalues = cv::Mat(numComponentPCA,1, CV_32FC1, cluster1_0_3_pcaEigenValLoad.data());
    pcaANN1_0_3.eigenvectors = cv::Mat(numComponentPCA,numVarANN, CV_32FC1, cluster1_0_3_pcaEigenVecLoad.data());
    
    std::string cluster1_0_4_pcaMeanLoadPATH = casePath + cluster1_0_4_Folder + "PCA_mean.txt";
    std::string cluster1_0_4_pcaEigenValLoadPATH = casePath + cluster1_0_4_Folder + "PCA_eigenValues.txt";
    std::string cluster1_0_4_pcaEigenVecLoadPATH = casePath + cluster1_0_4_Folder + "PCA_eigenVectors.txt";
    vector<float> cluster1_0_4_pcaMeanLoad;
    vector<float> cluster1_0_4_pcaEigenValLoad;
    vector<float> cluster1_0_4_pcaEigenVecLoad;
    readFromCommaDelimitedFile_Float(cluster1_0_4_pcaMeanLoadPATH, cluster1_0_4_pcaMeanLoad);
    readFromCommaDelimitedFile_Float(cluster1_0_4_pcaEigenValLoadPATH, cluster1_0_4_pcaEigenValLoad);
    readFromCommaDelimitedFile_Float(cluster1_0_4_pcaEigenVecLoadPATH, cluster1_0_4_pcaEigenVecLoad);
    cv::PCA pcaANN1_0_4;
    pcaANN1_0_4.mean = cv::Mat(1, numVarANN, CV_32FC1, cluster1_0_4_pcaMeanLoad.data());
    pcaANN1_0_4.eigenvalues = cv::Mat(numComponentPCA,1, CV_32FC1, cluster1_0_4_pcaEigenValLoad.data());
    pcaANN1_0_4.eigenvectors = cv::Mat(numComponentPCA,numVarANN, CV_32FC1, cluster1_0_4_pcaEigenVecLoad.data());
    
    std::string cluster1_0_5_pcaMeanLoadPATH = casePath + cluster1_0_5_Folder + "PCA_mean.txt";
    std::string cluster1_0_5_pcaEigenValLoadPATH = casePath + cluster1_0_5_Folder + "PCA_eigenValues.txt";
    std::string cluster1_0_5_pcaEigenVecLoadPATH = casePath + cluster1_0_5_Folder + "PCA_eigenVectors.txt";
    vector<float> cluster1_0_5_pcaMeanLoad;
    vector<float> cluster1_0_5_pcaEigenValLoad;
    vector<float> cluster1_0_5_pcaEigenVecLoad;
    readFromCommaDelimitedFile_Float(cluster1_0_5_pcaMeanLoadPATH, cluster1_0_5_pcaMeanLoad);
    readFromCommaDelimitedFile_Float(cluster1_0_5_pcaEigenValLoadPATH, cluster1_0_5_pcaEigenValLoad);
    readFromCommaDelimitedFile_Float(cluster1_0_5_pcaEigenVecLoadPATH, cluster1_0_5_pcaEigenVecLoad);
    cv::PCA pcaANN1_0_5;
    pcaANN1_0_5.mean = cv::Mat(1, numVarANN, CV_32FC1, cluster1_0_5_pcaMeanLoad.data());
    pcaANN1_0_5.eigenvalues = cv::Mat(numComponentPCA,1, CV_32FC1, cluster1_0_5_pcaEigenValLoad.data());
    pcaANN1_0_5.eigenvectors = cv::Mat(numComponentPCA,numVarANN, CV_32FC1, cluster1_0_5_pcaEigenVecLoad.data());
    
    std::string cluster1_0_6_pcaMeanLoadPATH = casePath + cluster1_0_6_Folder + "PCA_mean.txt";
    std::string cluster1_0_6_pcaEigenValLoadPATH = casePath + cluster1_0_6_Folder + "PCA_eigenValues.txt";
    std::string cluster1_0_6_pcaEigenVecLoadPATH = casePath + cluster1_0_6_Folder + "PCA_eigenVectors.txt";
    vector<float> cluster1_0_6_pcaMeanLoad;
    vector<float> cluster1_0_6_pcaEigenValLoad;
    vector<float> cluster1_0_6_pcaEigenVecLoad;
    readFromCommaDelimitedFile_Float(cluster1_0_6_pcaMeanLoadPATH, cluster1_0_6_pcaMeanLoad);
    readFromCommaDelimitedFile_Float(cluster1_0_6_pcaEigenValLoadPATH, cluster1_0_6_pcaEigenValLoad);
    readFromCommaDelimitedFile_Float(cluster1_0_6_pcaEigenVecLoadPATH, cluster1_0_6_pcaEigenVecLoad);
    cv::PCA pcaANN1_0_6;
    pcaANN1_0_6.mean = cv::Mat(1, numVarANN, CV_32FC1, cluster1_0_6_pcaMeanLoad.data());
    pcaANN1_0_6.eigenvalues = cv::Mat(numComponentPCA,1, CV_32FC1, cluster1_0_6_pcaEigenValLoad.data());
    pcaANN1_0_6.eigenvectors = cv::Mat(numComponentPCA,numVarANN, CV_32FC1, cluster1_0_6_pcaEigenVecLoad.data());
    
    std::string cluster1_0_7_pcaMeanLoadPATH = casePath + cluster1_0_7_Folder + "PCA_mean.txt";
    std::string cluster1_0_7_pcaEigenValLoadPATH = casePath + cluster1_0_7_Folder + "PCA_eigenValues.txt";
    std::string cluster1_0_7_pcaEigenVecLoadPATH = casePath + cluster1_0_7_Folder + "PCA_eigenVectors.txt";
    vector<float> cluster1_0_7_pcaMeanLoad;
    vector<float> cluster1_0_7_pcaEigenValLoad;
    vector<float> cluster1_0_7_pcaEigenVecLoad;
    readFromCommaDelimitedFile_Float(cluster1_0_7_pcaMeanLoadPATH, cluster1_0_7_pcaMeanLoad);
    readFromCommaDelimitedFile_Float(cluster1_0_7_pcaEigenValLoadPATH, cluster1_0_7_pcaEigenValLoad);
    readFromCommaDelimitedFile_Float(cluster1_0_7_pcaEigenVecLoadPATH, cluster1_0_7_pcaEigenVecLoad);
    cv::PCA pcaANN1_0_7;
    pcaANN1_0_7.mean = cv::Mat(1, numVarANN, CV_32FC1, cluster1_0_7_pcaMeanLoad.data());
    pcaANN1_0_7.eigenvalues = cv::Mat(numComponentPCA,1, CV_32FC1, cluster1_0_7_pcaEigenValLoad.data());
    pcaANN1_0_7.eigenvectors = cv::Mat(numComponentPCA,numVarANN, CV_32FC1, cluster1_0_7_pcaEigenVecLoad.data());
    
    std::string cluster1_0_8_pcaMeanLoadPATH = casePath + cluster1_0_8_Folder + "PCA_mean.txt";
    std::string cluster1_0_8_pcaEigenValLoadPATH = casePath + cluster1_0_8_Folder + "PCA_eigenValues.txt";
    std::string cluster1_0_8_pcaEigenVecLoadPATH = casePath + cluster1_0_8_Folder + "PCA_eigenVectors.txt";
    vector<float> cluster1_0_8_pcaMeanLoad;
    vector<float> cluster1_0_8_pcaEigenValLoad;
    vector<float> cluster1_0_8_pcaEigenVecLoad;
    readFromCommaDelimitedFile_Float(cluster1_0_8_pcaMeanLoadPATH, cluster1_0_8_pcaMeanLoad);
    readFromCommaDelimitedFile_Float(cluster1_0_8_pcaEigenValLoadPATH, cluster1_0_8_pcaEigenValLoad);
    readFromCommaDelimitedFile_Float(cluster1_0_8_pcaEigenVecLoadPATH, cluster1_0_8_pcaEigenVecLoad);
    cv::PCA pcaANN1_0_8;
    pcaANN1_0_8.mean = cv::Mat(1, numVarANN, CV_32FC1, cluster1_0_8_pcaMeanLoad.data());
    pcaANN1_0_8.eigenvalues = cv::Mat(numComponentPCA,1, CV_32FC1, cluster1_0_8_pcaEigenValLoad.data());
    pcaANN1_0_8.eigenvectors = cv::Mat(numComponentPCA,numVarANN, CV_32FC1, cluster1_0_8_pcaEigenVecLoad.data());
    
    std::string cluster1_0_9_pcaMeanLoadPATH = casePath + cluster1_0_9_Folder + "PCA_mean.txt";
    std::string cluster1_0_9_pcaEigenValLoadPATH = casePath + cluster1_0_9_Folder + "PCA_eigenValues.txt";
    std::string cluster1_0_9_pcaEigenVecLoadPATH = casePath + cluster1_0_9_Folder + "PCA_eigenVectors.txt";
    vector<float> cluster1_0_9_pcaMeanLoad;
    vector<float> cluster1_0_9_pcaEigenValLoad;
    vector<float> cluster1_0_9_pcaEigenVecLoad;
    readFromCommaDelimitedFile_Float(cluster1_0_9_pcaMeanLoadPATH, cluster1_0_9_pcaMeanLoad);
    readFromCommaDelimitedFile_Float(cluster1_0_9_pcaEigenValLoadPATH, cluster1_0_9_pcaEigenValLoad);
    readFromCommaDelimitedFile_Float(cluster1_0_9_pcaEigenVecLoadPATH, cluster1_0_9_pcaEigenVecLoad);
    cv::PCA pcaANN1_0_9;
    pcaANN1_0_9.mean = cv::Mat(1, numVarANN, CV_32FC1, cluster1_0_9_pcaMeanLoad.data());
    pcaANN1_0_9.eigenvalues = cv::Mat(numComponentPCA,1, CV_32FC1, cluster1_0_9_pcaEigenValLoad.data());
    pcaANN1_0_9.eigenvectors = cv::Mat(numComponentPCA,numVarANN, CV_32FC1, cluster1_0_9_pcaEigenVecLoad.data());
    
    std::string cluster1_0_10_pcaMeanLoadPATH = casePath + cluster1_0_10_Folder + "PCA_mean.txt";
    std::string cluster1_0_10_pcaEigenValLoadPATH = casePath + cluster1_0_10_Folder + "PCA_eigenValues.txt";
    std::string cluster1_0_10_pcaEigenVecLoadPATH = casePath + cluster1_0_10_Folder + "PCA_eigenVectors.txt";
    vector<float> cluster1_0_10_pcaMeanLoad;
    vector<float> cluster1_0_10_pcaEigenValLoad;
    vector<float> cluster1_0_10_pcaEigenVecLoad;
    readFromCommaDelimitedFile_Float(cluster1_0_10_pcaMeanLoadPATH, cluster1_0_10_pcaMeanLoad);
    readFromCommaDelimitedFile_Float(cluster1_0_10_pcaEigenValLoadPATH, cluster1_0_10_pcaEigenValLoad);
    readFromCommaDelimitedFile_Float(cluster1_0_10_pcaEigenVecLoadPATH, cluster1_0_10_pcaEigenVecLoad);
    cv::PCA pcaANN1_0_10;
    pcaANN1_0_10.mean = cv::Mat(1, numVarANN, CV_32FC1, cluster1_0_10_pcaMeanLoad.data());
    pcaANN1_0_10.eigenvalues = cv::Mat(numComponentPCA,1, CV_32FC1, cluster1_0_10_pcaEigenValLoad.data());
    pcaANN1_0_10.eigenvectors = cv::Mat(numComponentPCA,numVarANN, CV_32FC1, cluster1_0_10_pcaEigenVecLoad.data());
    
    std::string cluster1_0_11_pcaMeanLoadPATH = casePath + cluster1_0_11_Folder + "PCA_mean.txt";
    std::string cluster1_0_11_pcaEigenValLoadPATH = casePath + cluster1_0_11_Folder + "PCA_eigenValues.txt";
    std::string cluster1_0_11_pcaEigenVecLoadPATH = casePath + cluster1_0_11_Folder + "PCA_eigenVectors.txt";
    vector<float> cluster1_0_11_pcaMeanLoad;
    vector<float> cluster1_0_11_pcaEigenValLoad;
    vector<float> cluster1_0_11_pcaEigenVecLoad;
    readFromCommaDelimitedFile_Float(cluster1_0_11_pcaMeanLoadPATH, cluster1_0_11_pcaMeanLoad);
    readFromCommaDelimitedFile_Float(cluster1_0_11_pcaEigenValLoadPATH, cluster1_0_11_pcaEigenValLoad);
    readFromCommaDelimitedFile_Float(cluster1_0_11_pcaEigenVecLoadPATH, cluster1_0_11_pcaEigenVecLoad);
    cv::PCA pcaANN1_0_11;
    pcaANN1_0_11.mean = cv::Mat(1, numVarANN, CV_32FC1, cluster1_0_11_pcaMeanLoad.data());
    pcaANN1_0_11.eigenvalues = cv::Mat(numComponentPCA,1, CV_32FC1, cluster1_0_11_pcaEigenValLoad.data());
    pcaANN1_0_11.eigenvectors = cv::Mat(numComponentPCA,numVarANN, CV_32FC1, cluster1_0_11_pcaEigenVecLoad.data());
    
    std::string cluster1_0_12_pcaMeanLoadPATH = casePath + cluster1_0_12_Folder + "PCA_mean.txt";
    std::string cluster1_0_12_pcaEigenValLoadPATH = casePath + cluster1_0_12_Folder + "PCA_eigenValues.txt";
    std::string cluster1_0_12_pcaEigenVecLoadPATH = casePath + cluster1_0_12_Folder + "PCA_eigenVectors.txt";
    vector<float> cluster1_0_12_pcaMeanLoad;
    vector<float> cluster1_0_12_pcaEigenValLoad;
    vector<float> cluster1_0_12_pcaEigenVecLoad;
    readFromCommaDelimitedFile_Float(cluster1_0_12_pcaMeanLoadPATH, cluster1_0_12_pcaMeanLoad);
    readFromCommaDelimitedFile_Float(cluster1_0_12_pcaEigenValLoadPATH, cluster1_0_12_pcaEigenValLoad);
    readFromCommaDelimitedFile_Float(cluster1_0_12_pcaEigenVecLoadPATH, cluster1_0_12_pcaEigenVecLoad);
    cv::PCA pcaANN1_0_12;
    pcaANN1_0_12.mean = cv::Mat(1, numVarANN, CV_32FC1, cluster1_0_12_pcaMeanLoad.data());
    pcaANN1_0_12.eigenvalues = cv::Mat(numComponentPCA,1, CV_32FC1, cluster1_0_12_pcaEigenValLoad.data());
    pcaANN1_0_12.eigenvectors = cv::Mat(numComponentPCA,numVarANN, CV_32FC1, cluster1_0_12_pcaEigenVecLoad.data());
    
    std::string cluster1_0_13_pcaMeanLoadPATH = casePath + cluster1_0_13_Folder + "PCA_mean.txt";
    std::string cluster1_0_13_pcaEigenValLoadPATH = casePath + cluster1_0_13_Folder + "PCA_eigenValues.txt";
    std::string cluster1_0_13_pcaEigenVecLoadPATH = casePath + cluster1_0_13_Folder + "PCA_eigenVectors.txt";
    vector<float> cluster1_0_13_pcaMeanLoad;
    vector<float> cluster1_0_13_pcaEigenValLoad;
    vector<float> cluster1_0_13_pcaEigenVecLoad;
    readFromCommaDelimitedFile_Float(cluster1_0_13_pcaMeanLoadPATH, cluster1_0_13_pcaMeanLoad);
    readFromCommaDelimitedFile_Float(cluster1_0_13_pcaEigenValLoadPATH, cluster1_0_13_pcaEigenValLoad);
    readFromCommaDelimitedFile_Float(cluster1_0_13_pcaEigenVecLoadPATH, cluster1_0_13_pcaEigenVecLoad);
    cv::PCA pcaANN1_0_13;
    pcaANN1_0_13.mean = cv::Mat(1, numVarANN, CV_32FC1, cluster1_0_13_pcaMeanLoad.data());
    pcaANN1_0_13.eigenvalues = cv::Mat(numComponentPCA,1, CV_32FC1, cluster1_0_13_pcaEigenValLoad.data());
    pcaANN1_0_13.eigenvectors = cv::Mat(numComponentPCA,numVarANN, CV_32FC1, cluster1_0_13_pcaEigenVecLoad.data());
    
    std::string cluster1_0_14_pcaMeanLoadPATH = casePath + cluster1_0_14_Folder + "PCA_mean.txt";
    std::string cluster1_0_14_pcaEigenValLoadPATH = casePath + cluster1_0_14_Folder + "PCA_eigenValues.txt";
    std::string cluster1_0_14_pcaEigenVecLoadPATH = casePath + cluster1_0_14_Folder + "PCA_eigenVectors.txt";
    vector<float> cluster1_0_14_pcaMeanLoad;
    vector<float> cluster1_0_14_pcaEigenValLoad;
    vector<float> cluster1_0_14_pcaEigenVecLoad;
    readFromCommaDelimitedFile_Float(cluster1_0_14_pcaMeanLoadPATH, cluster1_0_14_pcaMeanLoad);
    readFromCommaDelimitedFile_Float(cluster1_0_14_pcaEigenValLoadPATH, cluster1_0_14_pcaEigenValLoad);
    readFromCommaDelimitedFile_Float(cluster1_0_14_pcaEigenVecLoadPATH, cluster1_0_14_pcaEigenVecLoad);
    cv::PCA pcaANN1_0_14;
    pcaANN1_0_14.mean = cv::Mat(1, numVarANN, CV_32FC1, cluster1_0_14_pcaMeanLoad.data());
    pcaANN1_0_14.eigenvalues = cv::Mat(numComponentPCA,1, CV_32FC1, cluster1_0_14_pcaEigenValLoad.data());
    pcaANN1_0_14.eigenvectors = cv::Mat(numComponentPCA,numVarANN, CV_32FC1, cluster1_0_14_pcaEigenVecLoad.data());
    
    std::string cluster1_0_15_pcaMeanLoadPATH = casePath + cluster1_0_15_Folder + "PCA_mean.txt";
    std::string cluster1_0_15_pcaEigenValLoadPATH = casePath + cluster1_0_15_Folder + "PCA_eigenValues.txt";
    std::string cluster1_0_15_pcaEigenVecLoadPATH = casePath + cluster1_0_15_Folder + "PCA_eigenVectors.txt";
    vector<float> cluster1_0_15_pcaMeanLoad;
    vector<float> cluster1_0_15_pcaEigenValLoad;
    vector<float> cluster1_0_15_pcaEigenVecLoad;
    readFromCommaDelimitedFile_Float(cluster1_0_15_pcaMeanLoadPATH, cluster1_0_15_pcaMeanLoad);
    readFromCommaDelimitedFile_Float(cluster1_0_15_pcaEigenValLoadPATH, cluster1_0_15_pcaEigenValLoad);
    readFromCommaDelimitedFile_Float(cluster1_0_15_pcaEigenVecLoadPATH, cluster1_0_15_pcaEigenVecLoad);
    cv::PCA pcaANN1_0_15;
    pcaANN1_0_15.mean = cv::Mat(1, numVarANN, CV_32FC1, cluster1_0_15_pcaMeanLoad.data());
    pcaANN1_0_15.eigenvalues = cv::Mat(numComponentPCA,1, CV_32FC1, cluster1_0_15_pcaEigenValLoad.data());
    pcaANN1_0_15.eigenvectors = cv::Mat(numComponentPCA,numVarANN, CV_32FC1, cluster1_0_15_pcaEigenVecLoad.data());
    
    std::string cluster1_0_16_pcaMeanLoadPATH = casePath + cluster1_0_16_Folder + "PCA_mean.txt";
    std::string cluster1_0_16_pcaEigenValLoadPATH = casePath + cluster1_0_16_Folder + "PCA_eigenValues.txt";
    std::string cluster1_0_16_pcaEigenVecLoadPATH = casePath + cluster1_0_16_Folder + "PCA_eigenVectors.txt";
    vector<float> cluster1_0_16_pcaMeanLoad;
    vector<float> cluster1_0_16_pcaEigenValLoad;
    vector<float> cluster1_0_16_pcaEigenVecLoad;
    readFromCommaDelimitedFile_Float(cluster1_0_16_pcaMeanLoadPATH, cluster1_0_16_pcaMeanLoad);
    readFromCommaDelimitedFile_Float(cluster1_0_16_pcaEigenValLoadPATH, cluster1_0_16_pcaEigenValLoad);
    readFromCommaDelimitedFile_Float(cluster1_0_16_pcaEigenVecLoadPATH, cluster1_0_16_pcaEigenVecLoad);
    cv::PCA pcaANN1_0_16;
    pcaANN1_0_16.mean = cv::Mat(1, numVarANN, CV_32FC1, cluster1_0_16_pcaMeanLoad.data());
    pcaANN1_0_16.eigenvalues = cv::Mat(numComponentPCA,1, CV_32FC1, cluster1_0_16_pcaEigenValLoad.data());
    pcaANN1_0_16.eigenvectors = cv::Mat(numComponentPCA,numVarANN, CV_32FC1, cluster1_0_16_pcaEigenVecLoad.data());
    
    std::string cluster1_0_17_pcaMeanLoadPATH = casePath + cluster1_0_17_Folder + "PCA_mean.txt";
    std::string cluster1_0_17_pcaEigenValLoadPATH = casePath + cluster1_0_17_Folder + "PCA_eigenValues.txt";
    std::string cluster1_0_17_pcaEigenVecLoadPATH = casePath + cluster1_0_17_Folder + "PCA_eigenVectors.txt";
    vector<float> cluster1_0_17_pcaMeanLoad;
    vector<float> cluster1_0_17_pcaEigenValLoad;
    vector<float> cluster1_0_17_pcaEigenVecLoad;
    readFromCommaDelimitedFile_Float(cluster1_0_17_pcaMeanLoadPATH, cluster1_0_17_pcaMeanLoad);
    readFromCommaDelimitedFile_Float(cluster1_0_17_pcaEigenValLoadPATH, cluster1_0_17_pcaEigenValLoad);
    readFromCommaDelimitedFile_Float(cluster1_0_17_pcaEigenVecLoadPATH, cluster1_0_17_pcaEigenVecLoad);
    cv::PCA pcaANN1_0_17;
    pcaANN1_0_17.mean = cv::Mat(1, numVarANN, CV_32FC1, cluster1_0_17_pcaMeanLoad.data());
    pcaANN1_0_17.eigenvalues = cv::Mat(numComponentPCA,1, CV_32FC1, cluster1_0_17_pcaEigenValLoad.data());
    pcaANN1_0_17.eigenvectors = cv::Mat(numComponentPCA,numVarANN, CV_32FC1, cluster1_0_17_pcaEigenVecLoad.data());
    
    std::string cluster1_0_18_pcaMeanLoadPATH = casePath + cluster1_0_18_Folder + "PCA_mean.txt";
    std::string cluster1_0_18_pcaEigenValLoadPATH = casePath + cluster1_0_18_Folder + "PCA_eigenValues.txt";
    std::string cluster1_0_18_pcaEigenVecLoadPATH = casePath + cluster1_0_18_Folder + "PCA_eigenVectors.txt";
    vector<float> cluster1_0_18_pcaMeanLoad;
    vector<float> cluster1_0_18_pcaEigenValLoad;
    vector<float> cluster1_0_18_pcaEigenVecLoad;
    readFromCommaDelimitedFile_Float(cluster1_0_18_pcaMeanLoadPATH, cluster1_0_18_pcaMeanLoad);
    readFromCommaDelimitedFile_Float(cluster1_0_18_pcaEigenValLoadPATH, cluster1_0_18_pcaEigenValLoad);
    readFromCommaDelimitedFile_Float(cluster1_0_18_pcaEigenVecLoadPATH, cluster1_0_18_pcaEigenVecLoad);
    cv::PCA pcaANN1_0_18;
    pcaANN1_0_18.mean = cv::Mat(1, numVarANN, CV_32FC1, cluster1_0_18_pcaMeanLoad.data());
    pcaANN1_0_18.eigenvalues = cv::Mat(numComponentPCA,1, CV_32FC1, cluster1_0_18_pcaEigenValLoad.data());
    pcaANN1_0_18.eigenvectors = cv::Mat(numComponentPCA,numVarANN, CV_32FC1, cluster1_0_18_pcaEigenVecLoad.data());
    
    std::string cluster1_0_19_pcaMeanLoadPATH = casePath + cluster1_0_19_Folder + "PCA_mean.txt";
    std::string cluster1_0_19_pcaEigenValLoadPATH = casePath + cluster1_0_19_Folder + "PCA_eigenValues.txt";
    std::string cluster1_0_19_pcaEigenVecLoadPATH = casePath + cluster1_0_19_Folder + "PCA_eigenVectors.txt";
    vector<float> cluster1_0_19_pcaMeanLoad;
    vector<float> cluster1_0_19_pcaEigenValLoad;
    vector<float> cluster1_0_19_pcaEigenVecLoad;
    readFromCommaDelimitedFile_Float(cluster1_0_19_pcaMeanLoadPATH, cluster1_0_19_pcaMeanLoad);
    readFromCommaDelimitedFile_Float(cluster1_0_19_pcaEigenValLoadPATH, cluster1_0_19_pcaEigenValLoad);
    readFromCommaDelimitedFile_Float(cluster1_0_19_pcaEigenVecLoadPATH, cluster1_0_19_pcaEigenVecLoad);
    cv::PCA pcaANN1_0_19;
    pcaANN1_0_19.mean = cv::Mat(1, numVarANN, CV_32FC1, cluster1_0_19_pcaMeanLoad.data());
    pcaANN1_0_19.eigenvalues = cv::Mat(numComponentPCA,1, CV_32FC1, cluster1_0_19_pcaEigenValLoad.data());
    pcaANN1_0_19.eigenvectors = cv::Mat(numComponentPCA,numVarANN, CV_32FC1, cluster1_0_19_pcaEigenVecLoad.data());
    
    std::string cluster1_0_20_pcaMeanLoadPATH = casePath + cluster1_0_20_Folder + "PCA_mean.txt";
    std::string cluster1_0_20_pcaEigenValLoadPATH = casePath + cluster1_0_20_Folder + "PCA_eigenValues.txt";
    std::string cluster1_0_20_pcaEigenVecLoadPATH = casePath + cluster1_0_20_Folder + "PCA_eigenVectors.txt";
    vector<float> cluster1_0_20_pcaMeanLoad;
    vector<float> cluster1_0_20_pcaEigenValLoad;
    vector<float> cluster1_0_20_pcaEigenVecLoad;
    readFromCommaDelimitedFile_Float(cluster1_0_20_pcaMeanLoadPATH, cluster1_0_20_pcaMeanLoad);
    readFromCommaDelimitedFile_Float(cluster1_0_20_pcaEigenValLoadPATH, cluster1_0_20_pcaEigenValLoad);
    readFromCommaDelimitedFile_Float(cluster1_0_20_pcaEigenVecLoadPATH, cluster1_0_20_pcaEigenVecLoad);
    cv::PCA pcaANN1_0_20;
    pcaANN1_0_20.mean = cv::Mat(1, numVarANN, CV_32FC1, cluster1_0_20_pcaMeanLoad.data());
    pcaANN1_0_20.eigenvalues = cv::Mat(numComponentPCA,1, CV_32FC1, cluster1_0_20_pcaEigenValLoad.data());
    pcaANN1_0_20.eigenvectors = cv::Mat(numComponentPCA,numVarANN, CV_32FC1, cluster1_0_20_pcaEigenVecLoad.data());
    
    std::string cluster1_0_21_pcaMeanLoadPATH = casePath + cluster1_0_21_Folder + "PCA_mean.txt";
    std::string cluster1_0_21_pcaEigenValLoadPATH = casePath + cluster1_0_21_Folder + "PCA_eigenValues.txt";
    std::string cluster1_0_21_pcaEigenVecLoadPATH = casePath + cluster1_0_21_Folder + "PCA_eigenVectors.txt";
    vector<float> cluster1_0_21_pcaMeanLoad;
    vector<float> cluster1_0_21_pcaEigenValLoad;
    vector<float> cluster1_0_21_pcaEigenVecLoad;
    readFromCommaDelimitedFile_Float(cluster1_0_21_pcaMeanLoadPATH, cluster1_0_21_pcaMeanLoad);
    readFromCommaDelimitedFile_Float(cluster1_0_21_pcaEigenValLoadPATH, cluster1_0_21_pcaEigenValLoad);
    readFromCommaDelimitedFile_Float(cluster1_0_21_pcaEigenVecLoadPATH, cluster1_0_21_pcaEigenVecLoad);
    cv::PCA pcaANN1_0_21;
    pcaANN1_0_21.mean = cv::Mat(1, numVarANN, CV_32FC1, cluster1_0_21_pcaMeanLoad.data());
    pcaANN1_0_21.eigenvalues = cv::Mat(numComponentPCA,1, CV_32FC1, cluster1_0_21_pcaEigenValLoad.data());
    pcaANN1_0_21.eigenvectors = cv::Mat(numComponentPCA,numVarANN, CV_32FC1, cluster1_0_21_pcaEigenVecLoad.data());
    
    std::string cluster1_0_22_pcaMeanLoadPATH = casePath + cluster1_0_22_Folder + "PCA_mean.txt";
    std::string cluster1_0_22_pcaEigenValLoadPATH = casePath + cluster1_0_22_Folder + "PCA_eigenValues.txt";
    std::string cluster1_0_22_pcaEigenVecLoadPATH = casePath + cluster1_0_22_Folder + "PCA_eigenVectors.txt";
    vector<float> cluster1_0_22_pcaMeanLoad;
    vector<float> cluster1_0_22_pcaEigenValLoad;
    vector<float> cluster1_0_22_pcaEigenVecLoad;
    readFromCommaDelimitedFile_Float(cluster1_0_22_pcaMeanLoadPATH, cluster1_0_22_pcaMeanLoad);
    readFromCommaDelimitedFile_Float(cluster1_0_22_pcaEigenValLoadPATH, cluster1_0_22_pcaEigenValLoad);
    readFromCommaDelimitedFile_Float(cluster1_0_22_pcaEigenVecLoadPATH, cluster1_0_22_pcaEigenVecLoad);
    cv::PCA pcaANN1_0_22;
    pcaANN1_0_22.mean = cv::Mat(1, numVarANN, CV_32FC1, cluster1_0_22_pcaMeanLoad.data());
    pcaANN1_0_22.eigenvalues = cv::Mat(numComponentPCA,1, CV_32FC1, cluster1_0_22_pcaEigenValLoad.data());
    pcaANN1_0_22.eigenvectors = cv::Mat(numComponentPCA,numVarANN, CV_32FC1, cluster1_0_22_pcaEigenVecLoad.data());
    
    std::string cluster1_0_23_pcaMeanLoadPATH = casePath + cluster1_0_23_Folder + "PCA_mean.txt";
    std::string cluster1_0_23_pcaEigenValLoadPATH = casePath + cluster1_0_23_Folder + "PCA_eigenValues.txt";
    std::string cluster1_0_23_pcaEigenVecLoadPATH = casePath + cluster1_0_23_Folder + "PCA_eigenVectors.txt";
    vector<float> cluster1_0_23_pcaMeanLoad;
    vector<float> cluster1_0_23_pcaEigenValLoad;
    vector<float> cluster1_0_23_pcaEigenVecLoad;
    readFromCommaDelimitedFile_Float(cluster1_0_23_pcaMeanLoadPATH, cluster1_0_23_pcaMeanLoad);
    readFromCommaDelimitedFile_Float(cluster1_0_23_pcaEigenValLoadPATH, cluster1_0_23_pcaEigenValLoad);
    readFromCommaDelimitedFile_Float(cluster1_0_23_pcaEigenVecLoadPATH, cluster1_0_23_pcaEigenVecLoad);
    cv::PCA pcaANN1_0_23;
    pcaANN1_0_23.mean = cv::Mat(1, numVarANN, CV_32FC1, cluster1_0_23_pcaMeanLoad.data());
    pcaANN1_0_23.eigenvalues = cv::Mat(numComponentPCA,1, CV_32FC1, cluster1_0_23_pcaEigenValLoad.data());
    pcaANN1_0_23.eigenvectors = cv::Mat(numComponentPCA,numVarANN, CV_32FC1, cluster1_0_23_pcaEigenVecLoad.data());
    
    std::string cluster1_0_24_pcaMeanLoadPATH = casePath + cluster1_0_24_Folder + "PCA_mean.txt";
    std::string cluster1_0_24_pcaEigenValLoadPATH = casePath + cluster1_0_24_Folder + "PCA_eigenValues.txt";
    std::string cluster1_0_24_pcaEigenVecLoadPATH = casePath + cluster1_0_24_Folder + "PCA_eigenVectors.txt";
    vector<float> cluster1_0_24_pcaMeanLoad;
    vector<float> cluster1_0_24_pcaEigenValLoad;
    vector<float> cluster1_0_24_pcaEigenVecLoad;
    readFromCommaDelimitedFile_Float(cluster1_0_24_pcaMeanLoadPATH, cluster1_0_24_pcaMeanLoad);
    readFromCommaDelimitedFile_Float(cluster1_0_24_pcaEigenValLoadPATH, cluster1_0_24_pcaEigenValLoad);
    readFromCommaDelimitedFile_Float(cluster1_0_24_pcaEigenVecLoadPATH, cluster1_0_24_pcaEigenVecLoad);
    cv::PCA pcaANN1_0_24;
    pcaANN1_0_24.mean = cv::Mat(1, numVarANN, CV_32FC1, cluster1_0_24_pcaMeanLoad.data());
    pcaANN1_0_24.eigenvalues = cv::Mat(numComponentPCA,1, CV_32FC1, cluster1_0_24_pcaEigenValLoad.data());
    pcaANN1_0_24.eigenvectors = cv::Mat(numComponentPCA,numVarANN, CV_32FC1, cluster1_0_24_pcaEigenVecLoad.data());
    
    std::string cluster1_0_25_pcaMeanLoadPATH = casePath + cluster1_0_25_Folder + "PCA_mean.txt";
    std::string cluster1_0_25_pcaEigenValLoadPATH = casePath + cluster1_0_25_Folder + "PCA_eigenValues.txt";
    std::string cluster1_0_25_pcaEigenVecLoadPATH = casePath + cluster1_0_25_Folder + "PCA_eigenVectors.txt";
    vector<float> cluster1_0_25_pcaMeanLoad;
    vector<float> cluster1_0_25_pcaEigenValLoad;
    vector<float> cluster1_0_25_pcaEigenVecLoad;
    readFromCommaDelimitedFile_Float(cluster1_0_25_pcaMeanLoadPATH, cluster1_0_25_pcaMeanLoad);
    readFromCommaDelimitedFile_Float(cluster1_0_25_pcaEigenValLoadPATH, cluster1_0_25_pcaEigenValLoad);
    readFromCommaDelimitedFile_Float(cluster1_0_25_pcaEigenVecLoadPATH, cluster1_0_25_pcaEigenVecLoad);
    cv::PCA pcaANN1_0_25;
    pcaANN1_0_25.mean = cv::Mat(1, numVarANN, CV_32FC1, cluster1_0_25_pcaMeanLoad.data());
    pcaANN1_0_25.eigenvalues = cv::Mat(numComponentPCA,1, CV_32FC1, cluster1_0_25_pcaEigenValLoad.data());
    pcaANN1_0_25.eigenvectors = cv::Mat(numComponentPCA,numVarANN, CV_32FC1, cluster1_0_25_pcaEigenVecLoad.data());
    
    std::string cluster1_0_26_pcaMeanLoadPATH = casePath + cluster1_0_26_Folder + "PCA_mean.txt";
    std::string cluster1_0_26_pcaEigenValLoadPATH = casePath + cluster1_0_26_Folder + "PCA_eigenValues.txt";
    std::string cluster1_0_26_pcaEigenVecLoadPATH = casePath + cluster1_0_26_Folder + "PCA_eigenVectors.txt";
    vector<float> cluster1_0_26_pcaMeanLoad;
    vector<float> cluster1_0_26_pcaEigenValLoad;
    vector<float> cluster1_0_26_pcaEigenVecLoad;
    readFromCommaDelimitedFile_Float(cluster1_0_26_pcaMeanLoadPATH, cluster1_0_26_pcaMeanLoad);
    readFromCommaDelimitedFile_Float(cluster1_0_26_pcaEigenValLoadPATH, cluster1_0_26_pcaEigenValLoad);
    readFromCommaDelimitedFile_Float(cluster1_0_26_pcaEigenVecLoadPATH, cluster1_0_26_pcaEigenVecLoad);
    cv::PCA pcaANN1_0_26;
    pcaANN1_0_26.mean = cv::Mat(1, numVarANN, CV_32FC1, cluster1_0_26_pcaMeanLoad.data());
    pcaANN1_0_26.eigenvalues = cv::Mat(numComponentPCA,1, CV_32FC1, cluster1_0_26_pcaEigenValLoad.data());
    pcaANN1_0_26.eigenvectors = cv::Mat(numComponentPCA,numVarANN, CV_32FC1, cluster1_0_26_pcaEigenVecLoad.data());
    
    std::string cluster1_0_27_pcaMeanLoadPATH = casePath + cluster1_0_27_Folder + "PCA_mean.txt";
    std::string cluster1_0_27_pcaEigenValLoadPATH = casePath + cluster1_0_27_Folder + "PCA_eigenValues.txt";
    std::string cluster1_0_27_pcaEigenVecLoadPATH = casePath + cluster1_0_27_Folder + "PCA_eigenVectors.txt";
    vector<float> cluster1_0_27_pcaMeanLoad;
    vector<float> cluster1_0_27_pcaEigenValLoad;
    vector<float> cluster1_0_27_pcaEigenVecLoad;
    readFromCommaDelimitedFile_Float(cluster1_0_27_pcaMeanLoadPATH, cluster1_0_27_pcaMeanLoad);
    readFromCommaDelimitedFile_Float(cluster1_0_27_pcaEigenValLoadPATH, cluster1_0_27_pcaEigenValLoad);
    readFromCommaDelimitedFile_Float(cluster1_0_27_pcaEigenVecLoadPATH, cluster1_0_27_pcaEigenVecLoad);
    cv::PCA pcaANN1_0_27;
    pcaANN1_0_27.mean = cv::Mat(1, numVarANN, CV_32FC1, cluster1_0_27_pcaMeanLoad.data());
    pcaANN1_0_27.eigenvalues = cv::Mat(numComponentPCA,1, CV_32FC1, cluster1_0_27_pcaEigenValLoad.data());
    pcaANN1_0_27.eigenvectors = cv::Mat(numComponentPCA,numVarANN, CV_32FC1, cluster1_0_27_pcaEigenVecLoad.data());
    
    std::string cluster1_0_28_pcaMeanLoadPATH = casePath + cluster1_0_28_Folder + "PCA_mean.txt";
    std::string cluster1_0_28_pcaEigenValLoadPATH = casePath + cluster1_0_28_Folder + "PCA_eigenValues.txt";
    std::string cluster1_0_28_pcaEigenVecLoadPATH = casePath + cluster1_0_28_Folder + "PCA_eigenVectors.txt";
    vector<float> cluster1_0_28_pcaMeanLoad;
    vector<float> cluster1_0_28_pcaEigenValLoad;
    vector<float> cluster1_0_28_pcaEigenVecLoad;
    readFromCommaDelimitedFile_Float(cluster1_0_28_pcaMeanLoadPATH, cluster1_0_28_pcaMeanLoad);
    readFromCommaDelimitedFile_Float(cluster1_0_28_pcaEigenValLoadPATH, cluster1_0_28_pcaEigenValLoad);
    readFromCommaDelimitedFile_Float(cluster1_0_28_pcaEigenVecLoadPATH, cluster1_0_28_pcaEigenVecLoad);
    cv::PCA pcaANN1_0_28;
    pcaANN1_0_28.mean = cv::Mat(1, numVarANN, CV_32FC1, cluster1_0_28_pcaMeanLoad.data());
    pcaANN1_0_28.eigenvalues = cv::Mat(numComponentPCA,1, CV_32FC1, cluster1_0_28_pcaEigenValLoad.data());
    pcaANN1_0_28.eigenvectors = cv::Mat(numComponentPCA,numVarANN, CV_32FC1, cluster1_0_28_pcaEigenVecLoad.data());
    
    std::string cluster1_0_29_pcaMeanLoadPATH = casePath + cluster1_0_29_Folder + "PCA_mean.txt";
    std::string cluster1_0_29_pcaEigenValLoadPATH = casePath + cluster1_0_29_Folder + "PCA_eigenValues.txt";
    std::string cluster1_0_29_pcaEigenVecLoadPATH = casePath + cluster1_0_29_Folder + "PCA_eigenVectors.txt";
    vector<float> cluster1_0_29_pcaMeanLoad;
    vector<float> cluster1_0_29_pcaEigenValLoad;
    vector<float> cluster1_0_29_pcaEigenVecLoad;
    readFromCommaDelimitedFile_Float(cluster1_0_29_pcaMeanLoadPATH, cluster1_0_29_pcaMeanLoad);
    readFromCommaDelimitedFile_Float(cluster1_0_29_pcaEigenValLoadPATH, cluster1_0_29_pcaEigenValLoad);
    readFromCommaDelimitedFile_Float(cluster1_0_29_pcaEigenVecLoadPATH, cluster1_0_29_pcaEigenVecLoad);
    cv::PCA pcaANN1_0_29;
    pcaANN1_0_29.mean = cv::Mat(1, numVarANN, CV_32FC1, cluster1_0_29_pcaMeanLoad.data());
    pcaANN1_0_29.eigenvalues = cv::Mat(numComponentPCA,1, CV_32FC1, cluster1_0_29_pcaEigenValLoad.data());
    pcaANN1_0_29.eigenvectors = cv::Mat(numComponentPCA,numVarANN, CV_32FC1, cluster1_0_29_pcaEigenVecLoad.data());
    
    std::string cluster1_0_30_pcaMeanLoadPATH = casePath + cluster1_0_30_Folder + "PCA_mean.txt";
    std::string cluster1_0_30_pcaEigenValLoadPATH = casePath + cluster1_0_30_Folder + "PCA_eigenValues.txt";
    std::string cluster1_0_30_pcaEigenVecLoadPATH = casePath + cluster1_0_30_Folder + "PCA_eigenVectors.txt";
    vector<float> cluster1_0_30_pcaMeanLoad;
    vector<float> cluster1_0_30_pcaEigenValLoad;
    vector<float> cluster1_0_30_pcaEigenVecLoad;
    readFromCommaDelimitedFile_Float(cluster1_0_30_pcaMeanLoadPATH, cluster1_0_30_pcaMeanLoad);
    readFromCommaDelimitedFile_Float(cluster1_0_30_pcaEigenValLoadPATH, cluster1_0_30_pcaEigenValLoad);
    readFromCommaDelimitedFile_Float(cluster1_0_30_pcaEigenVecLoadPATH, cluster1_0_30_pcaEigenVecLoad);
    cv::PCA pcaANN1_0_30;
    pcaANN1_0_30.mean = cv::Mat(1, numVarANN, CV_32FC1, cluster1_0_30_pcaMeanLoad.data());
    pcaANN1_0_30.eigenvalues = cv::Mat(numComponentPCA,1, CV_32FC1, cluster1_0_30_pcaEigenValLoad.data());
    pcaANN1_0_30.eigenvectors = cv::Mat(numComponentPCA,numVarANN, CV_32FC1, cluster1_0_30_pcaEigenVecLoad.data());
    
    std::string cluster1_0_31_pcaMeanLoadPATH = casePath + cluster1_0_31_Folder + "PCA_mean.txt";
    std::string cluster1_0_31_pcaEigenValLoadPATH = casePath + cluster1_0_31_Folder + "PCA_eigenValues.txt";
    std::string cluster1_0_31_pcaEigenVecLoadPATH = casePath + cluster1_0_31_Folder + "PCA_eigenVectors.txt";
    vector<float> cluster1_0_31_pcaMeanLoad;
    vector<float> cluster1_0_31_pcaEigenValLoad;
    vector<float> cluster1_0_31_pcaEigenVecLoad;
    readFromCommaDelimitedFile_Float(cluster1_0_31_pcaMeanLoadPATH, cluster1_0_31_pcaMeanLoad);
    readFromCommaDelimitedFile_Float(cluster1_0_31_pcaEigenValLoadPATH, cluster1_0_31_pcaEigenValLoad);
    readFromCommaDelimitedFile_Float(cluster1_0_31_pcaEigenVecLoadPATH, cluster1_0_31_pcaEigenVecLoad);
    cv::PCA pcaANN1_0_31;
    pcaANN1_0_31.mean = cv::Mat(1, numVarANN, CV_32FC1, cluster1_0_31_pcaMeanLoad.data());
    pcaANN1_0_31.eigenvalues = cv::Mat(numComponentPCA,1, CV_32FC1, cluster1_0_31_pcaEigenValLoad.data());
    pcaANN1_0_31.eigenvectors = cv::Mat(numComponentPCA,numVarANN, CV_32FC1, cluster1_0_31_pcaEigenVecLoad.data());
    
    std::string cluster1_0_32_pcaMeanLoadPATH = casePath + cluster1_0_32_Folder + "PCA_mean.txt";
    std::string cluster1_0_32_pcaEigenValLoadPATH = casePath + cluster1_0_32_Folder + "PCA_eigenValues.txt";
    std::string cluster1_0_32_pcaEigenVecLoadPATH = casePath + cluster1_0_32_Folder + "PCA_eigenVectors.txt";
    vector<float> cluster1_0_32_pcaMeanLoad;
    vector<float> cluster1_0_32_pcaEigenValLoad;
    vector<float> cluster1_0_32_pcaEigenVecLoad;
    readFromCommaDelimitedFile_Float(cluster1_0_32_pcaMeanLoadPATH, cluster1_0_32_pcaMeanLoad);
    readFromCommaDelimitedFile_Float(cluster1_0_32_pcaEigenValLoadPATH, cluster1_0_32_pcaEigenValLoad);
    readFromCommaDelimitedFile_Float(cluster1_0_32_pcaEigenVecLoadPATH, cluster1_0_32_pcaEigenVecLoad);
    cv::PCA pcaANN1_0_32;
    pcaANN1_0_32.mean = cv::Mat(1, numVarANN, CV_32FC1, cluster1_0_32_pcaMeanLoad.data());
    pcaANN1_0_32.eigenvalues = cv::Mat(numComponentPCA,1, CV_32FC1, cluster1_0_32_pcaEigenValLoad.data());
    pcaANN1_0_32.eigenvectors = cv::Mat(numComponentPCA,numVarANN, CV_32FC1, cluster1_0_32_pcaEigenVecLoad.data());
    
    std::string cluster1_0_33_pcaMeanLoadPATH = casePath + cluster1_0_33_Folder + "PCA_mean.txt";
    std::string cluster1_0_33_pcaEigenValLoadPATH = casePath + cluster1_0_33_Folder + "PCA_eigenValues.txt";
    std::string cluster1_0_33_pcaEigenVecLoadPATH = casePath + cluster1_0_33_Folder + "PCA_eigenVectors.txt";
    vector<float> cluster1_0_33_pcaMeanLoad;
    vector<float> cluster1_0_33_pcaEigenValLoad;
    vector<float> cluster1_0_33_pcaEigenVecLoad;
    readFromCommaDelimitedFile_Float(cluster1_0_33_pcaMeanLoadPATH, cluster1_0_33_pcaMeanLoad);
    readFromCommaDelimitedFile_Float(cluster1_0_33_pcaEigenValLoadPATH, cluster1_0_33_pcaEigenValLoad);
    readFromCommaDelimitedFile_Float(cluster1_0_33_pcaEigenVecLoadPATH, cluster1_0_33_pcaEigenVecLoad);
    cv::PCA pcaANN1_0_33;
    pcaANN1_0_33.mean = cv::Mat(1, numVarANN, CV_32FC1, cluster1_0_33_pcaMeanLoad.data());
    pcaANN1_0_33.eigenvalues = cv::Mat(numComponentPCA,1, CV_32FC1, cluster1_0_33_pcaEigenValLoad.data());
    pcaANN1_0_33.eigenvectors = cv::Mat(numComponentPCA,numVarANN, CV_32FC1, cluster1_0_33_pcaEigenVecLoad.data());
    
    std::string cluster1_0_34_pcaMeanLoadPATH = casePath + cluster1_0_34_Folder + "PCA_mean.txt";
    std::string cluster1_0_34_pcaEigenValLoadPATH = casePath + cluster1_0_34_Folder + "PCA_eigenValues.txt";
    std::string cluster1_0_34_pcaEigenVecLoadPATH = casePath + cluster1_0_34_Folder + "PCA_eigenVectors.txt";
    vector<float> cluster1_0_34_pcaMeanLoad;
    vector<float> cluster1_0_34_pcaEigenValLoad;
    vector<float> cluster1_0_34_pcaEigenVecLoad;
    readFromCommaDelimitedFile_Float(cluster1_0_34_pcaMeanLoadPATH, cluster1_0_34_pcaMeanLoad);
    readFromCommaDelimitedFile_Float(cluster1_0_34_pcaEigenValLoadPATH, cluster1_0_34_pcaEigenValLoad);
    readFromCommaDelimitedFile_Float(cluster1_0_34_pcaEigenVecLoadPATH, cluster1_0_34_pcaEigenVecLoad);
    cv::PCA pcaANN1_0_34;
    pcaANN1_0_34.mean = cv::Mat(1, numVarANN, CV_32FC1, cluster1_0_34_pcaMeanLoad.data());
    pcaANN1_0_34.eigenvalues = cv::Mat(numComponentPCA,1, CV_32FC1, cluster1_0_34_pcaEigenValLoad.data());
    pcaANN1_0_34.eigenvectors = cv::Mat(numComponentPCA,numVarANN, CV_32FC1, cluster1_0_34_pcaEigenVecLoad.data());
    
    std::string cluster1_0_35_pcaMeanLoadPATH = casePath + cluster1_0_35_Folder + "PCA_mean.txt";
    std::string cluster1_0_35_pcaEigenValLoadPATH = casePath + cluster1_0_35_Folder + "PCA_eigenValues.txt";
    std::string cluster1_0_35_pcaEigenVecLoadPATH = casePath + cluster1_0_35_Folder + "PCA_eigenVectors.txt";
    vector<float> cluster1_0_35_pcaMeanLoad;
    vector<float> cluster1_0_35_pcaEigenValLoad;
    vector<float> cluster1_0_35_pcaEigenVecLoad;
    readFromCommaDelimitedFile_Float(cluster1_0_35_pcaMeanLoadPATH, cluster1_0_35_pcaMeanLoad);
    readFromCommaDelimitedFile_Float(cluster1_0_35_pcaEigenValLoadPATH, cluster1_0_35_pcaEigenValLoad);
    readFromCommaDelimitedFile_Float(cluster1_0_35_pcaEigenVecLoadPATH, cluster1_0_35_pcaEigenVecLoad);
    cv::PCA pcaANN1_0_35;
    pcaANN1_0_35.mean = cv::Mat(1, numVarANN, CV_32FC1, cluster1_0_35_pcaMeanLoad.data());
    pcaANN1_0_35.eigenvalues = cv::Mat(numComponentPCA,1, CV_32FC1, cluster1_0_35_pcaEigenValLoad.data());
    pcaANN1_0_35.eigenvectors = cv::Mat(numComponentPCA,numVarANN, CV_32FC1, cluster1_0_35_pcaEigenVecLoad.data());
    
    std::string cluster1_0_36_pcaMeanLoadPATH = casePath + cluster1_0_36_Folder + "PCA_mean.txt";
    std::string cluster1_0_36_pcaEigenValLoadPATH = casePath + cluster1_0_36_Folder + "PCA_eigenValues.txt";
    std::string cluster1_0_36_pcaEigenVecLoadPATH = casePath + cluster1_0_36_Folder + "PCA_eigenVectors.txt";
    vector<float> cluster1_0_36_pcaMeanLoad;
    vector<float> cluster1_0_36_pcaEigenValLoad;
    vector<float> cluster1_0_36_pcaEigenVecLoad;
    readFromCommaDelimitedFile_Float(cluster1_0_36_pcaMeanLoadPATH, cluster1_0_36_pcaMeanLoad);
    readFromCommaDelimitedFile_Float(cluster1_0_36_pcaEigenValLoadPATH, cluster1_0_36_pcaEigenValLoad);
    readFromCommaDelimitedFile_Float(cluster1_0_36_pcaEigenVecLoadPATH, cluster1_0_36_pcaEigenVecLoad);
    cv::PCA pcaANN1_0_36;
    pcaANN1_0_36.mean = cv::Mat(1, numVarANN, CV_32FC1, cluster1_0_36_pcaMeanLoad.data());
    pcaANN1_0_36.eigenvalues = cv::Mat(numComponentPCA,1, CV_32FC1, cluster1_0_36_pcaEigenValLoad.data());
    pcaANN1_0_36.eigenvectors = cv::Mat(numComponentPCA,numVarANN, CV_32FC1, cluster1_0_36_pcaEigenVecLoad.data());
    
    std::string cluster1_0_37_pcaMeanLoadPATH = casePath + cluster1_0_37_Folder + "PCA_mean.txt";
    std::string cluster1_0_37_pcaEigenValLoadPATH = casePath + cluster1_0_37_Folder + "PCA_eigenValues.txt";
    std::string cluster1_0_37_pcaEigenVecLoadPATH = casePath + cluster1_0_37_Folder + "PCA_eigenVectors.txt";
    vector<float> cluster1_0_37_pcaMeanLoad;
    vector<float> cluster1_0_37_pcaEigenValLoad;
    vector<float> cluster1_0_37_pcaEigenVecLoad;
    readFromCommaDelimitedFile_Float(cluster1_0_37_pcaMeanLoadPATH, cluster1_0_37_pcaMeanLoad);
    readFromCommaDelimitedFile_Float(cluster1_0_37_pcaEigenValLoadPATH, cluster1_0_37_pcaEigenValLoad);
    readFromCommaDelimitedFile_Float(cluster1_0_37_pcaEigenVecLoadPATH, cluster1_0_37_pcaEigenVecLoad);
    cv::PCA pcaANN1_0_37;
    pcaANN1_0_37.mean = cv::Mat(1, numVarANN, CV_32FC1, cluster1_0_37_pcaMeanLoad.data());
    pcaANN1_0_37.eigenvalues = cv::Mat(numComponentPCA,1, CV_32FC1, cluster1_0_37_pcaEigenValLoad.data());
    pcaANN1_0_37.eigenvectors = cv::Mat(numComponentPCA,numVarANN, CV_32FC1, cluster1_0_37_pcaEigenVecLoad.data());
    
    std::string cluster1_0_38_pcaMeanLoadPATH = casePath + cluster1_0_38_Folder + "PCA_mean.txt";
    std::string cluster1_0_38_pcaEigenValLoadPATH = casePath + cluster1_0_38_Folder + "PCA_eigenValues.txt";
    std::string cluster1_0_38_pcaEigenVecLoadPATH = casePath + cluster1_0_38_Folder + "PCA_eigenVectors.txt";
    vector<float> cluster1_0_38_pcaMeanLoad;
    vector<float> cluster1_0_38_pcaEigenValLoad;
    vector<float> cluster1_0_38_pcaEigenVecLoad;
    readFromCommaDelimitedFile_Float(cluster1_0_38_pcaMeanLoadPATH, cluster1_0_38_pcaMeanLoad);
    readFromCommaDelimitedFile_Float(cluster1_0_38_pcaEigenValLoadPATH, cluster1_0_38_pcaEigenValLoad);
    readFromCommaDelimitedFile_Float(cluster1_0_38_pcaEigenVecLoadPATH, cluster1_0_38_pcaEigenVecLoad);
    cv::PCA pcaANN1_0_38;
    pcaANN1_0_38.mean = cv::Mat(1, numVarANN, CV_32FC1, cluster1_0_38_pcaMeanLoad.data());
    pcaANN1_0_38.eigenvalues = cv::Mat(numComponentPCA,1, CV_32FC1, cluster1_0_38_pcaEigenValLoad.data());
    pcaANN1_0_38.eigenvectors = cv::Mat(numComponentPCA,numVarANN, CV_32FC1, cluster1_0_38_pcaEigenVecLoad.data());
    
    std::string cluster1_0_39_pcaMeanLoadPATH = casePath + cluster1_0_39_Folder + "PCA_mean.txt";
    std::string cluster1_0_39_pcaEigenValLoadPATH = casePath + cluster1_0_39_Folder + "PCA_eigenValues.txt";
    std::string cluster1_0_39_pcaEigenVecLoadPATH = casePath + cluster1_0_39_Folder + "PCA_eigenVectors.txt";
    vector<float> cluster1_0_39_pcaMeanLoad;
    vector<float> cluster1_0_39_pcaEigenValLoad;
    vector<float> cluster1_0_39_pcaEigenVecLoad;
    readFromCommaDelimitedFile_Float(cluster1_0_39_pcaMeanLoadPATH, cluster1_0_39_pcaMeanLoad);
    readFromCommaDelimitedFile_Float(cluster1_0_39_pcaEigenValLoadPATH, cluster1_0_39_pcaEigenValLoad);
    readFromCommaDelimitedFile_Float(cluster1_0_39_pcaEigenVecLoadPATH, cluster1_0_39_pcaEigenVecLoad);
    cv::PCA pcaANN1_0_39;
    pcaANN1_0_39.mean = cv::Mat(1, numVarANN, CV_32FC1, cluster1_0_39_pcaMeanLoad.data());
    pcaANN1_0_39.eigenvalues = cv::Mat(numComponentPCA,1, CV_32FC1, cluster1_0_39_pcaEigenValLoad.data());
    pcaANN1_0_39.eigenvectors = cv::Mat(numComponentPCA,numVarANN, CV_32FC1, cluster1_0_39_pcaEigenVecLoad.data());

    // GrandChild of parentCluster1 - childCluster1 - grandChild
    std::string cluster1_1_0_pcaMeanLoadPATH = casePath + cluster1_1_0_Folder + "PCA_mean.txt";
    std::string cluster1_1_0_pcaEigenValLoadPATH = casePath + cluster1_1_0_Folder + "PCA_eigenValues.txt";
    std::string cluster1_1_0_pcaEigenVecLoadPATH = casePath + cluster1_1_0_Folder + "PCA_eigenVectors.txt";
    vector<float> cluster1_1_0_pcaMeanLoad;
    vector<float> cluster1_1_0_pcaEigenValLoad;
    vector<float> cluster1_1_0_pcaEigenVecLoad;
    readFromCommaDelimitedFile_Float(cluster1_1_0_pcaMeanLoadPATH, cluster1_1_0_pcaMeanLoad);
    readFromCommaDelimitedFile_Float(cluster1_1_0_pcaEigenValLoadPATH, cluster1_1_0_pcaEigenValLoad);
    readFromCommaDelimitedFile_Float(cluster1_1_0_pcaEigenVecLoadPATH, cluster1_1_0_pcaEigenVecLoad);
    cv::PCA pcaANN1_1_0;
    pcaANN1_1_0.mean = cv::Mat(1, numVarANN, CV_32FC1, cluster1_1_0_pcaMeanLoad.data());
    pcaANN1_1_0.eigenvalues = cv::Mat(numComponentPCA,1, CV_32FC1, cluster1_1_0_pcaEigenValLoad.data());
    pcaANN1_1_0.eigenvectors = cv::Mat(numComponentPCA,numVarANN, CV_32FC1, cluster1_1_0_pcaEigenVecLoad.data());
    
    std::string cluster1_1_1_pcaMeanLoadPATH = casePath + cluster1_1_1_Folder + "PCA_mean.txt";
    std::string cluster1_1_1_pcaEigenValLoadPATH = casePath + cluster1_1_1_Folder + "PCA_eigenValues.txt";
    std::string cluster1_1_1_pcaEigenVecLoadPATH = casePath + cluster1_1_1_Folder + "PCA_eigenVectors.txt";
    vector<float> cluster1_1_1_pcaMeanLoad;
    vector<float> cluster1_1_1_pcaEigenValLoad;
    vector<float> cluster1_1_1_pcaEigenVecLoad;
    readFromCommaDelimitedFile_Float(cluster1_1_1_pcaMeanLoadPATH, cluster1_1_1_pcaMeanLoad);
    readFromCommaDelimitedFile_Float(cluster1_1_1_pcaEigenValLoadPATH, cluster1_1_1_pcaEigenValLoad);
    readFromCommaDelimitedFile_Float(cluster1_1_1_pcaEigenVecLoadPATH, cluster1_1_1_pcaEigenVecLoad);
    cv::PCA pcaANN1_1_1;
    pcaANN1_1_1.mean = cv::Mat(1, numVarANN, CV_32FC1, cluster1_1_1_pcaMeanLoad.data());
    pcaANN1_1_1.eigenvalues = cv::Mat(numComponentPCA,1, CV_32FC1, cluster1_1_1_pcaEigenValLoad.data());
    pcaANN1_1_1.eigenvectors = cv::Mat(numComponentPCA,numVarANN, CV_32FC1, cluster1_1_1_pcaEigenVecLoad.data());
    
    std::string cluster1_1_2_pcaMeanLoadPATH = casePath + cluster1_1_2_Folder + "PCA_mean.txt";
    std::string cluster1_1_2_pcaEigenValLoadPATH = casePath + cluster1_1_2_Folder + "PCA_eigenValues.txt";
    std::string cluster1_1_2_pcaEigenVecLoadPATH = casePath + cluster1_1_2_Folder + "PCA_eigenVectors.txt";
    vector<float> cluster1_1_2_pcaMeanLoad;
    vector<float> cluster1_1_2_pcaEigenValLoad;
    vector<float> cluster1_1_2_pcaEigenVecLoad;
    readFromCommaDelimitedFile_Float(cluster1_1_2_pcaMeanLoadPATH, cluster1_1_2_pcaMeanLoad);
    readFromCommaDelimitedFile_Float(cluster1_1_2_pcaEigenValLoadPATH, cluster1_1_2_pcaEigenValLoad);
    readFromCommaDelimitedFile_Float(cluster1_1_2_pcaEigenVecLoadPATH, cluster1_1_2_pcaEigenVecLoad);
    cv::PCA pcaANN1_1_2;
    pcaANN1_1_2.mean = cv::Mat(1, numVarANN, CV_32FC1, cluster1_1_2_pcaMeanLoad.data());
    pcaANN1_1_2.eigenvalues = cv::Mat(numComponentPCA,1, CV_32FC1, cluster1_1_2_pcaEigenValLoad.data());
    pcaANN1_1_2.eigenvectors = cv::Mat(numComponentPCA,numVarANN, CV_32FC1, cluster1_1_2_pcaEigenVecLoad.data());
    
    std::string cluster1_1_3_pcaMeanLoadPATH = casePath + cluster1_1_3_Folder + "PCA_mean.txt";
    std::string cluster1_1_3_pcaEigenValLoadPATH = casePath + cluster1_1_3_Folder + "PCA_eigenValues.txt";
    std::string cluster1_1_3_pcaEigenVecLoadPATH = casePath + cluster1_1_3_Folder + "PCA_eigenVectors.txt";
    vector<float> cluster1_1_3_pcaMeanLoad;
    vector<float> cluster1_1_3_pcaEigenValLoad;
    vector<float> cluster1_1_3_pcaEigenVecLoad;
    readFromCommaDelimitedFile_Float(cluster1_1_3_pcaMeanLoadPATH, cluster1_1_3_pcaMeanLoad);
    readFromCommaDelimitedFile_Float(cluster1_1_3_pcaEigenValLoadPATH, cluster1_1_3_pcaEigenValLoad);
    readFromCommaDelimitedFile_Float(cluster1_1_3_pcaEigenVecLoadPATH, cluster1_1_3_pcaEigenVecLoad);
    cv::PCA pcaANN1_1_3;
    pcaANN1_1_3.mean = cv::Mat(1, numVarANN, CV_32FC1, cluster1_1_3_pcaMeanLoad.data());
    pcaANN1_1_3.eigenvalues = cv::Mat(numComponentPCA,1, CV_32FC1, cluster1_1_3_pcaEigenValLoad.data());
    pcaANN1_1_3.eigenvectors = cv::Mat(numComponentPCA,numVarANN, CV_32FC1, cluster1_1_3_pcaEigenVecLoad.data());
    
    std::string cluster1_1_4_pcaMeanLoadPATH = casePath + cluster1_1_4_Folder + "PCA_mean.txt";
    std::string cluster1_1_4_pcaEigenValLoadPATH = casePath + cluster1_1_4_Folder + "PCA_eigenValues.txt";
    std::string cluster1_1_4_pcaEigenVecLoadPATH = casePath + cluster1_1_4_Folder + "PCA_eigenVectors.txt";
    vector<float> cluster1_1_4_pcaMeanLoad;
    vector<float> cluster1_1_4_pcaEigenValLoad;
    vector<float> cluster1_1_4_pcaEigenVecLoad;
    readFromCommaDelimitedFile_Float(cluster1_1_4_pcaMeanLoadPATH, cluster1_1_4_pcaMeanLoad);
    readFromCommaDelimitedFile_Float(cluster1_1_4_pcaEigenValLoadPATH, cluster1_1_4_pcaEigenValLoad);
    readFromCommaDelimitedFile_Float(cluster1_1_4_pcaEigenVecLoadPATH, cluster1_1_4_pcaEigenVecLoad);
    cv::PCA pcaANN1_1_4;
    pcaANN1_1_4.mean = cv::Mat(1, numVarANN, CV_32FC1, cluster1_1_4_pcaMeanLoad.data());
    pcaANN1_1_4.eigenvalues = cv::Mat(numComponentPCA,1, CV_32FC1, cluster1_1_4_pcaEigenValLoad.data());
    pcaANN1_1_4.eigenvectors = cv::Mat(numComponentPCA,numVarANN, CV_32FC1, cluster1_1_4_pcaEigenVecLoad.data());
    
    std::string cluster1_1_5_pcaMeanLoadPATH = casePath + cluster1_1_5_Folder + "PCA_mean.txt";
    std::string cluster1_1_5_pcaEigenValLoadPATH = casePath + cluster1_1_5_Folder + "PCA_eigenValues.txt";
    std::string cluster1_1_5_pcaEigenVecLoadPATH = casePath + cluster1_1_5_Folder + "PCA_eigenVectors.txt";
    vector<float> cluster1_1_5_pcaMeanLoad;
    vector<float> cluster1_1_5_pcaEigenValLoad;
    vector<float> cluster1_1_5_pcaEigenVecLoad;
    readFromCommaDelimitedFile_Float(cluster1_1_5_pcaMeanLoadPATH, cluster1_1_5_pcaMeanLoad);
    readFromCommaDelimitedFile_Float(cluster1_1_5_pcaEigenValLoadPATH, cluster1_1_5_pcaEigenValLoad);
    readFromCommaDelimitedFile_Float(cluster1_1_5_pcaEigenVecLoadPATH, cluster1_1_5_pcaEigenVecLoad);
    cv::PCA pcaANN1_1_5;
    pcaANN1_1_5.mean = cv::Mat(1, numVarANN, CV_32FC1, cluster1_1_5_pcaMeanLoad.data());
    pcaANN1_1_5.eigenvalues = cv::Mat(numComponentPCA,1, CV_32FC1, cluster1_1_5_pcaEigenValLoad.data());
    pcaANN1_1_5.eigenvectors = cv::Mat(numComponentPCA,numVarANN, CV_32FC1, cluster1_1_5_pcaEigenVecLoad.data());
    
    std::string cluster1_1_6_pcaMeanLoadPATH = casePath + cluster1_1_6_Folder + "PCA_mean.txt";
    std::string cluster1_1_6_pcaEigenValLoadPATH = casePath + cluster1_1_6_Folder + "PCA_eigenValues.txt";
    std::string cluster1_1_6_pcaEigenVecLoadPATH = casePath + cluster1_1_6_Folder + "PCA_eigenVectors.txt";
    vector<float> cluster1_1_6_pcaMeanLoad;
    vector<float> cluster1_1_6_pcaEigenValLoad;
    vector<float> cluster1_1_6_pcaEigenVecLoad;
    readFromCommaDelimitedFile_Float(cluster1_1_6_pcaMeanLoadPATH, cluster1_1_6_pcaMeanLoad);
    readFromCommaDelimitedFile_Float(cluster1_1_6_pcaEigenValLoadPATH, cluster1_1_6_pcaEigenValLoad);
    readFromCommaDelimitedFile_Float(cluster1_1_6_pcaEigenVecLoadPATH, cluster1_1_6_pcaEigenVecLoad);
    cv::PCA pcaANN1_1_6;
    pcaANN1_1_6.mean = cv::Mat(1, numVarANN, CV_32FC1, cluster1_1_6_pcaMeanLoad.data());
    pcaANN1_1_6.eigenvalues = cv::Mat(numComponentPCA,1, CV_32FC1, cluster1_1_6_pcaEigenValLoad.data());
    pcaANN1_1_6.eigenvectors = cv::Mat(numComponentPCA,numVarANN, CV_32FC1, cluster1_1_6_pcaEigenVecLoad.data());
    
    std::string cluster1_1_7_pcaMeanLoadPATH = casePath + cluster1_1_7_Folder + "PCA_mean.txt";
    std::string cluster1_1_7_pcaEigenValLoadPATH = casePath + cluster1_1_7_Folder + "PCA_eigenValues.txt";
    std::string cluster1_1_7_pcaEigenVecLoadPATH = casePath + cluster1_1_7_Folder + "PCA_eigenVectors.txt";
    vector<float> cluster1_1_7_pcaMeanLoad;
    vector<float> cluster1_1_7_pcaEigenValLoad;
    vector<float> cluster1_1_7_pcaEigenVecLoad;
    readFromCommaDelimitedFile_Float(cluster1_1_7_pcaMeanLoadPATH, cluster1_1_7_pcaMeanLoad);
    readFromCommaDelimitedFile_Float(cluster1_1_7_pcaEigenValLoadPATH, cluster1_1_7_pcaEigenValLoad);
    readFromCommaDelimitedFile_Float(cluster1_1_7_pcaEigenVecLoadPATH, cluster1_1_7_pcaEigenVecLoad);
    cv::PCA pcaANN1_1_7;
    pcaANN1_1_7.mean = cv::Mat(1, numVarANN, CV_32FC1, cluster1_1_7_pcaMeanLoad.data());
    pcaANN1_1_7.eigenvalues = cv::Mat(numComponentPCA,1, CV_32FC1, cluster1_1_7_pcaEigenValLoad.data());
    pcaANN1_1_7.eigenvectors = cv::Mat(numComponentPCA,numVarANN, CV_32FC1, cluster1_1_7_pcaEigenVecLoad.data());
    
    std::string cluster1_1_8_pcaMeanLoadPATH = casePath + cluster1_1_8_Folder + "PCA_mean.txt";
    std::string cluster1_1_8_pcaEigenValLoadPATH = casePath + cluster1_1_8_Folder + "PCA_eigenValues.txt";
    std::string cluster1_1_8_pcaEigenVecLoadPATH = casePath + cluster1_1_8_Folder + "PCA_eigenVectors.txt";
    vector<float> cluster1_1_8_pcaMeanLoad;
    vector<float> cluster1_1_8_pcaEigenValLoad;
    vector<float> cluster1_1_8_pcaEigenVecLoad;
    readFromCommaDelimitedFile_Float(cluster1_1_8_pcaMeanLoadPATH, cluster1_1_8_pcaMeanLoad);
    readFromCommaDelimitedFile_Float(cluster1_1_8_pcaEigenValLoadPATH, cluster1_1_8_pcaEigenValLoad);
    readFromCommaDelimitedFile_Float(cluster1_1_8_pcaEigenVecLoadPATH, cluster1_1_8_pcaEigenVecLoad);
    cv::PCA pcaANN1_1_8;
    pcaANN1_1_8.mean = cv::Mat(1, numVarANN, CV_32FC1, cluster1_1_8_pcaMeanLoad.data());
    pcaANN1_1_8.eigenvalues = cv::Mat(numComponentPCA,1, CV_32FC1, cluster1_1_8_pcaEigenValLoad.data());
    pcaANN1_1_8.eigenvectors = cv::Mat(numComponentPCA,numVarANN, CV_32FC1, cluster1_1_8_pcaEigenVecLoad.data());
    
    std::string cluster1_1_9_pcaMeanLoadPATH = casePath + cluster1_1_9_Folder + "PCA_mean.txt";
    std::string cluster1_1_9_pcaEigenValLoadPATH = casePath + cluster1_1_9_Folder + "PCA_eigenValues.txt";
    std::string cluster1_1_9_pcaEigenVecLoadPATH = casePath + cluster1_1_9_Folder + "PCA_eigenVectors.txt";
    vector<float> cluster1_1_9_pcaMeanLoad;
    vector<float> cluster1_1_9_pcaEigenValLoad;
    vector<float> cluster1_1_9_pcaEigenVecLoad;
    readFromCommaDelimitedFile_Float(cluster1_1_9_pcaMeanLoadPATH, cluster1_1_9_pcaMeanLoad);
    readFromCommaDelimitedFile_Float(cluster1_1_9_pcaEigenValLoadPATH, cluster1_1_9_pcaEigenValLoad);
    readFromCommaDelimitedFile_Float(cluster1_1_9_pcaEigenVecLoadPATH, cluster1_1_9_pcaEigenVecLoad);
    cv::PCA pcaANN1_1_9;
    pcaANN1_1_9.mean = cv::Mat(1, numVarANN, CV_32FC1, cluster1_1_9_pcaMeanLoad.data());
    pcaANN1_1_9.eigenvalues = cv::Mat(numComponentPCA,1, CV_32FC1, cluster1_1_9_pcaEigenValLoad.data());
    pcaANN1_1_9.eigenvectors = cv::Mat(numComponentPCA,numVarANN, CV_32FC1, cluster1_1_9_pcaEigenVecLoad.data());
    
    std::string cluster1_1_10_pcaMeanLoadPATH = casePath + cluster1_1_10_Folder + "PCA_mean.txt";
    std::string cluster1_1_10_pcaEigenValLoadPATH = casePath + cluster1_1_10_Folder + "PCA_eigenValues.txt";
    std::string cluster1_1_10_pcaEigenVecLoadPATH = casePath + cluster1_1_10_Folder + "PCA_eigenVectors.txt";
    vector<float> cluster1_1_10_pcaMeanLoad;
    vector<float> cluster1_1_10_pcaEigenValLoad;
    vector<float> cluster1_1_10_pcaEigenVecLoad;
    readFromCommaDelimitedFile_Float(cluster1_1_10_pcaMeanLoadPATH, cluster1_1_10_pcaMeanLoad);
    readFromCommaDelimitedFile_Float(cluster1_1_10_pcaEigenValLoadPATH, cluster1_1_10_pcaEigenValLoad);
    readFromCommaDelimitedFile_Float(cluster1_1_10_pcaEigenVecLoadPATH, cluster1_1_10_pcaEigenVecLoad);
    cv::PCA pcaANN1_1_10;
    pcaANN1_1_10.mean = cv::Mat(1, numVarANN, CV_32FC1, cluster1_1_10_pcaMeanLoad.data());
    pcaANN1_1_10.eigenvalues = cv::Mat(numComponentPCA,1, CV_32FC1, cluster1_1_10_pcaEigenValLoad.data());
    pcaANN1_1_10.eigenvectors = cv::Mat(numComponentPCA,numVarANN, CV_32FC1, cluster1_1_10_pcaEigenVecLoad.data());
    
    std::string cluster1_1_11_pcaMeanLoadPATH = casePath + cluster1_1_11_Folder + "PCA_mean.txt";
    std::string cluster1_1_11_pcaEigenValLoadPATH = casePath + cluster1_1_11_Folder + "PCA_eigenValues.txt";
    std::string cluster1_1_11_pcaEigenVecLoadPATH = casePath + cluster1_1_11_Folder + "PCA_eigenVectors.txt";
    vector<float> cluster1_1_11_pcaMeanLoad;
    vector<float> cluster1_1_11_pcaEigenValLoad;
    vector<float> cluster1_1_11_pcaEigenVecLoad;
    readFromCommaDelimitedFile_Float(cluster1_1_11_pcaMeanLoadPATH, cluster1_1_11_pcaMeanLoad);
    readFromCommaDelimitedFile_Float(cluster1_1_11_pcaEigenValLoadPATH, cluster1_1_11_pcaEigenValLoad);
    readFromCommaDelimitedFile_Float(cluster1_1_11_pcaEigenVecLoadPATH, cluster1_1_11_pcaEigenVecLoad);
    cv::PCA pcaANN1_1_11;
    pcaANN1_1_11.mean = cv::Mat(1, numVarANN, CV_32FC1, cluster1_1_11_pcaMeanLoad.data());
    pcaANN1_1_11.eigenvalues = cv::Mat(numComponentPCA,1, CV_32FC1, cluster1_1_11_pcaEigenValLoad.data());
    pcaANN1_1_11.eigenvectors = cv::Mat(numComponentPCA,numVarANN, CV_32FC1, cluster1_1_11_pcaEigenVecLoad.data());
    
    std::string cluster1_1_12_pcaMeanLoadPATH = casePath + cluster1_1_12_Folder + "PCA_mean.txt";
    std::string cluster1_1_12_pcaEigenValLoadPATH = casePath + cluster1_1_12_Folder + "PCA_eigenValues.txt";
    std::string cluster1_1_12_pcaEigenVecLoadPATH = casePath + cluster1_1_12_Folder + "PCA_eigenVectors.txt";
    vector<float> cluster1_1_12_pcaMeanLoad;
    vector<float> cluster1_1_12_pcaEigenValLoad;
    vector<float> cluster1_1_12_pcaEigenVecLoad;
    readFromCommaDelimitedFile_Float(cluster1_1_12_pcaMeanLoadPATH, cluster1_1_12_pcaMeanLoad);
    readFromCommaDelimitedFile_Float(cluster1_1_12_pcaEigenValLoadPATH, cluster1_1_12_pcaEigenValLoad);
    readFromCommaDelimitedFile_Float(cluster1_1_12_pcaEigenVecLoadPATH, cluster1_1_12_pcaEigenVecLoad);
    cv::PCA pcaANN1_1_12;
    pcaANN1_1_12.mean = cv::Mat(1, numVarANN, CV_32FC1, cluster1_1_12_pcaMeanLoad.data());
    pcaANN1_1_12.eigenvalues = cv::Mat(numComponentPCA,1, CV_32FC1, cluster1_1_12_pcaEigenValLoad.data());
    pcaANN1_1_12.eigenvectors = cv::Mat(numComponentPCA,numVarANN, CV_32FC1, cluster1_1_12_pcaEigenVecLoad.data());
    
    std::string cluster1_1_13_pcaMeanLoadPATH = casePath + cluster1_1_13_Folder + "PCA_mean.txt";
    std::string cluster1_1_13_pcaEigenValLoadPATH = casePath + cluster1_1_13_Folder + "PCA_eigenValues.txt";
    std::string cluster1_1_13_pcaEigenVecLoadPATH = casePath + cluster1_1_13_Folder + "PCA_eigenVectors.txt";
    vector<float> cluster1_1_13_pcaMeanLoad;
    vector<float> cluster1_1_13_pcaEigenValLoad;
    vector<float> cluster1_1_13_pcaEigenVecLoad;
    readFromCommaDelimitedFile_Float(cluster1_1_13_pcaMeanLoadPATH, cluster1_1_13_pcaMeanLoad);
    readFromCommaDelimitedFile_Float(cluster1_1_13_pcaEigenValLoadPATH, cluster1_1_13_pcaEigenValLoad);
    readFromCommaDelimitedFile_Float(cluster1_1_13_pcaEigenVecLoadPATH, cluster1_1_13_pcaEigenVecLoad);
    cv::PCA pcaANN1_1_13;
    pcaANN1_1_13.mean = cv::Mat(1, numVarANN, CV_32FC1, cluster1_1_13_pcaMeanLoad.data());
    pcaANN1_1_13.eigenvalues = cv::Mat(numComponentPCA,1, CV_32FC1, cluster1_1_13_pcaEigenValLoad.data());
    pcaANN1_1_13.eigenvectors = cv::Mat(numComponentPCA,numVarANN, CV_32FC1, cluster1_1_13_pcaEigenVecLoad.data());
    
    std::string cluster1_1_14_pcaMeanLoadPATH = casePath + cluster1_1_14_Folder + "PCA_mean.txt";
    std::string cluster1_1_14_pcaEigenValLoadPATH = casePath + cluster1_1_14_Folder + "PCA_eigenValues.txt";
    std::string cluster1_1_14_pcaEigenVecLoadPATH = casePath + cluster1_1_14_Folder + "PCA_eigenVectors.txt";
    vector<float> cluster1_1_14_pcaMeanLoad;
    vector<float> cluster1_1_14_pcaEigenValLoad;
    vector<float> cluster1_1_14_pcaEigenVecLoad;
    readFromCommaDelimitedFile_Float(cluster1_1_14_pcaMeanLoadPATH, cluster1_1_14_pcaMeanLoad);
    readFromCommaDelimitedFile_Float(cluster1_1_14_pcaEigenValLoadPATH, cluster1_1_14_pcaEigenValLoad);
    readFromCommaDelimitedFile_Float(cluster1_1_14_pcaEigenVecLoadPATH, cluster1_1_14_pcaEigenVecLoad);
    cv::PCA pcaANN1_1_14;
    pcaANN1_1_14.mean = cv::Mat(1, numVarANN, CV_32FC1, cluster1_1_14_pcaMeanLoad.data());
    pcaANN1_1_14.eigenvalues = cv::Mat(numComponentPCA,1, CV_32FC1, cluster1_1_14_pcaEigenValLoad.data());
    pcaANN1_1_14.eigenvectors = cv::Mat(numComponentPCA,numVarANN, CV_32FC1, cluster1_1_14_pcaEigenVecLoad.data());
    
    
    // GrandChild of parentCluster1 - childCluster2
    std::string cluster1_2_0_pcaMeanLoadPATH = casePath + cluster1_2_0_Folder + "PCA_mean.txt";
    std::string cluster1_2_0_pcaEigenValLoadPATH = casePath + cluster1_2_0_Folder + "PCA_eigenValues.txt";
    std::string cluster1_2_0_pcaEigenVecLoadPATH = casePath + cluster1_2_0_Folder + "PCA_eigenVectors.txt";
    vector<float> cluster1_2_0_pcaMeanLoad;
    vector<float> cluster1_2_0_pcaEigenValLoad;
    vector<float> cluster1_2_0_pcaEigenVecLoad;
    readFromCommaDelimitedFile_Float(cluster1_2_0_pcaMeanLoadPATH, cluster1_2_0_pcaMeanLoad);
    readFromCommaDelimitedFile_Float(cluster1_2_0_pcaEigenValLoadPATH, cluster1_2_0_pcaEigenValLoad);
    readFromCommaDelimitedFile_Float(cluster1_2_0_pcaEigenVecLoadPATH, cluster1_2_0_pcaEigenVecLoad);
    cv::PCA pcaANN1_2_0;
    pcaANN1_2_0.mean = cv::Mat(1, numVarANN, CV_32FC1, cluster1_2_0_pcaMeanLoad.data());
    pcaANN1_2_0.eigenvalues = cv::Mat(numComponentPCA,1, CV_32FC1, cluster1_2_0_pcaEigenValLoad.data());
    pcaANN1_2_0.eigenvectors = cv::Mat(numComponentPCA,numVarANN, CV_32FC1, cluster1_2_0_pcaEigenVecLoad.data());
    
    std::string cluster1_2_1_pcaMeanLoadPATH = casePath + cluster1_2_1_Folder + "PCA_mean.txt";
    std::string cluster1_2_1_pcaEigenValLoadPATH = casePath + cluster1_2_1_Folder + "PCA_eigenValues.txt";
    std::string cluster1_2_1_pcaEigenVecLoadPATH = casePath + cluster1_2_1_Folder + "PCA_eigenVectors.txt";
    vector<float> cluster1_2_1_pcaMeanLoad;
    vector<float> cluster1_2_1_pcaEigenValLoad;
    vector<float> cluster1_2_1_pcaEigenVecLoad;
    readFromCommaDelimitedFile_Float(cluster1_2_1_pcaMeanLoadPATH, cluster1_2_1_pcaMeanLoad);
    readFromCommaDelimitedFile_Float(cluster1_2_1_pcaEigenValLoadPATH, cluster1_2_1_pcaEigenValLoad);
    readFromCommaDelimitedFile_Float(cluster1_2_1_pcaEigenVecLoadPATH, cluster1_2_1_pcaEigenVecLoad);
    cv::PCA pcaANN1_2_1;
    pcaANN1_2_1.mean = cv::Mat(1, numVarANN, CV_32FC1, cluster1_2_1_pcaMeanLoad.data());
    pcaANN1_2_1.eigenvalues = cv::Mat(numComponentPCA,1, CV_32FC1, cluster1_2_1_pcaEigenValLoad.data());
    pcaANN1_2_1.eigenvectors = cv::Mat(numComponentPCA,numVarANN, CV_32FC1, cluster1_2_1_pcaEigenVecLoad.data());
    
    std::string cluster1_2_2_pcaMeanLoadPATH = casePath + cluster1_2_2_Folder + "PCA_mean.txt";
    std::string cluster1_2_2_pcaEigenValLoadPATH = casePath + cluster1_2_2_Folder + "PCA_eigenValues.txt";
    std::string cluster1_2_2_pcaEigenVecLoadPATH = casePath + cluster1_2_2_Folder + "PCA_eigenVectors.txt";
    vector<float> cluster1_2_2_pcaMeanLoad;
    vector<float> cluster1_2_2_pcaEigenValLoad;
    vector<float> cluster1_2_2_pcaEigenVecLoad;
    readFromCommaDelimitedFile_Float(cluster1_2_2_pcaMeanLoadPATH, cluster1_2_2_pcaMeanLoad);
    readFromCommaDelimitedFile_Float(cluster1_2_2_pcaEigenValLoadPATH, cluster1_2_2_pcaEigenValLoad);
    readFromCommaDelimitedFile_Float(cluster1_2_2_pcaEigenVecLoadPATH, cluster1_2_2_pcaEigenVecLoad);
    cv::PCA pcaANN1_2_2;
    pcaANN1_2_2.mean = cv::Mat(1, numVarANN, CV_32FC1, cluster1_2_2_pcaMeanLoad.data());
    pcaANN1_2_2.eigenvalues = cv::Mat(numComponentPCA,1, CV_32FC1, cluster1_2_2_pcaEigenValLoad.data());
    pcaANN1_2_2.eigenvectors = cv::Mat(numComponentPCA,numVarANN, CV_32FC1, cluster1_2_2_pcaEigenVecLoad.data());
    
    std::string cluster1_2_3_pcaMeanLoadPATH = casePath + cluster1_2_3_Folder + "PCA_mean.txt";
    std::string cluster1_2_3_pcaEigenValLoadPATH = casePath + cluster1_2_3_Folder + "PCA_eigenValues.txt";
    std::string cluster1_2_3_pcaEigenVecLoadPATH = casePath + cluster1_2_3_Folder + "PCA_eigenVectors.txt";
    vector<float> cluster1_2_3_pcaMeanLoad;
    vector<float> cluster1_2_3_pcaEigenValLoad;
    vector<float> cluster1_2_3_pcaEigenVecLoad;
    readFromCommaDelimitedFile_Float(cluster1_2_3_pcaMeanLoadPATH, cluster1_2_3_pcaMeanLoad);
    readFromCommaDelimitedFile_Float(cluster1_2_3_pcaEigenValLoadPATH, cluster1_2_3_pcaEigenValLoad);
    readFromCommaDelimitedFile_Float(cluster1_2_3_pcaEigenVecLoadPATH, cluster1_2_3_pcaEigenVecLoad);
    cv::PCA pcaANN1_2_3;
    pcaANN1_2_3.mean = cv::Mat(1, numVarANN, CV_32FC1, cluster1_2_3_pcaMeanLoad.data());
    pcaANN1_2_3.eigenvalues = cv::Mat(numComponentPCA,1, CV_32FC1, cluster1_2_3_pcaEigenValLoad.data());
    pcaANN1_2_3.eigenvectors = cv::Mat(numComponentPCA,numVarANN, CV_32FC1, cluster1_2_3_pcaEigenVecLoad.data());
    
    std::string cluster1_2_4_pcaMeanLoadPATH = casePath + cluster1_2_4_Folder + "PCA_mean.txt";
    std::string cluster1_2_4_pcaEigenValLoadPATH = casePath + cluster1_2_4_Folder + "PCA_eigenValues.txt";
    std::string cluster1_2_4_pcaEigenVecLoadPATH = casePath + cluster1_2_4_Folder + "PCA_eigenVectors.txt";
    vector<float> cluster1_2_4_pcaMeanLoad;
    vector<float> cluster1_2_4_pcaEigenValLoad;
    vector<float> cluster1_2_4_pcaEigenVecLoad;
    readFromCommaDelimitedFile_Float(cluster1_2_4_pcaMeanLoadPATH, cluster1_2_4_pcaMeanLoad);
    readFromCommaDelimitedFile_Float(cluster1_2_4_pcaEigenValLoadPATH, cluster1_2_4_pcaEigenValLoad);
    readFromCommaDelimitedFile_Float(cluster1_2_4_pcaEigenVecLoadPATH, cluster1_2_4_pcaEigenVecLoad);
    cv::PCA pcaANN1_2_4;
    pcaANN1_2_4.mean = cv::Mat(1, numVarANN, CV_32FC1, cluster1_2_4_pcaMeanLoad.data());
    pcaANN1_2_4.eigenvalues = cv::Mat(numComponentPCA,1, CV_32FC1, cluster1_2_4_pcaEigenValLoad.data());
    pcaANN1_2_4.eigenvectors = cv::Mat(numComponentPCA,numVarANN, CV_32FC1, cluster1_2_4_pcaEigenVecLoad.data());
    
    std::string cluster1_2_5_pcaMeanLoadPATH = casePath + cluster1_2_5_Folder + "PCA_mean.txt";
    std::string cluster1_2_5_pcaEigenValLoadPATH = casePath + cluster1_2_5_Folder + "PCA_eigenValues.txt";
    std::string cluster1_2_5_pcaEigenVecLoadPATH = casePath + cluster1_2_5_Folder + "PCA_eigenVectors.txt";
    vector<float> cluster1_2_5_pcaMeanLoad;
    vector<float> cluster1_2_5_pcaEigenValLoad;
    vector<float> cluster1_2_5_pcaEigenVecLoad;
    readFromCommaDelimitedFile_Float(cluster1_2_5_pcaMeanLoadPATH, cluster1_2_5_pcaMeanLoad);
    readFromCommaDelimitedFile_Float(cluster1_2_5_pcaEigenValLoadPATH, cluster1_2_5_pcaEigenValLoad);
    readFromCommaDelimitedFile_Float(cluster1_2_5_pcaEigenVecLoadPATH, cluster1_2_5_pcaEigenVecLoad);
    cv::PCA pcaANN1_2_5;
    pcaANN1_2_5.mean = cv::Mat(1, numVarANN, CV_32FC1, cluster1_2_5_pcaMeanLoad.data());
    pcaANN1_2_5.eigenvalues = cv::Mat(numComponentPCA,1, CV_32FC1, cluster1_2_5_pcaEigenValLoad.data());
    pcaANN1_2_5.eigenvectors = cv::Mat(numComponentPCA,numVarANN, CV_32FC1, cluster1_2_5_pcaEigenVecLoad.data());
    
    std::string cluster1_2_6_pcaMeanLoadPATH = casePath + cluster1_2_6_Folder + "PCA_mean.txt";
    std::string cluster1_2_6_pcaEigenValLoadPATH = casePath + cluster1_2_6_Folder + "PCA_eigenValues.txt";
    std::string cluster1_2_6_pcaEigenVecLoadPATH = casePath + cluster1_2_6_Folder + "PCA_eigenVectors.txt";
    vector<float> cluster1_2_6_pcaMeanLoad;
    vector<float> cluster1_2_6_pcaEigenValLoad;
    vector<float> cluster1_2_6_pcaEigenVecLoad;
    readFromCommaDelimitedFile_Float(cluster1_2_6_pcaMeanLoadPATH, cluster1_2_6_pcaMeanLoad);
    readFromCommaDelimitedFile_Float(cluster1_2_6_pcaEigenValLoadPATH, cluster1_2_6_pcaEigenValLoad);
    readFromCommaDelimitedFile_Float(cluster1_2_6_pcaEigenVecLoadPATH, cluster1_2_6_pcaEigenVecLoad);
    cv::PCA pcaANN1_2_6;
    pcaANN1_2_6.mean = cv::Mat(1, numVarANN, CV_32FC1, cluster1_2_6_pcaMeanLoad.data());
    pcaANN1_2_6.eigenvalues = cv::Mat(numComponentPCA,1, CV_32FC1, cluster1_2_6_pcaEigenValLoad.data());
    pcaANN1_2_6.eigenvectors = cv::Mat(numComponentPCA,numVarANN, CV_32FC1, cluster1_2_6_pcaEigenVecLoad.data());
    
    std::string cluster1_2_7_pcaMeanLoadPATH = casePath + cluster1_2_7_Folder + "PCA_mean.txt";
    std::string cluster1_2_7_pcaEigenValLoadPATH = casePath + cluster1_2_7_Folder + "PCA_eigenValues.txt";
    std::string cluster1_2_7_pcaEigenVecLoadPATH = casePath + cluster1_2_7_Folder + "PCA_eigenVectors.txt";
    vector<float> cluster1_2_7_pcaMeanLoad;
    vector<float> cluster1_2_7_pcaEigenValLoad;
    vector<float> cluster1_2_7_pcaEigenVecLoad;
    readFromCommaDelimitedFile_Float(cluster1_2_7_pcaMeanLoadPATH, cluster1_2_7_pcaMeanLoad);
    readFromCommaDelimitedFile_Float(cluster1_2_7_pcaEigenValLoadPATH, cluster1_2_7_pcaEigenValLoad);
    readFromCommaDelimitedFile_Float(cluster1_2_7_pcaEigenVecLoadPATH, cluster1_2_7_pcaEigenVecLoad);
    cv::PCA pcaANN1_2_7;
    pcaANN1_2_7.mean = cv::Mat(1, numVarANN, CV_32FC1, cluster1_2_7_pcaMeanLoad.data());
    pcaANN1_2_7.eigenvalues = cv::Mat(numComponentPCA,1, CV_32FC1, cluster1_2_7_pcaEigenValLoad.data());
    pcaANN1_2_7.eigenvectors = cv::Mat(numComponentPCA,numVarANN, CV_32FC1, cluster1_2_7_pcaEigenVecLoad.data());
    
    std::string cluster1_2_8_pcaMeanLoadPATH = casePath + cluster1_2_8_Folder + "PCA_mean.txt";
    std::string cluster1_2_8_pcaEigenValLoadPATH = casePath + cluster1_2_8_Folder + "PCA_eigenValues.txt";
    std::string cluster1_2_8_pcaEigenVecLoadPATH = casePath + cluster1_2_8_Folder + "PCA_eigenVectors.txt";
    vector<float> cluster1_2_8_pcaMeanLoad;
    vector<float> cluster1_2_8_pcaEigenValLoad;
    vector<float> cluster1_2_8_pcaEigenVecLoad;
    readFromCommaDelimitedFile_Float(cluster1_2_8_pcaMeanLoadPATH, cluster1_2_8_pcaMeanLoad);
    readFromCommaDelimitedFile_Float(cluster1_2_8_pcaEigenValLoadPATH, cluster1_2_8_pcaEigenValLoad);
    readFromCommaDelimitedFile_Float(cluster1_2_8_pcaEigenVecLoadPATH, cluster1_2_8_pcaEigenVecLoad);
    cv::PCA pcaANN1_2_8;
    pcaANN1_2_8.mean = cv::Mat(1, numVarANN, CV_32FC1, cluster1_2_8_pcaMeanLoad.data());
    pcaANN1_2_8.eigenvalues = cv::Mat(numComponentPCA,1, CV_32FC1, cluster1_2_8_pcaEigenValLoad.data());
    pcaANN1_2_8.eigenvectors = cv::Mat(numComponentPCA,numVarANN, CV_32FC1, cluster1_2_8_pcaEigenVecLoad.data());
    
    std::string cluster1_2_9_pcaMeanLoadPATH = casePath + cluster1_2_9_Folder + "PCA_mean.txt";
    std::string cluster1_2_9_pcaEigenValLoadPATH = casePath + cluster1_2_9_Folder + "PCA_eigenValues.txt";
    std::string cluster1_2_9_pcaEigenVecLoadPATH = casePath + cluster1_2_9_Folder + "PCA_eigenVectors.txt";
    vector<float> cluster1_2_9_pcaMeanLoad;
    vector<float> cluster1_2_9_pcaEigenValLoad;
    vector<float> cluster1_2_9_pcaEigenVecLoad;
    readFromCommaDelimitedFile_Float(cluster1_2_9_pcaMeanLoadPATH, cluster1_2_9_pcaMeanLoad);
    readFromCommaDelimitedFile_Float(cluster1_2_9_pcaEigenValLoadPATH, cluster1_2_9_pcaEigenValLoad);
    readFromCommaDelimitedFile_Float(cluster1_2_9_pcaEigenVecLoadPATH, cluster1_2_9_pcaEigenVecLoad);
    cv::PCA pcaANN1_2_9;
    pcaANN1_2_9.mean = cv::Mat(1, numVarANN, CV_32FC1, cluster1_2_9_pcaMeanLoad.data());
    pcaANN1_2_9.eigenvalues = cv::Mat(numComponentPCA,1, CV_32FC1, cluster1_2_9_pcaEigenValLoad.data());
    pcaANN1_2_9.eigenvectors = cv::Mat(numComponentPCA,numVarANN, CV_32FC1, cluster1_2_9_pcaEigenVecLoad.data());
    
    std::string cluster1_2_10_pcaMeanLoadPATH = casePath + cluster1_2_10_Folder + "PCA_mean.txt";
    std::string cluster1_2_10_pcaEigenValLoadPATH = casePath + cluster1_2_10_Folder + "PCA_eigenValues.txt";
    std::string cluster1_2_10_pcaEigenVecLoadPATH = casePath + cluster1_2_10_Folder + "PCA_eigenVectors.txt";
    vector<float> cluster1_2_10_pcaMeanLoad;
    vector<float> cluster1_2_10_pcaEigenValLoad;
    vector<float> cluster1_2_10_pcaEigenVecLoad;
    readFromCommaDelimitedFile_Float(cluster1_2_10_pcaMeanLoadPATH, cluster1_2_10_pcaMeanLoad);
    readFromCommaDelimitedFile_Float(cluster1_2_10_pcaEigenValLoadPATH, cluster1_2_10_pcaEigenValLoad);
    readFromCommaDelimitedFile_Float(cluster1_2_10_pcaEigenVecLoadPATH, cluster1_2_10_pcaEigenVecLoad);
    cv::PCA pcaANN1_2_10;
    pcaANN1_2_10.mean = cv::Mat(1, numVarANN, CV_32FC1, cluster1_2_10_pcaMeanLoad.data());
    pcaANN1_2_10.eigenvalues = cv::Mat(numComponentPCA,1, CV_32FC1, cluster1_2_10_pcaEigenValLoad.data());
    pcaANN1_2_10.eigenvectors = cv::Mat(numComponentPCA,numVarANN, CV_32FC1, cluster1_2_10_pcaEigenVecLoad.data());
    
    std::string cluster1_2_11_pcaMeanLoadPATH = casePath + cluster1_2_11_Folder + "PCA_mean.txt";
    std::string cluster1_2_11_pcaEigenValLoadPATH = casePath + cluster1_2_11_Folder + "PCA_eigenValues.txt";
    std::string cluster1_2_11_pcaEigenVecLoadPATH = casePath + cluster1_2_11_Folder + "PCA_eigenVectors.txt";
    vector<float> cluster1_2_11_pcaMeanLoad;
    vector<float> cluster1_2_11_pcaEigenValLoad;
    vector<float> cluster1_2_11_pcaEigenVecLoad;
    readFromCommaDelimitedFile_Float(cluster1_2_11_pcaMeanLoadPATH, cluster1_2_11_pcaMeanLoad);
    readFromCommaDelimitedFile_Float(cluster1_2_11_pcaEigenValLoadPATH, cluster1_2_11_pcaEigenValLoad);
    readFromCommaDelimitedFile_Float(cluster1_2_11_pcaEigenVecLoadPATH, cluster1_2_11_pcaEigenVecLoad);
    cv::PCA pcaANN1_2_11;
    pcaANN1_2_11.mean = cv::Mat(1, numVarANN, CV_32FC1, cluster1_2_11_pcaMeanLoad.data());
    pcaANN1_2_11.eigenvalues = cv::Mat(numComponentPCA,1, CV_32FC1, cluster1_2_11_pcaEigenValLoad.data());
    pcaANN1_2_11.eigenvectors = cv::Mat(numComponentPCA,numVarANN, CV_32FC1, cluster1_2_11_pcaEigenVecLoad.data());
    
    std::string cluster1_2_12_pcaMeanLoadPATH = casePath + cluster1_2_12_Folder + "PCA_mean.txt";
    std::string cluster1_2_12_pcaEigenValLoadPATH = casePath + cluster1_2_12_Folder + "PCA_eigenValues.txt";
    std::string cluster1_2_12_pcaEigenVecLoadPATH = casePath + cluster1_2_12_Folder + "PCA_eigenVectors.txt";
    vector<float> cluster1_2_12_pcaMeanLoad;
    vector<float> cluster1_2_12_pcaEigenValLoad;
    vector<float> cluster1_2_12_pcaEigenVecLoad;
    readFromCommaDelimitedFile_Float(cluster1_2_12_pcaMeanLoadPATH, cluster1_2_12_pcaMeanLoad);
    readFromCommaDelimitedFile_Float(cluster1_2_12_pcaEigenValLoadPATH, cluster1_2_12_pcaEigenValLoad);
    readFromCommaDelimitedFile_Float(cluster1_2_12_pcaEigenVecLoadPATH, cluster1_2_12_pcaEigenVecLoad);
    cv::PCA pcaANN1_2_12;
    pcaANN1_2_12.mean = cv::Mat(1, numVarANN, CV_32FC1, cluster1_2_12_pcaMeanLoad.data());
    pcaANN1_2_12.eigenvalues = cv::Mat(numComponentPCA,1, CV_32FC1, cluster1_2_12_pcaEigenValLoad.data());
    pcaANN1_2_12.eigenvectors = cv::Mat(numComponentPCA,numVarANN, CV_32FC1, cluster1_2_12_pcaEigenVecLoad.data());
    
    std::string cluster1_2_13_pcaMeanLoadPATH = casePath + cluster1_2_13_Folder + "PCA_mean.txt";
    std::string cluster1_2_13_pcaEigenValLoadPATH = casePath + cluster1_2_13_Folder + "PCA_eigenValues.txt";
    std::string cluster1_2_13_pcaEigenVecLoadPATH = casePath + cluster1_2_13_Folder + "PCA_eigenVectors.txt";
    vector<float> cluster1_2_13_pcaMeanLoad;
    vector<float> cluster1_2_13_pcaEigenValLoad;
    vector<float> cluster1_2_13_pcaEigenVecLoad;
    readFromCommaDelimitedFile_Float(cluster1_2_13_pcaMeanLoadPATH, cluster1_2_13_pcaMeanLoad);
    readFromCommaDelimitedFile_Float(cluster1_2_13_pcaEigenValLoadPATH, cluster1_2_13_pcaEigenValLoad);
    readFromCommaDelimitedFile_Float(cluster1_2_13_pcaEigenVecLoadPATH, cluster1_2_13_pcaEigenVecLoad);
    cv::PCA pcaANN1_2_13;
    pcaANN1_2_13.mean = cv::Mat(1, numVarANN, CV_32FC1, cluster1_2_13_pcaMeanLoad.data());
    pcaANN1_2_13.eigenvalues = cv::Mat(numComponentPCA,1, CV_32FC1, cluster1_2_13_pcaEigenValLoad.data());
    pcaANN1_2_13.eigenvectors = cv::Mat(numComponentPCA,numVarANN, CV_32FC1, cluster1_2_13_pcaEigenVecLoad.data());
    
    std::string cluster1_2_14_pcaMeanLoadPATH = casePath + cluster1_2_14_Folder + "PCA_mean.txt";
    std::string cluster1_2_14_pcaEigenValLoadPATH = casePath + cluster1_2_14_Folder + "PCA_eigenValues.txt";
    std::string cluster1_2_14_pcaEigenVecLoadPATH = casePath + cluster1_2_14_Folder + "PCA_eigenVectors.txt";
    vector<float> cluster1_2_14_pcaMeanLoad;
    vector<float> cluster1_2_14_pcaEigenValLoad;
    vector<float> cluster1_2_14_pcaEigenVecLoad;
    readFromCommaDelimitedFile_Float(cluster1_2_14_pcaMeanLoadPATH, cluster1_2_14_pcaMeanLoad);
    readFromCommaDelimitedFile_Float(cluster1_2_14_pcaEigenValLoadPATH, cluster1_2_14_pcaEigenValLoad);
    readFromCommaDelimitedFile_Float(cluster1_2_14_pcaEigenVecLoadPATH, cluster1_2_14_pcaEigenVecLoad);
    cv::PCA pcaANN1_2_14;
    pcaANN1_2_14.mean = cv::Mat(1, numVarANN, CV_32FC1, cluster1_2_14_pcaMeanLoad.data());
    pcaANN1_2_14.eigenvalues = cv::Mat(numComponentPCA,1, CV_32FC1, cluster1_2_14_pcaEigenValLoad.data());
    pcaANN1_2_14.eigenvectors = cv::Mat(numComponentPCA,numVarANN, CV_32FC1, cluster1_2_14_pcaEigenVecLoad.data());
 */ //  HuuTri@20211006: Commented - ANN without PCA 

    // Declare ANN regression for each cluster
    // Child of parentCluster0
    std::string frozenPathReg_cluster0_0 = casePath + cluster0_0_Folder + "frozen_bestModel.pb";
    Model modelANNReg0_0(frozenPathReg_cluster0_0);
    std::vector<std::string> operNameANN0_0;
    operNameANN0_0 = modelANNReg0_0.get_operations();
    Tensor inputPCAANNReg0_0{modelANNReg0_0, operNameANN0_0[0]}; //First tensor = Input (Standardized values)
    Tensor outputStandardizeANNReg0_0{modelANNReg0_0, operNameANN0_0[operNameANN0_0.size()-1]};    //Last tensor = Output (Standardized values)
    
    std::string frozenPathReg_cluster0_1 = casePath + cluster0_1_Folder + "frozen_bestModel.pb";
    Model modelANNReg0_1(frozenPathReg_cluster0_1);
    std::vector<std::string> operNameANN0_1;
    operNameANN0_1 = modelANNReg0_1.get_operations();
    Tensor inputPCAANNReg0_1{modelANNReg0_1, operNameANN0_1[0]}; //First tensor = Input (Standardized values)
    Tensor outputStandardizeANNReg0_1{modelANNReg0_1, operNameANN0_1[operNameANN0_1.size()-1]};    //Last tensor = Output (Standardized values)
    
    std::string frozenPathReg_cluster0_2 = casePath + cluster0_2_Folder + "frozen_bestModel.pb";
    Model modelANNReg0_2(frozenPathReg_cluster0_2);
    std::vector<std::string> operNameANN0_2;
    operNameANN0_2 = modelANNReg0_2.get_operations();
    Tensor inputPCAANNReg0_2{modelANNReg0_2, operNameANN0_2[0]}; //First tensor = Input (Standardized values)
    Tensor outputStandardizeANNReg0_2{modelANNReg0_2, operNameANN0_2[operNameANN0_2.size()-1]};    //Last tensor = Output (Standardized values)
    
    std::string frozenPathReg_cluster0_3 = casePath + cluster0_3_Folder + "frozen_bestModel.pb";
    Model modelANNReg0_3(frozenPathReg_cluster0_3);
    std::vector<std::string> operNameANN0_3;
    operNameANN0_3 = modelANNReg0_3.get_operations();
    Tensor inputPCAANNReg0_3{modelANNReg0_3, operNameANN0_3[0]}; //First tensor = Input (Standardized values)
    Tensor outputStandardizeANNReg0_3{modelANNReg0_3, operNameANN0_3[operNameANN0_3.size()-1]};    //Last tensor = Output (Standardized values)
    
    std::string frozenPathReg_cluster0_4 = casePath + cluster0_4_Folder + "frozen_bestModel.pb";
    Model modelANNReg0_4(frozenPathReg_cluster0_4);
    std::vector<std::string> operNameANN0_4;
    operNameANN0_4 = modelANNReg0_4.get_operations();
    Tensor inputPCAANNReg0_4{modelANNReg0_4, operNameANN0_4[0]}; //First tensor = Input (Standardized values)
    Tensor outputStandardizeANNReg0_4{modelANNReg0_4, operNameANN0_4[operNameANN0_4.size()-1]};    //Last tensor = Output (Standardized values)
    
    std::string frozenPathReg_cluster0_5 = casePath + cluster0_5_Folder + "frozen_bestModel.pb";
    Model modelANNReg0_5(frozenPathReg_cluster0_5);
    std::vector<std::string> operNameANN0_5;
    operNameANN0_5 = modelANNReg0_5.get_operations();
    Tensor inputPCAANNReg0_5{modelANNReg0_5, operNameANN0_5[0]}; //First tensor = Input (Standardized values)
    Tensor outputStandardizeANNReg0_5{modelANNReg0_5, operNameANN0_5[operNameANN0_5.size()-1]};    //Last tensor = Output (Standardized values)
    
    std::string frozenPathReg_cluster0_6 = casePath + cluster0_6_Folder + "frozen_bestModel.pb";
    Model modelANNReg0_6(frozenPathReg_cluster0_6);
    std::vector<std::string> operNameANN0_6;
    operNameANN0_6 = modelANNReg0_6.get_operations();
    Tensor inputPCAANNReg0_6{modelANNReg0_6, operNameANN0_6[0]}; //First tensor = Input (Standardized values)
    Tensor outputStandardizeANNReg0_6{modelANNReg0_6, operNameANN0_6[operNameANN0_6.size()-1]};    //Last tensor = Output (Standardized values)
    
    std::string frozenPathReg_cluster0_7 = casePath + cluster0_7_Folder + "frozen_bestModel.pb";
    Model modelANNReg0_7(frozenPathReg_cluster0_7);
    std::vector<std::string> operNameANN0_7;
    operNameANN0_7 = modelANNReg0_7.get_operations();
    Tensor inputPCAANNReg0_7{modelANNReg0_7, operNameANN0_7[0]}; //First tensor = Input (Standardized values)
    Tensor outputStandardizeANNReg0_7{modelANNReg0_7, operNameANN0_7[operNameANN0_7.size()-1]};    //Last tensor = Output (Standardized values)
    
    std::string frozenPathReg_cluster0_8 = casePath + cluster0_8_Folder + "frozen_bestModel.pb";
    Model modelANNReg0_8(frozenPathReg_cluster0_8);
    std::vector<std::string> operNameANN0_8;
    operNameANN0_8 = modelANNReg0_8.get_operations();
    Tensor inputPCAANNReg0_8{modelANNReg0_8, operNameANN0_8[0]}; //First tensor = Input (Standardized values)
    Tensor outputStandardizeANNReg0_8{modelANNReg0_8, operNameANN0_8[operNameANN0_8.size()-1]};    //Last tensor = Output (Standardized values)
    
    std::string frozenPathReg_cluster0_9 = casePath + cluster0_9_Folder + "frozen_bestModel.pb";
    Model modelANNReg0_9(frozenPathReg_cluster0_9);
    std::vector<std::string> operNameANN0_9;
    operNameANN0_9 = modelANNReg0_9.get_operations();
    Tensor inputPCAANNReg0_9{modelANNReg0_9, operNameANN0_9[0]}; //First tensor = Input (Standardized values)
    Tensor outputStandardizeANNReg0_9{modelANNReg0_9, operNameANN0_9[operNameANN0_9.size()-1]};    //Last tensor = Output (Standardized values)
    
    std::string frozenPathReg_cluster0_10 = casePath + cluster0_10_Folder + "frozen_bestModel.pb";
    Model modelANNReg0_10(frozenPathReg_cluster0_10);
    std::vector<std::string> operNameANN0_10;
    operNameANN0_10 = modelANNReg0_10.get_operations();
    Tensor inputPCAANNReg0_10{modelANNReg0_10, operNameANN0_10[0]}; //First tensor = Input (Standardized values)
    Tensor outputStandardizeANNReg0_10{modelANNReg0_10, operNameANN0_10[operNameANN0_10.size()-1]};    //Last tensor = Output (Standardized values)
    
    std::string frozenPathReg_cluster0_11 = casePath + cluster0_11_Folder + "frozen_bestModel.pb";
    Model modelANNReg0_11(frozenPathReg_cluster0_11);
    std::vector<std::string> operNameANN0_11;
    operNameANN0_11 = modelANNReg0_11.get_operations();
    Tensor inputPCAANNReg0_11{modelANNReg0_11, operNameANN0_11[0]}; //First tensor = Input (Standardized values)
    Tensor outputStandardizeANNReg0_11{modelANNReg0_11, operNameANN0_11[operNameANN0_11.size()-1]};    //Last tensor = Output (Standardized values)
    
    std::string frozenPathReg_cluster0_12 = casePath + cluster0_12_Folder + "frozen_bestModel.pb";
    Model modelANNReg0_12(frozenPathReg_cluster0_12);
    std::vector<std::string> operNameANN0_12;
    operNameANN0_12 = modelANNReg0_12.get_operations();
    Tensor inputPCAANNReg0_12{modelANNReg0_12, operNameANN0_12[0]}; //First tensor = Input (Standardized values)
    Tensor outputStandardizeANNReg0_12{modelANNReg0_12, operNameANN0_12[operNameANN0_12.size()-1]};    //Last tensor = Output (Standardized values)
    
    std::string frozenPathReg_cluster0_13 = casePath + cluster0_13_Folder + "frozen_bestModel.pb";
    Model modelANNReg0_13(frozenPathReg_cluster0_13);
    std::vector<std::string> operNameANN0_13;
    operNameANN0_13 = modelANNReg0_13.get_operations();
    Tensor inputPCAANNReg0_13{modelANNReg0_13, operNameANN0_13[0]}; //First tensor = Input (Standardized values)
    Tensor outputStandardizeANNReg0_13{modelANNReg0_13, operNameANN0_13[operNameANN0_13.size()-1]};    //Last tensor = Output (Standardized values)
    
    std::string frozenPathReg_cluster0_14 = casePath + cluster0_14_Folder + "frozen_bestModel.pb";
    Model modelANNReg0_14(frozenPathReg_cluster0_14);
    std::vector<std::string> operNameANN0_14;
    operNameANN0_14 = modelANNReg0_14.get_operations();
    Tensor inputPCAANNReg0_14{modelANNReg0_14, operNameANN0_14[0]}; //First tensor = Input (Standardized values)
    Tensor outputStandardizeANNReg0_14{modelANNReg0_14, operNameANN0_14[operNameANN0_14.size()-1]};    //Last tensor = Output (Standardized values)
    
    std::string frozenPathReg_cluster0_15 = casePath + cluster0_15_Folder + "frozen_bestModel.pb";
    Model modelANNReg0_15(frozenPathReg_cluster0_15);
    std::vector<std::string> operNameANN0_15;
    operNameANN0_15 = modelANNReg0_15.get_operations();
    Tensor inputPCAANNReg0_15{modelANNReg0_15, operNameANN0_15[0]}; //First tensor = Input (Standardized values)
    Tensor outputStandardizeANNReg0_15{modelANNReg0_15, operNameANN0_15[operNameANN0_15.size()-1]};    //Last tensor = Output (Standardized values)
    
    std::string frozenPathReg_cluster0_16 = casePath + cluster0_16_Folder + "frozen_bestModel.pb";
    Model modelANNReg0_16(frozenPathReg_cluster0_16);
    std::vector<std::string> operNameANN0_16;
    operNameANN0_16 = modelANNReg0_16.get_operations();
    Tensor inputPCAANNReg0_16{modelANNReg0_16, operNameANN0_16[0]}; //First tensor = Input (Standardized values)
    Tensor outputStandardizeANNReg0_16{modelANNReg0_16, operNameANN0_16[operNameANN0_16.size()-1]};    //Last tensor = Output (Standardized values)
    
    std::string frozenPathReg_cluster0_17 = casePath + cluster0_17_Folder + "frozen_bestModel.pb";
    Model modelANNReg0_17(frozenPathReg_cluster0_17);
    std::vector<std::string> operNameANN0_17;
    operNameANN0_17 = modelANNReg0_17.get_operations();
    Tensor inputPCAANNReg0_17{modelANNReg0_17, operNameANN0_17[0]}; //First tensor = Input (Standardized values)
    Tensor outputStandardizeANNReg0_17{modelANNReg0_17, operNameANN0_17[operNameANN0_17.size()-1]};    //Last tensor = Output (Standardized values)
    
    std::string frozenPathReg_cluster0_18 = casePath + cluster0_18_Folder + "frozen_bestModel.pb";
    Model modelANNReg0_18(frozenPathReg_cluster0_18);
    std::vector<std::string> operNameANN0_18;
    operNameANN0_18 = modelANNReg0_18.get_operations();
    Tensor inputPCAANNReg0_18{modelANNReg0_18, operNameANN0_18[0]}; //First tensor = Input (Standardized values)
    Tensor outputStandardizeANNReg0_18{modelANNReg0_18, operNameANN0_18[operNameANN0_18.size()-1]};    //Last tensor = Output (Standardized values)
    
    std::string frozenPathReg_cluster0_19 = casePath + cluster0_19_Folder + "frozen_bestModel.pb";
    Model modelANNReg0_19(frozenPathReg_cluster0_19);
    std::vector<std::string> operNameANN0_19;
    operNameANN0_19 = modelANNReg0_19.get_operations();
    Tensor inputPCAANNReg0_19{modelANNReg0_19, operNameANN0_19[0]}; //First tensor = Input (Standardized values)
    Tensor outputStandardizeANNReg0_19{modelANNReg0_19, operNameANN0_19[operNameANN0_19.size()-1]};    //Last tensor = Output (Standardized values)
    
    // parentCluster1 - childCluster 0 - grandChild
    std::string frozenPathReg_cluster1_0_0 = casePath + cluster1_0_0_Folder + "frozen_bestModel.pb";
    Model modelANNReg1_0_0(frozenPathReg_cluster1_0_0);
    std::vector<std::string> operNameANN1_0_0;
    operNameANN1_0_0 = modelANNReg1_0_0.get_operations();
    Tensor inputPCAANNReg1_0_0{modelANNReg1_0_0, operNameANN1_0_0[0]}; //First tensor = Input (Standardized values)
    Tensor outputStandardizeANNReg1_0_0{modelANNReg1_0_0, operNameANN1_0_0[operNameANN1_0_0.size()-1]};    //Last tensor = Output (Standardized values)
    
    std::string frozenPathReg_cluster1_0_1 = casePath + cluster1_0_1_Folder + "frozen_bestModel.pb";
    Model modelANNReg1_0_1(frozenPathReg_cluster1_0_1);
    std::vector<std::string> operNameANN1_0_1;
    operNameANN1_0_1 = modelANNReg1_0_1.get_operations();
    Tensor inputPCAANNReg1_0_1{modelANNReg1_0_1, operNameANN1_0_1[0]}; //First tensor = Input (Standardized values)
    Tensor outputStandardizeANNReg1_0_1{modelANNReg1_0_1, operNameANN1_0_1[operNameANN1_0_1.size()-1]};    //Last tensor = Output (Standardized values)
    
    std::string frozenPathReg_cluster1_0_2 = casePath + cluster1_0_2_Folder + "frozen_bestModel.pb";
    Model modelANNReg1_0_2(frozenPathReg_cluster1_0_2);
    std::vector<std::string> operNameANN1_0_2;
    operNameANN1_0_2 = modelANNReg1_0_2.get_operations();
    Tensor inputPCAANNReg1_0_2{modelANNReg1_0_2, operNameANN1_0_2[0]}; //First tensor = Input (Standardized values)
    Tensor outputStandardizeANNReg1_0_2{modelANNReg1_0_2, operNameANN1_0_2[operNameANN1_0_2.size()-1]};    //Last tensor = Output (Standardized values)
    
    std::string frozenPathReg_cluster1_0_3 = casePath + cluster1_0_3_Folder + "frozen_bestModel.pb";
    Model modelANNReg1_0_3(frozenPathReg_cluster1_0_3);
    std::vector<std::string> operNameANN1_0_3;
    operNameANN1_0_3 = modelANNReg1_0_3.get_operations();
    Tensor inputPCAANNReg1_0_3{modelANNReg1_0_3, operNameANN1_0_3[0]}; //First tensor = Input (Standardized values)
    Tensor outputStandardizeANNReg1_0_3{modelANNReg1_0_3, operNameANN1_0_3[operNameANN1_0_3.size()-1]};    //Last tensor = Output (Standardized values)
    
    std::string frozenPathReg_cluster1_0_4 = casePath + cluster1_0_4_Folder + "frozen_bestModel.pb";
    Model modelANNReg1_0_4(frozenPathReg_cluster1_0_4);
    std::vector<std::string> operNameANN1_0_4;
    operNameANN1_0_4 = modelANNReg1_0_4.get_operations();
    Tensor inputPCAANNReg1_0_4{modelANNReg1_0_4, operNameANN1_0_4[0]}; //First tensor = Input (Standardized values)
    Tensor outputStandardizeANNReg1_0_4{modelANNReg1_0_4, operNameANN1_0_4[operNameANN1_0_4.size()-1]};    //Last tensor = Output (Standardized values)
    
    std::string frozenPathReg_cluster1_0_5 = casePath + cluster1_0_5_Folder + "frozen_bestModel.pb";
    Model modelANNReg1_0_5(frozenPathReg_cluster1_0_5);
    std::vector<std::string> operNameANN1_0_5;
    operNameANN1_0_5 = modelANNReg1_0_5.get_operations();
    Tensor inputPCAANNReg1_0_5{modelANNReg1_0_5, operNameANN1_0_5[0]}; //First tensor = Input (Standardized values)
    Tensor outputStandardizeANNReg1_0_5{modelANNReg1_0_5, operNameANN1_0_5[operNameANN1_0_5.size()-1]};    //Last tensor = Output (Standardized values)
    
    std::string frozenPathReg_cluster1_0_6 = casePath + cluster1_0_6_Folder + "frozen_bestModel.pb";
    Model modelANNReg1_0_6(frozenPathReg_cluster1_0_6);
    std::vector<std::string> operNameANN1_0_6;
    operNameANN1_0_6 = modelANNReg1_0_6.get_operations();
    Tensor inputPCAANNReg1_0_6{modelANNReg1_0_6, operNameANN1_0_6[0]}; //First tensor = Input (Standardized values)
    Tensor outputStandardizeANNReg1_0_6{modelANNReg1_0_6, operNameANN1_0_6[operNameANN1_0_6.size()-1]};    //Last tensor = Output (Standardized values)
    
    std::string frozenPathReg_cluster1_0_7 = casePath + cluster1_0_7_Folder + "frozen_bestModel.pb";
    Model modelANNReg1_0_7(frozenPathReg_cluster1_0_7);
    std::vector<std::string> operNameANN1_0_7;
    operNameANN1_0_7 = modelANNReg1_0_7.get_operations();
    Tensor inputPCAANNReg1_0_7{modelANNReg1_0_7, operNameANN1_0_7[0]}; //First tensor = Input (Standardized values)
    Tensor outputStandardizeANNReg1_0_7{modelANNReg1_0_7, operNameANN1_0_7[operNameANN1_0_7.size()-1]};    //Last tensor = Output (Standardized values)
    
    std::string frozenPathReg_cluster1_0_8 = casePath + cluster1_0_8_Folder + "frozen_bestModel.pb";
    Model modelANNReg1_0_8(frozenPathReg_cluster1_0_8);
    std::vector<std::string> operNameANN1_0_8;
    operNameANN1_0_8 = modelANNReg1_0_8.get_operations();
    Tensor inputPCAANNReg1_0_8{modelANNReg1_0_8, operNameANN1_0_8[0]}; //First tensor = Input (Standardized values)
    Tensor outputStandardizeANNReg1_0_8{modelANNReg1_0_8, operNameANN1_0_8[operNameANN1_0_8.size()-1]};    //Last tensor = Output (Standardized values)
    
    std::string frozenPathReg_cluster1_0_9 = casePath + cluster1_0_9_Folder + "frozen_bestModel.pb";
    Model modelANNReg1_0_9(frozenPathReg_cluster1_0_9);
    std::vector<std::string> operNameANN1_0_9;
    operNameANN1_0_9 = modelANNReg1_0_9.get_operations();
    Tensor inputPCAANNReg1_0_9{modelANNReg1_0_9, operNameANN1_0_9[0]}; //First tensor = Input (Standardized values)
    Tensor outputStandardizeANNReg1_0_9{modelANNReg1_0_9, operNameANN1_0_9[operNameANN1_0_9.size()-1]};    //Last tensor = Output (Standardized values)
    
    std::string frozenPathReg_cluster1_0_10 = casePath + cluster1_0_10_Folder + "frozen_bestModel.pb";
    Model modelANNReg1_0_10(frozenPathReg_cluster1_0_10);
    std::vector<std::string> operNameANN1_0_10;
    operNameANN1_0_10 = modelANNReg1_0_10.get_operations();
    Tensor inputPCAANNReg1_0_10{modelANNReg1_0_10, operNameANN1_0_10[0]}; //First tensor = Input (Standardized values)
    Tensor outputStandardizeANNReg1_0_10{modelANNReg1_0_10, operNameANN1_0_10[operNameANN1_0_10.size()-1]};    //Last tensor = Output (Standardized values)
    
    std::string frozenPathReg_cluster1_0_11 = casePath + cluster1_0_11_Folder + "frozen_bestModel.pb";
    Model modelANNReg1_0_11(frozenPathReg_cluster1_0_11);
    std::vector<std::string> operNameANN1_0_11;
    operNameANN1_0_11 = modelANNReg1_0_11.get_operations();
    Tensor inputPCAANNReg1_0_11{modelANNReg1_0_11, operNameANN1_0_11[0]}; //First tensor = Input (Standardized values)
    Tensor outputStandardizeANNReg1_0_11{modelANNReg1_0_11, operNameANN1_0_11[operNameANN1_0_11.size()-1]};    //Last tensor = Output (Standardized values)
    
    std::string frozenPathReg_cluster1_0_12 = casePath + cluster1_0_12_Folder + "frozen_bestModel.pb";
    Model modelANNReg1_0_12(frozenPathReg_cluster1_0_12);
    std::vector<std::string> operNameANN1_0_12;
    operNameANN1_0_12 = modelANNReg1_0_12.get_operations();
    Tensor inputPCAANNReg1_0_12{modelANNReg1_0_12, operNameANN1_0_12[0]}; //First tensor = Input (Standardized values)
    Tensor outputStandardizeANNReg1_0_12{modelANNReg1_0_12, operNameANN1_0_12[operNameANN1_0_12.size()-1]};    //Last tensor = Output (Standardized values)
    
    std::string frozenPathReg_cluster1_0_13 = casePath + cluster1_0_13_Folder + "frozen_bestModel.pb";
    Model modelANNReg1_0_13(frozenPathReg_cluster1_0_13);
    std::vector<std::string> operNameANN1_0_13;
    operNameANN1_0_13 = modelANNReg1_0_13.get_operations();
    Tensor inputPCAANNReg1_0_13{modelANNReg1_0_13, operNameANN1_0_13[0]}; //First tensor = Input (Standardized values)
    Tensor outputStandardizeANNReg1_0_13{modelANNReg1_0_13, operNameANN1_0_13[operNameANN1_0_13.size()-1]};    //Last tensor = Output (Standardized values)
    
    std::string frozenPathReg_cluster1_0_14 = casePath + cluster1_0_14_Folder + "frozen_bestModel.pb";
    Model modelANNReg1_0_14(frozenPathReg_cluster1_0_14);
    std::vector<std::string> operNameANN1_0_14;
    operNameANN1_0_14 = modelANNReg1_0_14.get_operations();
    Tensor inputPCAANNReg1_0_14{modelANNReg1_0_14, operNameANN1_0_14[0]}; //First tensor = Input (Standardized values)
    Tensor outputStandardizeANNReg1_0_14{modelANNReg1_0_14, operNameANN1_0_14[operNameANN1_0_14.size()-1]};    //Last tensor = Output (Standardized values)
    
    std::string frozenPathReg_cluster1_0_15 = casePath + cluster1_0_15_Folder + "frozen_bestModel.pb";
    Model modelANNReg1_0_15(frozenPathReg_cluster1_0_15);
    std::vector<std::string> operNameANN1_0_15;
    operNameANN1_0_15 = modelANNReg1_0_15.get_operations();
    Tensor inputPCAANNReg1_0_15{modelANNReg1_0_15, operNameANN1_0_15[0]}; //First tensor = Input (Standardized values)
    Tensor outputStandardizeANNReg1_0_15{modelANNReg1_0_15, operNameANN1_0_15[operNameANN1_0_15.size()-1]};    //Last tensor = Output (Standardized values)
    
    std::string frozenPathReg_cluster1_0_16 = casePath + cluster1_0_16_Folder + "frozen_bestModel.pb";
    Model modelANNReg1_0_16(frozenPathReg_cluster1_0_16);
    std::vector<std::string> operNameANN1_0_16;
    operNameANN1_0_16 = modelANNReg1_0_16.get_operations();
    Tensor inputPCAANNReg1_0_16{modelANNReg1_0_16, operNameANN1_0_16[0]}; //First tensor = Input (Standardized values)
    Tensor outputStandardizeANNReg1_0_16{modelANNReg1_0_16, operNameANN1_0_16[operNameANN1_0_16.size()-1]};    //Last tensor = Output (Standardized values)
    
    std::string frozenPathReg_cluster1_0_17 = casePath + cluster1_0_17_Folder + "frozen_bestModel.pb";
    Model modelANNReg1_0_17(frozenPathReg_cluster1_0_17);
    std::vector<std::string> operNameANN1_0_17;
    operNameANN1_0_17 = modelANNReg1_0_17.get_operations();
    Tensor inputPCAANNReg1_0_17{modelANNReg1_0_17, operNameANN1_0_17[0]}; //First tensor = Input (Standardized values)
    Tensor outputStandardizeANNReg1_0_17{modelANNReg1_0_17, operNameANN1_0_17[operNameANN1_0_17.size()-1]};    //Last tensor = Output (Standardized values)
    
    std::string frozenPathReg_cluster1_0_18 = casePath + cluster1_0_18_Folder + "frozen_bestModel.pb";
    Model modelANNReg1_0_18(frozenPathReg_cluster1_0_18);
    std::vector<std::string> operNameANN1_0_18;
    operNameANN1_0_18 = modelANNReg1_0_18.get_operations();
    Tensor inputPCAANNReg1_0_18{modelANNReg1_0_18, operNameANN1_0_18[0]}; //First tensor = Input (Standardized values)
    Tensor outputStandardizeANNReg1_0_18{modelANNReg1_0_18, operNameANN1_0_18[operNameANN1_0_18.size()-1]};    //Last tensor = Output (Standardized values)
    
    std::string frozenPathReg_cluster1_0_19 = casePath + cluster1_0_19_Folder + "frozen_bestModel.pb";
    Model modelANNReg1_0_19(frozenPathReg_cluster1_0_19);
    std::vector<std::string> operNameANN1_0_19;
    operNameANN1_0_19 = modelANNReg1_0_19.get_operations();
    Tensor inputPCAANNReg1_0_19{modelANNReg1_0_19, operNameANN1_0_19[0]}; //First tensor = Input (Standardized values)
    Tensor outputStandardizeANNReg1_0_19{modelANNReg1_0_19, operNameANN1_0_19[operNameANN1_0_19.size()-1]};    //Last tensor = Output (Standardized values)
    

    // parentCluster1 - childCluster 1 - grandChild
    std::string frozenPathReg_cluster1_1_0 = casePath + cluster1_1_0_Folder + "frozen_bestModel.pb";
    Model modelANNReg1_1_0(frozenPathReg_cluster1_1_0);
    std::vector<std::string> operNameANN1_1_0;
    operNameANN1_1_0 = modelANNReg1_1_0.get_operations();
    Tensor inputPCAANNReg1_1_0{modelANNReg1_1_0, operNameANN1_1_0[0]}; //First tensor = Input (Standardized values)
    Tensor outputStandardizeANNReg1_1_0{modelANNReg1_1_0, operNameANN1_1_0[operNameANN1_1_0.size()-1]};    //Last tensor = Output (Standardized values)
    
    std::string frozenPathReg_cluster1_1_1 = casePath + cluster1_1_1_Folder + "frozen_bestModel.pb";
    Model modelANNReg1_1_1(frozenPathReg_cluster1_1_1);
    std::vector<std::string> operNameANN1_1_1;
    operNameANN1_1_1 = modelANNReg1_1_1.get_operations();
    Tensor inputPCAANNReg1_1_1{modelANNReg1_1_1, operNameANN1_1_1[0]}; //First tensor = Input (Standardized values)
    Tensor outputStandardizeANNReg1_1_1{modelANNReg1_1_1, operNameANN1_1_1[operNameANN1_1_1.size()-1]};    //Last tensor = Output (Standardized values)
    
    std::string frozenPathReg_cluster1_1_2 = casePath + cluster1_1_2_Folder + "frozen_bestModel.pb";
    Model modelANNReg1_1_2(frozenPathReg_cluster1_1_2);
    std::vector<std::string> operNameANN1_1_2;
    operNameANN1_1_2 = modelANNReg1_1_2.get_operations();
    Tensor inputPCAANNReg1_1_2{modelANNReg1_1_2, operNameANN1_1_2[0]}; //First tensor = Input (Standardized values)
    Tensor outputStandardizeANNReg1_1_2{modelANNReg1_1_2, operNameANN1_1_2[operNameANN1_1_2.size()-1]};    //Last tensor = Output (Standardized values)
    
    std::string frozenPathReg_cluster1_1_3 = casePath + cluster1_1_3_Folder + "frozen_bestModel.pb";
    Model modelANNReg1_1_3(frozenPathReg_cluster1_1_3);
    std::vector<std::string> operNameANN1_1_3;
    operNameANN1_1_3 = modelANNReg1_1_3.get_operations();
    Tensor inputPCAANNReg1_1_3{modelANNReg1_1_3, operNameANN1_1_3[0]}; //First tensor = Input (Standardized values)
    Tensor outputStandardizeANNReg1_1_3{modelANNReg1_1_3, operNameANN1_1_3[operNameANN1_1_3.size()-1]};    //Last tensor = Output (Standardized values)
    
    std::string frozenPathReg_cluster1_1_4 = casePath + cluster1_1_4_Folder + "frozen_bestModel.pb";
    Model modelANNReg1_1_4(frozenPathReg_cluster1_1_4);
    std::vector<std::string> operNameANN1_1_4;
    operNameANN1_1_4 = modelANNReg1_1_4.get_operations();
    Tensor inputPCAANNReg1_1_4{modelANNReg1_1_4, operNameANN1_1_4[0]}; //First tensor = Input (Standardized values)
    Tensor outputStandardizeANNReg1_1_4{modelANNReg1_1_4, operNameANN1_1_4[operNameANN1_1_4.size()-1]};    //Last tensor = Output (Standardized values)
    
    std::string frozenPathReg_cluster1_1_5 = casePath + cluster1_1_5_Folder + "frozen_bestModel.pb";
    Model modelANNReg1_1_5(frozenPathReg_cluster1_1_5);
    std::vector<std::string> operNameANN1_1_5;
    operNameANN1_1_5 = modelANNReg1_1_5.get_operations();
    Tensor inputPCAANNReg1_1_5{modelANNReg1_1_5, operNameANN1_1_5[0]}; //First tensor = Input (Standardized values)
    Tensor outputStandardizeANNReg1_1_5{modelANNReg1_1_5, operNameANN1_1_5[operNameANN1_1_5.size()-1]};    //Last tensor = Output (Standardized values)
    
    std::string frozenPathReg_cluster1_1_6 = casePath + cluster1_1_6_Folder + "frozen_bestModel.pb";
    Model modelANNReg1_1_6(frozenPathReg_cluster1_1_6);
    std::vector<std::string> operNameANN1_1_6;
    operNameANN1_1_6 = modelANNReg1_1_6.get_operations();
    Tensor inputPCAANNReg1_1_6{modelANNReg1_1_6, operNameANN1_1_6[0]}; //First tensor = Input (Standardized values)
    Tensor outputStandardizeANNReg1_1_6{modelANNReg1_1_6, operNameANN1_1_6[operNameANN1_1_6.size()-1]};    //Last tensor = Output (Standardized values)
    
    std::string frozenPathReg_cluster1_1_7 = casePath + cluster1_1_7_Folder + "frozen_bestModel.pb";
    Model modelANNReg1_1_7(frozenPathReg_cluster1_1_7);
    std::vector<std::string> operNameANN1_1_7;
    operNameANN1_1_7 = modelANNReg1_1_7.get_operations();
    Tensor inputPCAANNReg1_1_7{modelANNReg1_1_7, operNameANN1_1_7[0]}; //First tensor = Input (Standardized values)
    Tensor outputStandardizeANNReg1_1_7{modelANNReg1_1_7, operNameANN1_1_7[operNameANN1_1_7.size()-1]};    //Last tensor = Output (Standardized values)
    
    std::string frozenPathReg_cluster1_1_8 = casePath + cluster1_1_8_Folder + "frozen_bestModel.pb";
    Model modelANNReg1_1_8(frozenPathReg_cluster1_1_8);
    std::vector<std::string> operNameANN1_1_8;
    operNameANN1_1_8 = modelANNReg1_1_8.get_operations();
    Tensor inputPCAANNReg1_1_8{modelANNReg1_1_8, operNameANN1_1_8[0]}; //First tensor = Input (Standardized values)
    Tensor outputStandardizeANNReg1_1_8{modelANNReg1_1_8, operNameANN1_1_8[operNameANN1_1_8.size()-1]};    //Last tensor = Output (Standardized values)
    
    std::string frozenPathReg_cluster1_1_9 = casePath + cluster1_1_9_Folder + "frozen_bestModel.pb";
    Model modelANNReg1_1_9(frozenPathReg_cluster1_1_9);
    std::vector<std::string> operNameANN1_1_9;
    operNameANN1_1_9 = modelANNReg1_1_9.get_operations();
    Tensor inputPCAANNReg1_1_9{modelANNReg1_1_9, operNameANN1_1_9[0]}; //First tensor = Input (Standardized values)
    Tensor outputStandardizeANNReg1_1_9{modelANNReg1_1_9, operNameANN1_1_9[operNameANN1_1_9.size()-1]};    //Last tensor = Output (Standardized values)
    

    // parentCluster1 - childCluster 2 - grandChild
    std::string frozenPathReg_cluster1_2_0 = casePath + cluster1_2_0_Folder + "frozen_bestModel.pb";
    Model modelANNReg1_2_0(frozenPathReg_cluster1_2_0);
    std::vector<std::string> operNameANN1_2_0;
    operNameANN1_2_0 = modelANNReg1_2_0.get_operations();
    Tensor inputPCAANNReg1_2_0{modelANNReg1_2_0, operNameANN1_2_0[0]}; //First tensor = Input (Standardized values)
    Tensor outputStandardizeANNReg1_2_0{modelANNReg1_2_0, operNameANN1_2_0[operNameANN1_2_0.size()-1]};    //Last tensor = Output (Standardized values)
    
    std::string frozenPathReg_cluster1_2_1 = casePath + cluster1_2_1_Folder + "frozen_bestModel.pb";
    Model modelANNReg1_2_1(frozenPathReg_cluster1_2_1);
    std::vector<std::string> operNameANN1_2_1;
    operNameANN1_2_1 = modelANNReg1_2_1.get_operations();
    Tensor inputPCAANNReg1_2_1{modelANNReg1_2_1, operNameANN1_2_1[0]}; //First tensor = Input (Standardized values)
    Tensor outputStandardizeANNReg1_2_1{modelANNReg1_2_1, operNameANN1_2_1[operNameANN1_2_1.size()-1]};    //Last tensor = Output (Standardized values)
    
    std::string frozenPathReg_cluster1_2_2 = casePath + cluster1_2_2_Folder + "frozen_bestModel.pb";
    Model modelANNReg1_2_2(frozenPathReg_cluster1_2_2);
    std::vector<std::string> operNameANN1_2_2;
    operNameANN1_2_2 = modelANNReg1_2_2.get_operations();
    Tensor inputPCAANNReg1_2_2{modelANNReg1_2_2, operNameANN1_2_2[0]}; //First tensor = Input (Standardized values)
    Tensor outputStandardizeANNReg1_2_2{modelANNReg1_2_2, operNameANN1_2_2[operNameANN1_2_2.size()-1]};    //Last tensor = Output (Standardized values)
    
    std::string frozenPathReg_cluster1_2_3 = casePath + cluster1_2_3_Folder + "frozen_bestModel.pb";
    Model modelANNReg1_2_3(frozenPathReg_cluster1_2_3);
    std::vector<std::string> operNameANN1_2_3;
    operNameANN1_2_3 = modelANNReg1_2_3.get_operations();
    Tensor inputPCAANNReg1_2_3{modelANNReg1_2_3, operNameANN1_2_3[0]}; //First tensor = Input (Standardized values)
    Tensor outputStandardizeANNReg1_2_3{modelANNReg1_2_3, operNameANN1_2_3[operNameANN1_2_3.size()-1]};    //Last tensor = Output (Standardized values)
    
    std::string frozenPathReg_cluster1_2_4 = casePath + cluster1_2_4_Folder + "frozen_bestModel.pb";
    Model modelANNReg1_2_4(frozenPathReg_cluster1_2_4);
    std::vector<std::string> operNameANN1_2_4;
    operNameANN1_2_4 = modelANNReg1_2_4.get_operations();
    Tensor inputPCAANNReg1_2_4{modelANNReg1_2_4, operNameANN1_2_4[0]}; //First tensor = Input (Standardized values)
    Tensor outputStandardizeANNReg1_2_4{modelANNReg1_2_4, operNameANN1_2_4[operNameANN1_2_4.size()-1]};    //Last tensor = Output (Standardized values)
    
    std::string frozenPathReg_cluster1_2_5 = casePath + cluster1_2_5_Folder + "frozen_bestModel.pb";
    Model modelANNReg1_2_5(frozenPathReg_cluster1_2_5);
    std::vector<std::string> operNameANN1_2_5;
    operNameANN1_2_5 = modelANNReg1_2_5.get_operations();
    Tensor inputPCAANNReg1_2_5{modelANNReg1_2_5, operNameANN1_2_5[0]}; //First tensor = Input (Standardized values)
    Tensor outputStandardizeANNReg1_2_5{modelANNReg1_2_5, operNameANN1_2_5[operNameANN1_2_5.size()-1]};    //Last tensor = Output (Standardized values)
    
    std::string frozenPathReg_cluster1_2_6 = casePath + cluster1_2_6_Folder + "frozen_bestModel.pb";
    Model modelANNReg1_2_6(frozenPathReg_cluster1_2_6);
    std::vector<std::string> operNameANN1_2_6;
    operNameANN1_2_6 = modelANNReg1_2_6.get_operations();
    Tensor inputPCAANNReg1_2_6{modelANNReg1_2_6, operNameANN1_2_6[0]}; //First tensor = Input (Standardized values)
    Tensor outputStandardizeANNReg1_2_6{modelANNReg1_2_6, operNameANN1_2_6[operNameANN1_2_6.size()-1]};    //Last tensor = Output (Standardized values)
    
    std::string frozenPathReg_cluster1_2_7 = casePath + cluster1_2_7_Folder + "frozen_bestModel.pb";
    Model modelANNReg1_2_7(frozenPathReg_cluster1_2_7);
    std::vector<std::string> operNameANN1_2_7;
    operNameANN1_2_7 = modelANNReg1_2_7.get_operations();
    Tensor inputPCAANNReg1_2_7{modelANNReg1_2_7, operNameANN1_2_7[0]}; //First tensor = Input (Standardized values)
    Tensor outputStandardizeANNReg1_2_7{modelANNReg1_2_7, operNameANN1_2_7[operNameANN1_2_7.size()-1]};    //Last tensor = Output (Standardized values)
    
    std::string frozenPathReg_cluster1_2_8 = casePath + cluster1_2_8_Folder + "frozen_bestModel.pb";
    Model modelANNReg1_2_8(frozenPathReg_cluster1_2_8);
    std::vector<std::string> operNameANN1_2_8;
    operNameANN1_2_8 = modelANNReg1_2_8.get_operations();
    Tensor inputPCAANNReg1_2_8{modelANNReg1_2_8, operNameANN1_2_8[0]}; //First tensor = Input (Standardized values)
    Tensor outputStandardizeANNReg1_2_8{modelANNReg1_2_8, operNameANN1_2_8[operNameANN1_2_8.size()-1]};    //Last tensor = Output (Standardized values)
    
    std::string frozenPathReg_cluster1_2_9 = casePath + cluster1_2_9_Folder + "frozen_bestModel.pb";
    Model modelANNReg1_2_9(frozenPathReg_cluster1_2_9);
    std::vector<std::string> operNameANN1_2_9;
    operNameANN1_2_9 = modelANNReg1_2_9.get_operations();
    Tensor inputPCAANNReg1_2_9{modelANNReg1_2_9, operNameANN1_2_9[0]}; //First tensor = Input (Standardized values)
    Tensor outputStandardizeANNReg1_2_9{modelANNReg1_2_9, operNameANN1_2_9[operNameANN1_2_9.size()-1]};    //Last tensor = Output (Standardized values)
    


    // Create the temporary vectors to store the intermediary result of each particle
    //Create vectors to GLOBAL standardize Y T from Y (without N2)  and T of particle >> Kmeans or ANN Classifier
    std::vector<float> input_GlobalStandardized(numVarANN); //Expected 11 type of data : 10Y T
    
    //Create vectors to PARENT LOCAL standardize Y T from Y (without N2)  and T of particle >> Kmeans or ANN Classifier
    std::vector<float> input_parentLocalStandardized(numVarANN); //Expected 11 type of data : 10Y T
    
    //Create vectors to CHILD LOCAL standardize Y T from Y (without N2)  and T of particle >> ANNRegression or Kmeans (parentCluster0)
    std::vector<float> input_childLocalStandardized(numVarANN); //Expected 11 type of data : 10Y T
    
    //Create vectors to PARENT LOCAL standardize Y T from Y (without N2)  and T of particle >> ANNRegression (parentCluster1)
    std::vector<float> input_grandChildLocalStandardized(numVarANN); //Expected 11 type of data : 10Y T
    
    //Create vector to store inputPCA (inputANN in PCA space)
    std::vector<float> inputPCA_Vec(numComponentPCA);
    
    //Create vector to store outputANNStandardize_Vec (outputANN Standardize in composition space)
    std::vector<float> outputStandardizeANN_Vec(numVarANN);
    
    
    // Commun variable to load
    vector<float> child_MeanINPUTLoad = {};
    vector<float> child_StdINPUTLoad = {};
    vector<float> child_MeanOUTPUTLoad = {};
    vector<float> child_StdOUTPUTLoad = {};
    
    vector<float> child_MaxOUTPUTLoad = {}; //Max Variation (output ANN - composition space)
    vector<float> child_MinOUTPUTLoad = {}; //Min Variation (output ANN - composition space)
    
    vector<double> maxAftREACTORLoad = {}; //double because std::min need double
    vector<double> minAftREACTORLoad = {}; // Max/Min Y,T inside a cluster (composition space)
    
/* Huu-Tri@20211006: Commented - ANN without PCA
    //Create matrices for PCA transform (input & output)
    cv::PCA pcaANN;
    cv::Mat input_childLocalStandardized_Mat; // Mat: Composition space = input of PCA transform
    cv::Mat input_grandChildLocalStandardized_Mat; // Mat: Composition space = input of PCA transform
    cv::Mat inputPCA;	// Mat: LPCA space = output of PCA transform = input of ANN in LPCA space
*/ //Huu-Tri@20211006: Commented - ANN without PCA
    
    // To count EMST mixed particles
    int modifEMSTParticle = 0; // To count the number of particle modified by EMST
    
    // GLOBAL Standardize vector
    cv::Mat x;
    
    // Distance to parent Cluster
    double d0 = 0;
    double d1 = 0;
    
    // minLoc for parentCluster
    int indexParentCluster;
    cv::Point minLoc;
    
    // LOCAL PARENT Standardize vector
    cv::Mat xParent;
    
    // minLoc for childCluster
    cv::Point minLocParent;
    int indexChildCluster;
    
    // LOCAL CHILD Standardize vector (childCluster1_0, 1_1, 1_2)
    cv::Mat xChild;
    
    // minLoc for childCluster
    cv::Point minLocChild;
    int indexGrandChildCluster;
    
    // Distance to childCluster
    //parentCluster0
    double d0_0 = 0; // Calculate Euclidean distance from x to cluster 0
    double d0_1 = 0;
    double d0_2 = 0;
    double d0_3 = 0;
    double d0_4 = 0;
    double d0_5 = 0;
    double d0_6 = 0;
    double d0_7 = 0;
    double d0_8 = 0;
    double d0_9 = 0;
    double d0_10 = 0;
    double d0_11 = 0;
    double d0_12 = 0;
    double d0_13 = 0;
    double d0_14 = 0;
    double d0_15 = 0;
    double d0_16 = 0;
    double d0_17 = 0;
    double d0_18 = 0;
    double d0_19 = 0;
    
    //parentCluster1
    double d1_0 = 0; // Calculate Euclidean distance from x to cluster 0
    double d1_1 = 0;
    double d1_2 = 0;
    
    //parentCluster1 - childCluster0
    double d1_0_0 = 0;
    double d1_0_1 = 0;
    double d1_0_2 = 0;
    double d1_0_3 = 0;
    double d1_0_4 = 0;
    double d1_0_5 = 0;
    double d1_0_6 = 0;
    double d1_0_7 = 0;
    double d1_0_8 = 0;
    double d1_0_9 = 0;
    double d1_0_10 = 0;
    double d1_0_11 = 0;
    double d1_0_12 = 0;
    double d1_0_13 = 0;
    double d1_0_14 = 0;
    double d1_0_15 = 0;
    double d1_0_16 = 0;
    double d1_0_17 = 0;
    double d1_0_18 = 0;
    double d1_0_19 = 0;
    double d1_0_20 = 0;
    double d1_0_21 = 0;
    double d1_0_22 = 0;
    double d1_0_23 = 0;
    double d1_0_24 = 0;
    double d1_0_25 = 0;
    double d1_0_26 = 0;
    double d1_0_27 = 0;
    double d1_0_28 = 0;
    double d1_0_29 = 0;
    double d1_0_30 = 0;
    double d1_0_31 = 0;
    double d1_0_32 = 0;
    double d1_0_33 = 0;
    double d1_0_34 = 0;
    double d1_0_35 = 0;
    double d1_0_36 = 0;
    double d1_0_37 = 0;
    double d1_0_38 = 0;
    double d1_0_39 = 0;
    
    //parentCluster1 - childCluster1
    double d1_1_0 = 0;
    double d1_1_1 = 0;
    double d1_1_2 = 0;
    double d1_1_3 = 0;
    double d1_1_4 = 0;
    double d1_1_5 = 0;
    double d1_1_6 = 0;
    double d1_1_7 = 0;
    double d1_1_8 = 0;
    double d1_1_9 = 0;
    double d1_1_10 = 0;
    double d1_1_11 = 0;
    double d1_1_12 = 0;
    double d1_1_13 = 0;
    double d1_1_14 = 0;
    
    //parentCluster1 - childCluster2
    double d1_2_0 = 0;
    double d1_2_1 = 0;
    double d1_2_2 = 0;
    double d1_2_3 = 0;
    double d1_2_4 = 0;
    double d1_2_5 = 0;
    double d1_2_6 = 0;
    double d1_2_7 = 0;
    double d1_2_8 = 0;
    double d1_2_9 = 0;
    double d1_2_10 = 0;
    double d1_2_11 = 0;
    double d1_2_12 = 0;
    double d1_2_13 = 0;
    double d1_2_14 = 0;
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
    string Case = "120_EMST_ORChData.FourUMONS_NewVersionORCh_Stochastic_8Juillet2020_Qx10_NewConditions_outletGB_GRI30_ObtainFullDataParticle/"; // Be careful, add / at the end

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
     bool flagHeatLossStochas = false;
     double T_wall = 1600.0;
     double alpha_loss = 0.0;	//1.5e+06;
     double beta_heat = 0.0;
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
    bool flagANN = true;

     bool activateIndexCluster = true;	//flag to activate Full data (write full_dataParticle.dat)
     if(flagANN && activateIndexCluster && i==0)
     {
	if(rank==0)	//Parallel MPI stuff to prevent print multiple lines, rank = 0 <=> First processor
	{
		cout << "*indexClusterFile saving is activated" << endl;
		// Check if meanSpeciesProdRate.txt file exists at the first step
		// If yes, clear the file content
		if(file_exists("outputs/indexClusterFile.dat"))
		{
		
			cout << " -------------- Warrning --------------" << endl;
			cout << "outputs/indexClusterFile.dat exists. Clearing file ... " << endl;	
			
			ofstream indexClusterFile_clear("outputs/indexClusterFile.dat", ios::out | ios::trunc);	//open file in trunc mode to clear the content
			indexClusterFile_clear.close(); //close the file
	
		}	
		
		ofstream indexClusterFile("outputs/indexCluster.dat",ios::app); //ios::app = append at the end of the file
		if(indexClusterFile)
		{ 
   			indexClusterFile << "#Time	";
			indexClusterFile << "Particle_number	";
            indexClusterFile << "parentCluster	";
            indexClusterFile << "childCluster	";
			indexClusterFile << "grandChildCluster" << endl;
            
		}
        indexClusterFile.close();
	} //END if(rank==0)
     } //ENd  if(activateIndexCluster) 

    //if(i>0) // Mix first step
    //{
	//flagANN = true;
    //}
    if (flagANN == true) //Initialization
    {

	// The SavedModel format (a folder) should be frozen to a single .pb file
	// Tensorflow v2 : https://leimao.github.io/blog/Save-Load-Inference-From-TF2-Frozen-Graph/ 
		
	// Collecte the maximum values of Y H T on entire full_dataParticle to normalize
		// Copy these vectors from ORCh_ANN_train.ipynb
		// ATTENTION : Must change these vectors when using another model (delta_t and conditions are differents between each case)
		// Species order is the same as ORCh_ANN_train and also the reduced scheme: Have to check !!!!
	
    //Declare parameters above before LOOP i
        
        // Create kRun
        int kRun=0;
        
        
     // BEGIN the loop for all particles
        for (int p=0; p<nTot; p++)
        {
            
            if(state_emst[p]>0) // EMST mixed particles count
            {
                modifEMSTParticle++;
            }
          
        // GLOBAL Standardize Y ,T: standardizedX = (X - meanScale) / standardDeviationScale
            // Standardize T
            input_GlobalStandardized[numVarANN-1] = listParticles[p]->m_T_gas - GLOBAL_meanScale[numVarANN-1];
            input_GlobalStandardized[numVarANN-1] /= GLOBAL_stdScale[numVarANN-1];
            // Standardize Y
            for (int kANN=0; kANN<nsp-1; kANN++) // Don't take N2 (last species)
            {
                input_GlobalStandardized[kANN] = listParticles[p]->m_Yk_gas[kANN] - GLOBAL_meanScale[kANN];
                input_GlobalStandardized[kANN] /= GLOBAL_stdScale[kANN];
                //Cross-check if there is N2 inside, if yes, should move N2 to the end of the scheme
                if(mixture->speciesName(kANN) == "N2") {cout << "Warning : N2 in the ANN!!!" << endl;}
            }

        // CLASSIFIER GLOBAL: Caculate the distance : Eucledian (Kmeans) or Projection parition (LPCA)
            // >> Determine indexParentCluster
            // LPCA example: cv::Mat d0 = (x-r0).t() * eigenVec0.t() * eigenVec0 * (x-r0);
            // Eucledean: Can be calculated by std::vector
            
            // Declare vector x
            x = cv::Mat(numVarANN, 1, CV_32F, input_GlobalStandardized.data());
            d0 = cv::norm(x,r0,cv::NORM_L2SQR); // Calculate Euclidean distance from x to cluster 0
            d1 = cv::norm(x,r1,cv::NORM_L2SQR); // NORM_L2SQRSQR can be used to calculate distance without square root => accelerate the calculation
            // Group distance to d
            cv::Mat d;
            d.push_back(d0);
            d.push_back(d1);
            // Take position of mean (argmin) = indexParentCluster
            cv::minMaxLoc(d, NULL, NULL, &minLoc, NULL); // d, minVal, maxVal, minLoc, maxLoc
            indexParentCluster = minLoc.y; //Point has 2 coordinates (x,y)
            
            if(indexParentCluster==0) //parentCluster0 - 14 child clusters
            {

                // LOCAL PARENT Standardize Y ,T for Kmeans LOCAL PARENT: standardizedX = (X - meanScale) / standardDeviationScale
                    // Standardize T
                input_parentLocalStandardized[numVarANN-1] = listParticles[p]->m_T_gas - parentCluster0_meanScaleINPUT[numVarANN-1];
                input_parentLocalStandardized[numVarANN-1] /= parentCluster0_stdScaleINPUT[numVarANN-1];
                    // Standardize Y
                for (int kANN=0; kANN<nsp-1; kANN++) // Don't take N2 (last species)
                {
                    input_parentLocalStandardized[kANN] = listParticles[p]->m_Yk_gas[kANN] - parentCluster0_meanScaleINPUT[kANN];
                    input_parentLocalStandardized[kANN] /= parentCluster0_stdScaleINPUT[kANN];
                    //Cross-check if there is N2 inside, if yes, should move N2 to the end of the scheme
                    if(mixture->speciesName(kANN) == "N2") {cout << "Warning : N2 in the ANN!!!" << endl;}
                }
                
                // CLASSIFIER LOCAL PARENT: Caculate the distance : Eucledian (Kmeans) or Projection parition (LPCA)
                    // >> Determine indexChildCluster
                    // LPCA example: cv::Mat d0 = (x-r0).t() * eigenVec0.t() * eigenVec0 * (x-r0);
                    // Eucledean: Can be calculated by std::vector
                
                    // Declare vector x
                xParent = cv::Mat(numVarANN, 1, CV_32F, input_parentLocalStandardized.data());
                d0_0 = cv::norm(xParent,r0_0,cv::NORM_L2SQR); // Calculate Euclidean distance from x to cluster 0
                d0_1 = cv::norm(xParent,r0_1,cv::NORM_L2SQR);
                d0_2 = cv::norm(xParent,r0_2,cv::NORM_L2SQR);
                d0_3 = cv::norm(xParent,r0_3,cv::NORM_L2SQR);
                d0_4 = cv::norm(xParent,r0_4,cv::NORM_L2SQR);
                d0_5 = cv::norm(xParent,r0_5,cv::NORM_L2SQR);
                d0_6 = cv::norm(xParent,r0_6,cv::NORM_L2SQR);
                d0_7 = cv::norm(xParent,r0_7,cv::NORM_L2SQR);
                d0_8 = cv::norm(xParent,r0_8,cv::NORM_L2SQR);
                d0_9 = cv::norm(xParent,r0_9,cv::NORM_L2SQR);
                d0_10 = cv::norm(xParent,r0_10,cv::NORM_L2SQR);
                d0_11 = cv::norm(xParent,r0_11,cv::NORM_L2SQR);
                d0_12 = cv::norm(xParent,r0_12,cv::NORM_L2SQR);
                d0_13 = cv::norm(xParent,r0_13,cv::NORM_L2SQR);
                d0_14 = cv::norm(xParent,r0_14,cv::NORM_L2SQR);
                d0_15 = cv::norm(xParent,r0_15,cv::NORM_L2SQR);
                d0_16 = cv::norm(xParent,r0_16,cv::NORM_L2SQR);
                d0_17 = cv::norm(xParent,r0_17,cv::NORM_L2SQR);
                d0_18 = cv::norm(xParent,r0_18,cv::NORM_L2SQR);
                d0_19 = cv::norm(xParent,r0_19,cv::NORM_L2SQR);
                // Group distance to d
                cv::Mat dParent;
                dParent.push_back(d0_0);
                dParent.push_back(d0_1);
                dParent.push_back(d0_2);
                dParent.push_back(d0_3);
                dParent.push_back(d0_4);
                dParent.push_back(d0_5);
                dParent.push_back(d0_6);
                dParent.push_back(d0_7);
                dParent.push_back(d0_8);
                dParent.push_back(d0_9);
                dParent.push_back(d0_10);
                dParent.push_back(d0_11);
                dParent.push_back(d0_12);
                dParent.push_back(d0_13);
                dParent.push_back(d0_14);
                dParent.push_back(d0_15);
                dParent.push_back(d0_16);
                dParent.push_back(d0_17);
                dParent.push_back(d0_18);
                dParent.push_back(d0_19);
                    // Take position of mean (argmin) = indexParentCluster
                cv::minMaxLoc(dParent, NULL, NULL, &minLocParent, NULL); // d, minVal, maxVal, minLoc, maxLoc
                indexChildCluster = minLocParent.y; //Point has 2 coordinates (x,y)

                if(rank==0) // Print to check
                {
                    cout << "Particles " << p << " || parentCluster = " << indexParentCluster << " || childCluster = " << indexChildCluster << " ";
                }
                
                // Save indexCluster of this particle
                if(activateIndexCluster)
                {
                    if(rank==0)
                    {
                        ofstream indexClusterFile("outputs/indexCluster.dat",ios::app); //ios::app = append at the end of the file
                        indexClusterFile << t << "	";
                        indexClusterFile << p << "	";
                        indexClusterFile << indexParentCluster <<  "	";
                        indexClusterFile << indexChildCluster <<  "	";
                        indexClusterFile << "" << endl; // parentCluster0 has not grand child
                        indexClusterFile.close();
                    }
                }

 
                
                switch(indexChildCluster) //parentCluster0
                {
                    // Cluster 0_0
                    case 0:
                        if(rank==0)
                        {
                            cout << " Case0_0"<<  endl;
                        }
                        
                        advanceANN_withoutPCA(cluster0_0_meanScaleINPUT, cluster0_0_stdScaleINPUT,
                                   cluster0_0_meanScaleOUTPUT, cluster0_0_stdScaleOUTPUT,
                                   cluster0_0_maxOUTPUT, cluster0_0_minOUTPUT,
                                   child0_0_maxVal, child0_0_minVal,
                                   inputPCAANNReg0_0, modelANNReg0_0, outputStandardizeANNReg0_0,
                                   listParticles[p]->m_T_gas, listParticles[p]->m_Yk_gas,
                                   numVarANN, nsp,
                                   input_childLocalStandardized, 
                                   outputStandardizeANN_Vec, mixture);
                        break;
                        
                        
                    // Cluster 0_1
                    case 1:
                        if(rank==0)
                        {
                            cout << " Case0_1"<<  endl;
                        }
                        advanceANN_withoutPCA(cluster0_1_meanScaleINPUT, cluster0_1_stdScaleINPUT,
                                   cluster0_1_meanScaleOUTPUT, cluster0_1_stdScaleOUTPUT,
                                   cluster0_1_maxOUTPUT, cluster0_1_minOUTPUT,
                                   child0_1_maxVal, child0_1_minVal,
                                   inputPCAANNReg0_1, modelANNReg0_1, outputStandardizeANNReg0_1,
                                   listParticles[p]->m_T_gas, listParticles[p]->m_Yk_gas,
                                   numVarANN, nsp,
                                   input_childLocalStandardized, 
                                   outputStandardizeANN_Vec, mixture);
                        break;
                        
                    // Cluster 0_2
                    case 2:
                        if(rank==0)
                        {
                            cout << " Case0_2"<<  endl;
                        }
                        advanceANN_withoutPCA(cluster0_2_meanScaleINPUT, cluster0_2_stdScaleINPUT,
                                   cluster0_2_meanScaleOUTPUT, cluster0_2_stdScaleOUTPUT,
                                   cluster0_2_maxOUTPUT, cluster0_2_minOUTPUT,
                                   child0_2_maxVal, child0_2_minVal,
                                   inputPCAANNReg0_2, modelANNReg0_2, outputStandardizeANNReg0_2,
                                   listParticles[p]->m_T_gas, listParticles[p]->m_Yk_gas,
                                   numVarANN, nsp,
                                   input_childLocalStandardized, 
                                   outputStandardizeANN_Vec, mixture);
                        break;
                        
                    // Cluster 0_3
                    case 3:
                        if(rank==0)
                        {
                            cout << " Case0_3"<<  endl;
                        }
                        advanceANN_withoutPCA(cluster0_3_meanScaleINPUT, cluster0_3_stdScaleINPUT,
                                   cluster0_3_meanScaleOUTPUT, cluster0_3_stdScaleOUTPUT,
                                   cluster0_3_maxOUTPUT, cluster0_3_minOUTPUT,
                                   child0_3_maxVal, child0_3_minVal,
                                   inputPCAANNReg0_3, modelANNReg0_3, outputStandardizeANNReg0_3,
                                   listParticles[p]->m_T_gas, listParticles[p]->m_Yk_gas,
                                   numVarANN, nsp,
                                   input_childLocalStandardized, 
                                   outputStandardizeANN_Vec, mixture);
                        break;

                    // Cluster 0_4
                    case 4:
                        if(rank==0)
                        {
                            cout << " Case0_4"<<  endl;
                        }
                        advanceANN_withoutPCA(cluster0_4_meanScaleINPUT, cluster0_4_stdScaleINPUT,
                                   cluster0_4_meanScaleOUTPUT, cluster0_4_stdScaleOUTPUT,
                                   cluster0_4_maxOUTPUT, cluster0_4_minOUTPUT,
                                   child0_4_maxVal, child0_4_minVal,
                                   inputPCAANNReg0_4, modelANNReg0_4, outputStandardizeANNReg0_4,
                                   listParticles[p]->m_T_gas, listParticles[p]->m_Yk_gas,
                                   numVarANN, nsp,
                                   input_childLocalStandardized, 
                                   outputStandardizeANN_Vec, mixture);
                        break;
                        
                    // Cluster 0_5
                    case 5:
                        if(rank==0)
                        {
                            cout << " Case0_5"<<  endl;
                        }
                        advanceANN_withoutPCA(cluster0_5_meanScaleINPUT, cluster0_5_stdScaleINPUT,
                                   cluster0_5_meanScaleOUTPUT, cluster0_5_stdScaleOUTPUT,
                                   cluster0_5_maxOUTPUT, cluster0_5_minOUTPUT,
                                   child0_5_maxVal, child0_5_minVal,
                                   inputPCAANNReg0_5, modelANNReg0_5, outputStandardizeANNReg0_5,
                                   listParticles[p]->m_T_gas, listParticles[p]->m_Yk_gas,
                                   numVarANN, nsp,
                                   input_childLocalStandardized, 
                                   outputStandardizeANN_Vec, mixture);
                        break;
                        
                    // Cluster 0_6
                    case 6:
                        if(rank==0)
                        {
                            cout << " Case0_6"<<  endl;
                        }
                        advanceANN_withoutPCA(cluster0_6_meanScaleINPUT, cluster0_6_stdScaleINPUT,
                                   cluster0_6_meanScaleOUTPUT, cluster0_6_stdScaleOUTPUT,
                                   cluster0_6_maxOUTPUT, cluster0_6_minOUTPUT,
                                   child0_6_maxVal, child0_6_minVal,
                                   inputPCAANNReg0_6, modelANNReg0_6, outputStandardizeANNReg0_6,
                                   listParticles[p]->m_T_gas, listParticles[p]->m_Yk_gas,
                                   numVarANN, nsp,
                                   input_childLocalStandardized, 
                                   outputStandardizeANN_Vec, mixture);
                        break;
                        
                    // Cluster 0_7
                    case 7:
                        if(rank==0)
                        {
                            cout << " Case0_7"<<  endl;
                        }
                        advanceANN_withoutPCA(cluster0_7_meanScaleINPUT, cluster0_7_stdScaleINPUT,
                                   cluster0_7_meanScaleOUTPUT, cluster0_7_stdScaleOUTPUT,
                                   cluster0_7_maxOUTPUT, cluster0_7_minOUTPUT,
                                   child0_7_maxVal, child0_7_minVal,
                                   inputPCAANNReg0_7, modelANNReg0_7, outputStandardizeANNReg0_7,
                                   listParticles[p]->m_T_gas, listParticles[p]->m_Yk_gas,
                                   numVarANN, nsp,
                                   input_childLocalStandardized, 
                                   outputStandardizeANN_Vec, mixture);
                        break;
                        
                    // Cluster 0_8
                    case 8:
                        if(rank==0)
                        {
                            cout << " Case0_8"<<  endl;
                        }
                        advanceANN_withoutPCA(cluster0_8_meanScaleINPUT, cluster0_8_stdScaleINPUT,
                                   cluster0_8_meanScaleOUTPUT, cluster0_8_stdScaleOUTPUT,
                                   cluster0_8_maxOUTPUT, cluster0_8_minOUTPUT,
                                   child0_8_maxVal, child0_8_minVal,
                                   inputPCAANNReg0_8, modelANNReg0_8, outputStandardizeANNReg0_8,
                                   listParticles[p]->m_T_gas, listParticles[p]->m_Yk_gas,
                                   numVarANN, nsp,
                                   input_childLocalStandardized, 
                                   outputStandardizeANN_Vec, mixture);
                        break;
                        
                    // Cluster 0_9
                    case 9:
                        if(rank==0)
                        {
                            cout << " Case0_9"<<  endl;
                        }
                        advanceANN_withoutPCA(cluster0_9_meanScaleINPUT, cluster0_9_stdScaleINPUT,
                                   cluster0_9_meanScaleOUTPUT, cluster0_9_stdScaleOUTPUT,
                                   cluster0_9_maxOUTPUT, cluster0_9_minOUTPUT,
                                   child0_9_maxVal, child0_9_minVal,
                                   inputPCAANNReg0_9, modelANNReg0_9, outputStandardizeANNReg0_9,
                                   listParticles[p]->m_T_gas, listParticles[p]->m_Yk_gas,
                                   numVarANN, nsp,
                                   input_childLocalStandardized, 
                                   outputStandardizeANN_Vec, mixture);
                        break;
                        
                    // Cluster 0_10
                    case 10:
                        if(rank==0)
                        {
                            cout << " Case0_10"<<  endl;
                        }
                        advanceANN_withoutPCA(cluster0_10_meanScaleINPUT, cluster0_10_stdScaleINPUT,
                                   cluster0_10_meanScaleOUTPUT, cluster0_10_stdScaleOUTPUT,
                                   cluster0_10_maxOUTPUT, cluster0_10_minOUTPUT,
                                   child0_10_maxVal, child0_10_minVal,
                                   inputPCAANNReg0_10, modelANNReg0_10, outputStandardizeANNReg0_10,
                                   listParticles[p]->m_T_gas, listParticles[p]->m_Yk_gas,
                                   numVarANN, nsp,
                                   input_childLocalStandardized, 
                                   outputStandardizeANN_Vec, mixture);
                        break;
                        
                    // Cluster 0_11
                    case 11:
                        if(rank==0)
                        {
                            cout << " Case0_11"<<  endl;
                        }
                        advanceANN_withoutPCA(cluster0_11_meanScaleINPUT, cluster0_11_stdScaleINPUT,
                                   cluster0_11_meanScaleOUTPUT, cluster0_11_stdScaleOUTPUT,
                                   cluster0_11_maxOUTPUT, cluster0_11_minOUTPUT,
                                   child0_11_maxVal, child0_11_minVal,
                                   inputPCAANNReg0_11, modelANNReg0_11, outputStandardizeANNReg0_11,
                                   listParticles[p]->m_T_gas, listParticles[p]->m_Yk_gas,
                                   numVarANN, nsp,
                                   input_childLocalStandardized, 
                                   outputStandardizeANN_Vec, mixture);
                        break;
                        
                    // Cluster 0_12
                    case 12:
                        if(rank==0)
                        {
                            cout << " Case0_12"<<  endl;
                        }
                        advanceANN_withoutPCA(cluster0_12_meanScaleINPUT, cluster0_12_stdScaleINPUT,
                                   cluster0_12_meanScaleOUTPUT, cluster0_12_stdScaleOUTPUT,
                                   cluster0_12_maxOUTPUT, cluster0_12_minOUTPUT,
                                   child0_12_maxVal, child0_12_minVal,
                                   inputPCAANNReg0_12, modelANNReg0_12, outputStandardizeANNReg0_12,
                                   listParticles[p]->m_T_gas, listParticles[p]->m_Yk_gas,
                                   numVarANN, nsp,
                                   input_childLocalStandardized, 
                                   outputStandardizeANN_Vec, mixture);
                        break;
                        
                        // Cluster 0_13
                    case 13:
                        if(rank==0)
                        {
                            cout << " Case0_13"<<  endl;
                        }
                        advanceANN_withoutPCA(cluster0_13_meanScaleINPUT, cluster0_13_stdScaleINPUT,
                                   cluster0_13_meanScaleOUTPUT, cluster0_13_stdScaleOUTPUT,
                                   cluster0_13_maxOUTPUT, cluster0_13_minOUTPUT,
                                   child0_13_maxVal, child0_13_minVal,
                                   inputPCAANNReg0_13, modelANNReg0_13, outputStandardizeANNReg0_13,
                                   listParticles[p]->m_T_gas, listParticles[p]->m_Yk_gas,
                                   numVarANN, nsp,
                                   input_childLocalStandardized, 
                                   outputStandardizeANN_Vec, mixture);
                        break;
                        
                        // Cluster 0_14
                    case 14:
                        if(rank==0)
                        {
                            cout << " Case0_14"<<  endl;
                        }
                        advanceANN_withoutPCA(cluster0_14_meanScaleINPUT, cluster0_14_stdScaleINPUT,
                                   cluster0_14_meanScaleOUTPUT, cluster0_14_stdScaleOUTPUT,
                                   cluster0_14_maxOUTPUT, cluster0_14_minOUTPUT,
                                   child0_14_maxVal, child0_14_minVal,
                                   inputPCAANNReg0_14, modelANNReg0_14, outputStandardizeANNReg0_14,
                                   listParticles[p]->m_T_gas, listParticles[p]->m_Yk_gas,
                                   numVarANN, nsp,
                                   input_childLocalStandardized, 
                                   outputStandardizeANN_Vec, mixture);
                        break;
                        
                        // Cluster 0_15
                    case 15:
                        if(rank==0)
                        {
                            cout << " Case0_15"<<  endl;
                        }
                        advanceANN_withoutPCA(cluster0_15_meanScaleINPUT, cluster0_15_stdScaleINPUT,
                                   cluster0_15_meanScaleOUTPUT, cluster0_15_stdScaleOUTPUT,
                                   cluster0_15_maxOUTPUT, cluster0_15_minOUTPUT,
                                   child0_15_maxVal, child0_15_minVal,
                                   inputPCAANNReg0_15, modelANNReg0_15, outputStandardizeANNReg0_15,
                                   listParticles[p]->m_T_gas, listParticles[p]->m_Yk_gas,
                                   numVarANN, nsp,
                                   input_childLocalStandardized, 
                                   outputStandardizeANN_Vec, mixture);
                        break;
                        
                        // Cluster 0_16
                    case 16:
                        if(rank==0)
                        {
                            cout << " Case0_16"<<  endl;
                        }
                        advanceANN_withoutPCA(cluster0_16_meanScaleINPUT, cluster0_16_stdScaleINPUT,
                                   cluster0_16_meanScaleOUTPUT, cluster0_16_stdScaleOUTPUT,
                                   cluster0_16_maxOUTPUT, cluster0_16_minOUTPUT,
                                   child0_16_maxVal, child0_16_minVal,
                                   inputPCAANNReg0_16, modelANNReg0_16, outputStandardizeANNReg0_16,
                                   listParticles[p]->m_T_gas, listParticles[p]->m_Yk_gas,
                                   numVarANN, nsp,
                                   input_childLocalStandardized, 
                                   outputStandardizeANN_Vec, mixture);
                        break;
                        
                        // Cluster 0_17
                    case 17:
                        if(rank==0)
                        {
                            cout << " Case0_17"<<  endl;
                        }
                        advanceANN_withoutPCA(cluster0_17_meanScaleINPUT, cluster0_17_stdScaleINPUT,
                                   cluster0_17_meanScaleOUTPUT, cluster0_17_stdScaleOUTPUT,
                                   cluster0_17_maxOUTPUT, cluster0_17_minOUTPUT,
                                   child0_17_maxVal, child0_17_minVal,
                                   inputPCAANNReg0_17, modelANNReg0_17, outputStandardizeANNReg0_17,
                                   listParticles[p]->m_T_gas, listParticles[p]->m_Yk_gas,
                                   numVarANN, nsp,
                                   input_childLocalStandardized, 
                                   outputStandardizeANN_Vec, mixture);
                        break;
                  
                        // Cluster 0_18
                    case 18:
                        if(rank==0)
                        {
                            cout << " Case0_18"<<  endl;
                        }
                        advanceANN_withoutPCA(cluster0_18_meanScaleINPUT, cluster0_18_stdScaleINPUT,
                                   cluster0_18_meanScaleOUTPUT, cluster0_18_stdScaleOUTPUT,
                                   cluster0_18_maxOUTPUT, cluster0_18_minOUTPUT,
                                   child0_18_maxVal, child0_18_minVal,
                                   inputPCAANNReg0_18, modelANNReg0_18, outputStandardizeANNReg0_18,
                                   listParticles[p]->m_T_gas, listParticles[p]->m_Yk_gas,
                                   numVarANN, nsp,
                                   input_childLocalStandardized, 
                                   outputStandardizeANN_Vec, mixture);
                        break;
                  
                        // Cluster 0_19
                    case 19:
                        if(rank==0)
                        {
                            cout << " Case0_19"<<  endl;
                        }
                        advanceANN_withoutPCA(cluster0_19_meanScaleINPUT, cluster0_19_stdScaleINPUT,
                                   cluster0_19_meanScaleOUTPUT, cluster0_19_stdScaleOUTPUT,
                                   cluster0_19_maxOUTPUT, cluster0_19_minOUTPUT,
                                   child0_19_maxVal, child0_19_minVal,
                                   inputPCAANNReg0_19, modelANNReg0_19, outputStandardizeANNReg0_19,
                                   listParticles[p]->m_T_gas, listParticles[p]->m_Yk_gas,
                                   numVarANN, nsp,
                                   input_childLocalStandardized, 
                                   outputStandardizeANN_Vec, mixture);
                        break;
                        
                } //END switch(indexChildCluster) parentCluster0
                
                // Release OPENCV MAT
                dParent.release();
            } // END parentCluster0
            
            else //parentCluster1
            {
                // LOCAL PARENT Standardize Y ,T for Kmeans LOCAL PARENT: standardizedX = (X - meanScale) / standardDeviationScale
                // Standardize T
                input_parentLocalStandardized[numVarANN-1] = listParticles[p]->m_T_gas - parentCluster1_meanScaleINPUT[numVarANN-1];
                input_parentLocalStandardized[numVarANN-1] /= parentCluster1_stdScaleINPUT[numVarANN-1];
                // Standardize Y
                for (int kANN=0; kANN<nsp-1; kANN++) // Don't take N2 (last species)
                {
                    input_parentLocalStandardized[kANN] = listParticles[p]->m_Yk_gas[kANN] - parentCluster1_meanScaleINPUT[kANN];
                    input_parentLocalStandardized[kANN] /= parentCluster1_stdScaleINPUT[kANN];
                    //Cross-check if there is N2 inside, if yes, should move N2 to the end of the scheme
                    if(mixture->speciesName(kANN) == "N2") {cout << "Warning : N2 in the ANN!!!" << endl;}
                }
                
                // CLASSIFIER LOCAL PARENT: Caculate the distance : Eucledian (Kmeans) or Projection parition (LPCA)
                // >> Determine indexChildCluster
                // LPCA example: cv::Mat d0 = (x-r0).t() * eigenVec0.t() * eigenVec0 * (x-r0);
                // Eucledean: Can be calculated by std::vector
                
                // Declare vector x
                xParent = cv::Mat(numVarANN, 1, CV_32F, input_parentLocalStandardized.data());
                d1_0 = cv::norm(xParent,r1_0,cv::NORM_L2SQR); // Calculate Euclidean distance from x to cluster 0
                d1_1 = cv::norm(xParent,r1_1,cv::NORM_L2SQR);
                d1_2 = cv::norm(xParent,r1_2,cv::NORM_L2SQR);


                // Group distance to d
                cv::Mat dParent;
                dParent.push_back(d1_0);
                dParent.push_back(d1_1);
                dParent.push_back(d1_2);
                
                // Take position of mean (argmin) = indexParentCluster
                cv::minMaxLoc(dParent, NULL, NULL, &minLocParent, NULL); // d, minVal, maxVal, minLoc, maxLoc
                indexChildCluster = minLocParent.y; //Point has 2 coordinates (x,y)
                
                
                if(indexChildCluster==0) // childCluster 1_0
                {
                    // Classifier grand child
                    // LOCAL CHILD Standardize Y ,T for Kmeans LOCAL CHILD: standardizedX = (X - meanScale) / standardDeviationScale
                    // Standardize T
                    input_childLocalStandardized[numVarANN-1] = listParticles[p]->m_T_gas - cluster1_0_meanScaleINPUT[numVarANN-1];
                    input_childLocalStandardized[numVarANN-1] /= cluster1_0_stdScaleINPUT[numVarANN-1];
                    // Standardize Y
                    for (int kANN=0; kANN<nsp-1; kANN++) // Don't take N2 (last species)
                    {
                        input_childLocalStandardized[kANN] = listParticles[p]->m_Yk_gas[kANN] - cluster1_0_meanScaleINPUT[kANN];
                        input_childLocalStandardized[kANN] /= cluster1_0_stdScaleINPUT[kANN];
                        //Cross-check if there is N2 inside, if yes, should move N2 to the end of the scheme
                        if(mixture->speciesName(kANN) == "N2") {cout << "Warning : N2 in the ANN!!!" << endl;}
                    }

                    
                    // CLASSIFIER LOCAL CHILD: Caculate the distance : Eucledian (Kmeans) or Projection parition (LPCA)
                    // >> Determine indexGrandChildCluster
                    // LPCA example: cv::Mat d0 = (x-r0).t() * eigenVec0.t() * eigenVec0 * (x-r0);
                    // Eucledean: Can be calculated by std::vector
                    
                    // Declare vector x
                    xChild = cv::Mat(numVarANN, 1, CV_32F, input_childLocalStandardized.data());
                    d1_0_0 = cv::norm(xChild,r1_0_0,cv::NORM_L2SQR); // Calculate Euclidean distance from x to cluster 0
                    d1_0_1 = cv::norm(xChild,r1_0_1,cv::NORM_L2SQR);
                    d1_0_2 = cv::norm(xChild,r1_0_2,cv::NORM_L2SQR);
                    d1_0_3 = cv::norm(xChild,r1_0_3,cv::NORM_L2SQR);
                    d1_0_4 = cv::norm(xChild,r1_0_4,cv::NORM_L2SQR);
                    d1_0_5 = cv::norm(xChild,r1_0_5,cv::NORM_L2SQR);
                    d1_0_6 = cv::norm(xChild,r1_0_6,cv::NORM_L2SQR);
                    d1_0_7 = cv::norm(xChild,r1_0_7,cv::NORM_L2SQR);
                    d1_0_8 = cv::norm(xChild,r1_0_8,cv::NORM_L2SQR);
                    d1_0_9 = cv::norm(xChild,r1_0_9,cv::NORM_L2SQR);
                    d1_0_10 = cv::norm(xChild,r1_0_10,cv::NORM_L2SQR);
                    d1_0_11 = cv::norm(xChild,r1_0_11,cv::NORM_L2SQR);
                    d1_0_12 = cv::norm(xChild,r1_0_12,cv::NORM_L2SQR);
                    d1_0_13 = cv::norm(xChild,r1_0_13,cv::NORM_L2SQR);
                    d1_0_14 = cv::norm(xChild,r1_0_14,cv::NORM_L2SQR);
                    d1_0_15 = cv::norm(xChild,r1_0_15,cv::NORM_L2SQR);
                    d1_0_16 = cv::norm(xChild,r1_0_16,cv::NORM_L2SQR);
                    d1_0_17 = cv::norm(xChild,r1_0_17,cv::NORM_L2SQR);
                    d1_0_18 = cv::norm(xChild,r1_0_18,cv::NORM_L2SQR);
                    d1_0_19 = cv::norm(xChild,r1_0_19,cv::NORM_L2SQR);

                
                    // Group distance to d
                    cv::Mat dChild;
                    dChild.push_back(d1_0_0);
                    dChild.push_back(d1_0_1);
                    dChild.push_back(d1_0_2);
                    dChild.push_back(d1_0_3);
                    dChild.push_back(d1_0_4);
                    dChild.push_back(d1_0_5);
                    dChild.push_back(d1_0_6);
                    dChild.push_back(d1_0_7);
                    dChild.push_back(d1_0_8);
                    dChild.push_back(d1_0_9);
                    dChild.push_back(d1_0_10);
                    dChild.push_back(d1_0_11);
                    dChild.push_back(d1_0_12);
                    dChild.push_back(d1_0_13);
                    dChild.push_back(d1_0_14);
                    dChild.push_back(d1_0_15);
                    dChild.push_back(d1_0_16);
                    dChild.push_back(d1_0_17);
                    dChild.push_back(d1_0_18);
                    dChild.push_back(d1_0_19);

                    // Take position of mean (argmin) = indexParentCluster
                    cv::minMaxLoc(dChild, NULL, NULL, &minLocChild, NULL); // d, minVal, maxVal, minLoc, maxLoc
                    indexGrandChildCluster = minLocChild.y; //Point has 2 coordinates (x,y)
                    
                    if(rank==0) // Print to check
                    {
                        cout << "Particles " << p << " || parentCluster = " << indexParentCluster;
                        cout << " || childCluster = " << indexChildCluster << " || grandChildCluster = " << indexGrandChildCluster << " ";
                    }
                    
                    // Save indexCluster of this particle
                    if(activateIndexCluster)
                    {
                        if(rank==0)
                        {
                            ofstream indexClusterFile("outputs/indexCluster.dat",ios::app); //ios::app = append at the end of the file
                            indexClusterFile << t << "	";
                            indexClusterFile << p << "	";
                            indexClusterFile << indexParentCluster <<  "	";
                            indexClusterFile << indexChildCluster <<  "	";
                            indexClusterFile << indexGrandChildCluster << endl;
                            indexClusterFile.close();
                        }
                    }
                    
                    switch(indexGrandChildCluster) // cluster1_0_X
                    {
                        // Cluster 1_0_0
                        case 0:
                            if(rank==0)
                            {
                                cout << " Case1_0_0"<<  endl;
                            }
                            advanceANN_withoutPCA(cluster1_0_0_meanScaleINPUT, cluster1_0_0_stdScaleINPUT,
                                       cluster1_0_0_meanScaleOUTPUT, cluster1_0_0_stdScaleOUTPUT,
                                       cluster1_0_0_maxOUTPUT, cluster1_0_0_minOUTPUT,
                                       grandChild1_0_0_maxVal, grandChild1_0_0_minVal, // grandChild
                                       inputPCAANNReg1_0_0, modelANNReg1_0_0, outputStandardizeANNReg1_0_0,
                                       listParticles[p]->m_T_gas, listParticles[p]->m_Yk_gas,
                                       numVarANN, nsp,
                                       input_childLocalStandardized, 
                                       outputStandardizeANN_Vec, mixture);
                            break;
                            
                        // Cluster 1_0_1
                        case 1:
                            if(rank==0)
                            {
                                cout << " Case1_0_1"<<  endl;
                            }
                            advanceANN_withoutPCA(cluster1_0_1_meanScaleINPUT, cluster1_0_1_stdScaleINPUT,
                                       cluster1_0_1_meanScaleOUTPUT, cluster1_0_1_stdScaleOUTPUT,
                                       cluster1_0_1_maxOUTPUT, cluster1_0_1_minOUTPUT,
                                       grandChild1_0_1_maxVal, grandChild1_0_1_minVal, // grandChild
                                       inputPCAANNReg1_0_1, modelANNReg1_0_1, outputStandardizeANNReg1_0_1,
                                       listParticles[p]->m_T_gas, listParticles[p]->m_Yk_gas,
                                       numVarANN, nsp,
                                       input_childLocalStandardized, 
                                       outputStandardizeANN_Vec, mixture);
                            break;
                            
                        // Cluster 1_0_2
                        case 2:
                            if(rank==0)
                            {
                                cout << " Case1_0_2"<<  endl;
                            }
                            advanceANN_withoutPCA(cluster1_0_2_meanScaleINPUT, cluster1_0_2_stdScaleINPUT,
                                       cluster1_0_2_meanScaleOUTPUT, cluster1_0_2_stdScaleOUTPUT,
                                       cluster1_0_2_maxOUTPUT, cluster1_0_2_minOUTPUT,
                                       grandChild1_0_2_maxVal, grandChild1_0_2_minVal, // grandChild
                                       inputPCAANNReg1_0_2, modelANNReg1_0_2, outputStandardizeANNReg1_0_2,
                                       listParticles[p]->m_T_gas, listParticles[p]->m_Yk_gas,
                                       numVarANN, nsp,
                                       input_childLocalStandardized, 
                                       outputStandardizeANN_Vec, mixture);
                            break;
                            
                        // Cluster 1_0_3
                        case 3:
                            if(rank==0)
                            {
                                cout << " Case1_0_3"<<  endl;
                            }
                            advanceANN_withoutPCA(cluster1_0_3_meanScaleINPUT, cluster1_0_3_stdScaleINPUT,
                                       cluster1_0_3_meanScaleOUTPUT, cluster1_0_3_stdScaleOUTPUT,
                                       cluster1_0_3_maxOUTPUT, cluster1_0_3_minOUTPUT,
                                       grandChild1_0_3_maxVal, grandChild1_0_3_minVal, // grandChild
                                       inputPCAANNReg1_0_3, modelANNReg1_0_3, outputStandardizeANNReg1_0_3,
                                       listParticles[p]->m_T_gas, listParticles[p]->m_Yk_gas,
                                       numVarANN, nsp,
                                       input_childLocalStandardized, 
                                       outputStandardizeANN_Vec, mixture);
                            break;
                            
                        // Cluster 1_0_4
                        case 4:
                            if(rank==0)
                            {
                                cout << " Case1_0_4"<<  endl;
                            }
                            advanceANN_withoutPCA(cluster1_0_4_meanScaleINPUT, cluster1_0_4_stdScaleINPUT,
                                       cluster1_0_4_meanScaleOUTPUT, cluster1_0_4_stdScaleOUTPUT,
                                       cluster1_0_4_maxOUTPUT, cluster1_0_4_minOUTPUT,
                                       grandChild1_0_4_maxVal, grandChild1_0_4_minVal, // grandChild
                                       inputPCAANNReg1_0_4, modelANNReg1_0_4, outputStandardizeANNReg1_0_4,
                                       listParticles[p]->m_T_gas, listParticles[p]->m_Yk_gas,
                                       numVarANN, nsp,
                                       input_childLocalStandardized, 
                                       outputStandardizeANN_Vec, mixture);
                            break;
                            
                            // Cluster 1_0_5
                        case 5:
                            if(rank==0)
                            {
                                cout << " Case1_0_5"<<  endl;
                            }
                            advanceANN_withoutPCA(cluster1_0_5_meanScaleINPUT, cluster1_0_5_stdScaleINPUT,
                                       cluster1_0_5_meanScaleOUTPUT, cluster1_0_5_stdScaleOUTPUT,
                                       cluster1_0_5_maxOUTPUT, cluster1_0_5_minOUTPUT,
                                       grandChild1_0_5_maxVal, grandChild1_0_5_minVal, // grandChild
                                       inputPCAANNReg1_0_5, modelANNReg1_0_5, outputStandardizeANNReg1_0_5,
                                       listParticles[p]->m_T_gas, listParticles[p]->m_Yk_gas,
                                       numVarANN, nsp,
                                       input_childLocalStandardized, 
                                       outputStandardizeANN_Vec, mixture);
                            break;
                            
                            // Cluster 1_0_6
                        case 6:
                            if(rank==0)
                            {
                                cout << " Case1_0_6"<<  endl;
                            }
                            advanceANN_withoutPCA(cluster1_0_6_meanScaleINPUT, cluster1_0_6_stdScaleINPUT,
                                       cluster1_0_6_meanScaleOUTPUT, cluster1_0_6_stdScaleOUTPUT,
                                       cluster1_0_6_maxOUTPUT, cluster1_0_6_minOUTPUT,
                                       grandChild1_0_6_maxVal, grandChild1_0_6_minVal, // grandChild,
                                       inputPCAANNReg1_0_6, modelANNReg1_0_6, outputStandardizeANNReg1_0_6,
                                       listParticles[p]->m_T_gas, listParticles[p]->m_Yk_gas,
                                       numVarANN, nsp,
                                       input_childLocalStandardized, 
                                       outputStandardizeANN_Vec, mixture);
                            break;
                            
                            // Cluster 1_0_7
                        case 7:
                            if(rank==0)
                            {
                                cout << " Case1_0_7"<<  endl;
                            }
                            advanceANN_withoutPCA(cluster1_0_7_meanScaleINPUT, cluster1_0_7_stdScaleINPUT,
                                       cluster1_0_7_meanScaleOUTPUT, cluster1_0_7_stdScaleOUTPUT,
                                       cluster1_0_7_maxOUTPUT, cluster1_0_7_minOUTPUT,
                                       grandChild1_0_7_maxVal, grandChild1_0_7_minVal, // grandChild
                                       inputPCAANNReg1_0_7, modelANNReg1_0_7, outputStandardizeANNReg1_0_7,
                                       listParticles[p]->m_T_gas, listParticles[p]->m_Yk_gas,
                                       numVarANN, nsp,
                                       input_childLocalStandardized, 
                                       outputStandardizeANN_Vec, mixture);
                            break;
                            
                            // Cluster 1_0_8
                        case 8:
                            if(rank==0)
                            {
                                cout << " Case1_0_8"<<  endl;
                            }
                            advanceANN_withoutPCA(cluster1_0_8_meanScaleINPUT, cluster1_0_8_stdScaleINPUT,
                                       cluster1_0_8_meanScaleOUTPUT, cluster1_0_8_stdScaleOUTPUT,
                                       cluster1_0_8_maxOUTPUT, cluster1_0_8_minOUTPUT,
                                       grandChild1_0_8_maxVal, grandChild1_0_8_minVal, // grandChild
                                       inputPCAANNReg1_0_8, modelANNReg1_0_8, outputStandardizeANNReg1_0_8,
                                       listParticles[p]->m_T_gas, listParticles[p]->m_Yk_gas,
                                       numVarANN, nsp,
                                       input_childLocalStandardized, 
                                       outputStandardizeANN_Vec, mixture);
                            break;
                            
                            // Cluster 1_0_9
                        case 9:
                            if(rank==0)
                            {
                                cout << " Case1_0_9"<<  endl;
                            }
                            advanceANN_withoutPCA(cluster1_0_9_meanScaleINPUT, cluster1_0_9_stdScaleINPUT,
                                       cluster1_0_9_meanScaleOUTPUT, cluster1_0_9_stdScaleOUTPUT,
                                       cluster1_0_9_maxOUTPUT, cluster1_0_9_minOUTPUT,
                                       grandChild1_0_9_maxVal, grandChild1_0_9_minVal, // grandChild
                                       inputPCAANNReg1_0_9, modelANNReg1_0_9, outputStandardizeANNReg1_0_9,
                                       listParticles[p]->m_T_gas, listParticles[p]->m_Yk_gas,
                                       numVarANN, nsp,
                                       input_childLocalStandardized, 
                                       outputStandardizeANN_Vec, mixture);
                            break;
                            
                            // Cluster 1_0_10
                        case 10:
                            if(rank==0)
                            {
                                cout << " Case1_0_10"<<  endl;
                            }
                            advanceANN_withoutPCA(cluster1_0_10_meanScaleINPUT, cluster1_0_10_stdScaleINPUT,
                                       cluster1_0_10_meanScaleOUTPUT, cluster1_0_10_stdScaleOUTPUT,
                                       cluster1_0_10_maxOUTPUT, cluster1_0_10_minOUTPUT,
                                       grandChild1_0_10_maxVal, grandChild1_0_10_minVal, // grandChild
                                       inputPCAANNReg1_0_10, modelANNReg1_0_10, outputStandardizeANNReg1_0_10,
                                       listParticles[p]->m_T_gas, listParticles[p]->m_Yk_gas,
                                       numVarANN, nsp,
                                       input_childLocalStandardized, 
                                       outputStandardizeANN_Vec, mixture);
                            break;
                            
                            // Cluster 1_0_11
                        case 11:
                            if(rank==0)
                            {
                                cout << " Case1_0_11"<<  endl;
                            }
                            advanceANN_withoutPCA(cluster1_0_11_meanScaleINPUT, cluster1_0_11_stdScaleINPUT,
                                       cluster1_0_11_meanScaleOUTPUT, cluster1_0_11_stdScaleOUTPUT,
                                       cluster1_0_11_maxOUTPUT, cluster1_0_11_minOUTPUT,
                                       grandChild1_0_11_maxVal, grandChild1_0_11_minVal, // grandChild
                                       inputPCAANNReg1_0_11, modelANNReg1_0_11, outputStandardizeANNReg1_0_11,
                                       listParticles[p]->m_T_gas, listParticles[p]->m_Yk_gas,
                                       numVarANN, nsp,
                                       input_childLocalStandardized, 
                                       outputStandardizeANN_Vec, mixture);
                            break;
                            
                            // Cluster 1_0_12
                        case 12:
                            if(rank==0)
                            {
                                cout << " Case1_0_12"<<  endl;
                            }
                            advanceANN_withoutPCA(cluster1_0_12_meanScaleINPUT, cluster1_0_12_stdScaleINPUT,
                                       cluster1_0_12_meanScaleOUTPUT, cluster1_0_12_stdScaleOUTPUT,
                                       cluster1_0_12_maxOUTPUT, cluster1_0_12_minOUTPUT,
                                       grandChild1_0_12_maxVal, grandChild1_0_12_minVal, // grandChild
                                       inputPCAANNReg1_0_12, modelANNReg1_0_12, outputStandardizeANNReg1_0_12,
                                       listParticles[p]->m_T_gas, listParticles[p]->m_Yk_gas,
                                       numVarANN, nsp,
                                       input_childLocalStandardized, 
                                       outputStandardizeANN_Vec, mixture);
                            break;
                            
                            // Cluster 1_0_13
                        case 13:
                            if(rank==0)
                            {
                                cout << " Case1_0_13"<<  endl;
                            }
                            advanceANN_withoutPCA(cluster1_0_13_meanScaleINPUT, cluster1_0_13_stdScaleINPUT,
                                       cluster1_0_13_meanScaleOUTPUT, cluster1_0_13_stdScaleOUTPUT,
                                       cluster1_0_13_maxOUTPUT, cluster1_0_13_minOUTPUT,
                                       grandChild1_0_13_maxVal, grandChild1_0_13_minVal, // grandChild
                                       inputPCAANNReg1_0_13, modelANNReg1_0_13, outputStandardizeANNReg1_0_13,
                                       listParticles[p]->m_T_gas, listParticles[p]->m_Yk_gas,
                                       numVarANN, nsp,
                                       input_childLocalStandardized, 
                                       outputStandardizeANN_Vec, mixture);
                            break;
                            
                            // Cluster 1_0_14
                        case 14:
                            if(rank==0)
                            {
                                cout << " Case1_0_14"<<  endl;
                            }
                            advanceANN_withoutPCA(cluster1_0_14_meanScaleINPUT, cluster1_0_14_stdScaleINPUT,
                                       cluster1_0_14_meanScaleOUTPUT, cluster1_0_14_stdScaleOUTPUT,
                                       cluster1_0_14_maxOUTPUT, cluster1_0_14_minOUTPUT,
                                       grandChild1_0_14_maxVal, grandChild1_0_14_minVal, // grandChild
                                       inputPCAANNReg1_0_14, modelANNReg1_0_14, outputStandardizeANNReg1_0_14,
                                       listParticles[p]->m_T_gas, listParticles[p]->m_Yk_gas,
                                       numVarANN, nsp,
                                       input_childLocalStandardized, 
                                       outputStandardizeANN_Vec, mixture);
                            break;
                            
                            // Cluster 1_0_15
                        case 15:
                            if(rank==0)
                            {
                                cout << " Case1_0_15"<<  endl;
                            }
                            advanceANN_withoutPCA(cluster1_0_15_meanScaleINPUT, cluster1_0_15_stdScaleINPUT,
                                       cluster1_0_15_meanScaleOUTPUT, cluster1_0_15_stdScaleOUTPUT,
                                       cluster1_0_15_maxOUTPUT, cluster1_0_15_minOUTPUT,
                                       grandChild1_0_15_maxVal, grandChild1_0_15_minVal, // grandChild
                                       inputPCAANNReg1_0_15, modelANNReg1_0_15, outputStandardizeANNReg1_0_15,
                                       listParticles[p]->m_T_gas, listParticles[p]->m_Yk_gas,
                                       numVarANN, nsp,
                                       input_childLocalStandardized, 
                                       outputStandardizeANN_Vec, mixture);
                            break;
                            
                            // Cluster 1_0_16
                        case 16:
                            if(rank==0)
                            {
                                cout << " Case1_0_16"<<  endl;
                            }
                            advanceANN_withoutPCA(cluster1_0_16_meanScaleINPUT, cluster1_0_16_stdScaleINPUT,
                                       cluster1_0_16_meanScaleOUTPUT, cluster1_0_16_stdScaleOUTPUT,
                                       cluster1_0_16_maxOUTPUT, cluster1_0_16_minOUTPUT,
                                       grandChild1_0_16_maxVal, grandChild1_0_16_minVal, // grandChild
                                       inputPCAANNReg1_0_16, modelANNReg1_0_16, outputStandardizeANNReg1_0_16,
                                       listParticles[p]->m_T_gas, listParticles[p]->m_Yk_gas,
                                       numVarANN, nsp,
                                       input_childLocalStandardized, 
                                       outputStandardizeANN_Vec, mixture);
                            break;
                            
                            // Cluster 1_0_17
                        case 17:
                            if(rank==0)
                            {
                                cout << " Case1_0_17"<<  endl;
                            }
                            advanceANN_withoutPCA(cluster1_0_17_meanScaleINPUT, cluster1_0_17_stdScaleINPUT,
                                       cluster1_0_17_meanScaleOUTPUT, cluster1_0_17_stdScaleOUTPUT,
                                       cluster1_0_17_maxOUTPUT, cluster1_0_17_minOUTPUT,
                                       grandChild1_0_17_maxVal, grandChild1_0_17_minVal, // grandChild
                                       inputPCAANNReg1_0_17, modelANNReg1_0_17, outputStandardizeANNReg1_0_17,
                                       listParticles[p]->m_T_gas, listParticles[p]->m_Yk_gas,
                                       numVarANN, nsp,
                                       input_childLocalStandardized, 
                                       outputStandardizeANN_Vec, mixture);
                            break;
                            
                            // Cluster 1_0_18
                        case 18:
                            if(rank==0)
                            {
                                cout << " Case1_0_18"<<  endl;
                            }
                            advanceANN_withoutPCA(cluster1_0_18_meanScaleINPUT, cluster1_0_18_stdScaleINPUT,
                                       cluster1_0_18_meanScaleOUTPUT, cluster1_0_18_stdScaleOUTPUT,
                                       cluster1_0_18_maxOUTPUT, cluster1_0_18_minOUTPUT,
                                       grandChild1_0_18_maxVal, grandChild1_0_18_minVal, // grandChild
                                       inputPCAANNReg1_0_18, modelANNReg1_0_18, outputStandardizeANNReg1_0_18,
                                       listParticles[p]->m_T_gas, listParticles[p]->m_Yk_gas,
                                       numVarANN, nsp,
                                       input_childLocalStandardized, 
                                       outputStandardizeANN_Vec, mixture);
                            break;
                            
                            // Cluster 1_0_19
                        case 19:
                            if(rank==0)
                            {
                                cout << " Case1_0_19"<<  endl;
                            }
                            advanceANN_withoutPCA(cluster1_0_19_meanScaleINPUT, cluster1_0_19_stdScaleINPUT,
                                       cluster1_0_19_meanScaleOUTPUT, cluster1_0_19_stdScaleOUTPUT,
                                       cluster1_0_19_maxOUTPUT, cluster1_0_19_minOUTPUT,
                                       grandChild1_0_19_maxVal, grandChild1_0_19_minVal, // grandChild
                                       inputPCAANNReg1_0_19, modelANNReg1_0_19, outputStandardizeANNReg1_0_19,
                                       listParticles[p]->m_T_gas, listParticles[p]->m_Yk_gas,
                                       numVarANN, nsp,
                                       input_childLocalStandardized, 
                                       outputStandardizeANN_Vec, mixture);
                            break;
                            
                            
                    } // END switch(indexGrandChildCluster) // cluster1_0_X
                    
                } // END if(indexChildCluster==0)
                else if(indexChildCluster==1) // childCluster1_1
                {
                    // Classifier grand child
                    // LOCAL CHILD Standardize Y ,T for Kmeans LOCAL CHILD: standardizedX = (X - meanScale) / standardDeviationScale
                    // Standardize T
                    input_childLocalStandardized[numVarANN-1] = listParticles[p]->m_T_gas - cluster1_1_meanScaleINPUT[numVarANN-1];
                    input_childLocalStandardized[numVarANN-1] /= cluster1_1_stdScaleINPUT[numVarANN-1];
                    // Standardize Y
                    for (int kANN=0; kANN<nsp-1; kANN++) // Don't take N2 (last species)
                    {
                        input_childLocalStandardized[kANN] = listParticles[p]->m_Yk_gas[kANN] - cluster1_1_meanScaleINPUT[kANN];
                        input_childLocalStandardized[kANN] /= cluster1_1_stdScaleINPUT[kANN];
                        //Cross-check if there is N2 inside, if yes, should move N2 to the end of the scheme
                        if(mixture->speciesName(kANN) == "N2") {cout << "Warning : N2 in the ANN!!!" << endl;}
                    }
                    
                    
                    // CLASSIFIER LOCAL CHILD: Caculate the distance : Eucledian (Kmeans) or Projection parition (LPCA)
                    // >> Determine indexGrandChildCluster
                    // LPCA example: cv::Mat d0 = (x-r0).t() * eigenVec0.t() * eigenVec0 * (x-r0);
                    // Eucledean: Can be calculated by std::vector
                    
                    // Declare vector x
                    xChild = cv::Mat(numVarANN, 1, CV_32F, input_childLocalStandardized.data());
                    d1_1_0 = cv::norm(xChild,r1_1_0,cv::NORM_L2SQR); // Calculate Euclidean distance from x to cluster 0
                    d1_1_1 = cv::norm(xChild,r1_1_1,cv::NORM_L2SQR);
                    d1_1_2 = cv::norm(xChild,r1_1_2,cv::NORM_L2SQR);
                    d1_1_3 = cv::norm(xChild,r1_1_3,cv::NORM_L2SQR);
                    d1_1_4 = cv::norm(xChild,r1_1_4,cv::NORM_L2SQR);
                    d1_1_5 = cv::norm(xChild,r1_1_5,cv::NORM_L2SQR);
                    d1_1_6 = cv::norm(xChild,r1_1_6,cv::NORM_L2SQR);
                    d1_1_7 = cv::norm(xChild,r1_1_7,cv::NORM_L2SQR);
                    d1_1_8 = cv::norm(xChild,r1_1_8,cv::NORM_L2SQR);
                    d1_1_9 = cv::norm(xChild,r1_1_9,cv::NORM_L2SQR);


                    // Group distance to d
                    cv::Mat dChild;
                    dChild.push_back(d1_1_0);
                    dChild.push_back(d1_1_1);
                    dChild.push_back(d1_1_2);
                    dChild.push_back(d1_1_3);
                    dChild.push_back(d1_1_4);
                    dChild.push_back(d1_1_5);
                    dChild.push_back(d1_1_6);
                    dChild.push_back(d1_1_7);
                    dChild.push_back(d1_1_8);
                    dChild.push_back(d1_1_9);

                    
                    // Take position of mean (argmin) = indexParentCluster
                    cv::minMaxLoc(dChild, NULL, NULL, &minLocChild, NULL); // d, minVal, maxVal, minLoc, maxLoc
                    indexGrandChildCluster = minLocChild.y; //Point has 2 coordinates (x,y)
                    
                    if(rank==0) // Print to check
                    {
                        cout << "Particles " << p << " || parentCluster = " << indexParentCluster;
                        cout << " || childCluster = " << indexChildCluster << " || grandChildCluster = " << indexGrandChildCluster << " ";
                    }
                    
                    // Save indexCluster of this particle
                    if(activateIndexCluster)
                    {
                        if(rank==0)
                        {
                            ofstream indexClusterFile("outputs/indexCluster.dat",ios::app); //ios::app = append at the end of the file
                            indexClusterFile << t << "	";
                            indexClusterFile << p << "	";
                            indexClusterFile << indexParentCluster <<  "	";
                            indexClusterFile << indexChildCluster <<  "	";
                            indexClusterFile << indexGrandChildCluster << endl;
                            indexClusterFile.close();
                        }
                    }
                    
                    switch(indexGrandChildCluster) // cluster1_1_X
                    {
                            
                        // Cluster 1_1_0
                        case 0:
                            if(rank==0)
                            {
                                cout << " Case1_1_0"<<  endl;
                            }
                            advanceANN_withoutPCA(cluster1_1_0_meanScaleINPUT, cluster1_1_0_stdScaleINPUT,
                                       cluster1_1_0_meanScaleOUTPUT, cluster1_1_0_stdScaleOUTPUT,
                                       cluster1_1_0_maxOUTPUT, cluster1_1_0_minOUTPUT,
                                       grandChild1_1_0_maxVal, grandChild1_1_0_minVal, // grandChild
                                       inputPCAANNReg1_1_0, modelANNReg1_1_0, outputStandardizeANNReg1_1_0,
                                       listParticles[p]->m_T_gas, listParticles[p]->m_Yk_gas,
                                       numVarANN, nsp,
                                       input_childLocalStandardized, 
                                       outputStandardizeANN_Vec, mixture);
                            break;
                            
                        // Cluster 1_1_1
                        case 1:
                            if(rank==0)
                            {
                                cout << " Case1_1_1"<<  endl;
                            }
                            advanceANN_withoutPCA(cluster1_1_1_meanScaleINPUT, cluster1_1_1_stdScaleINPUT,
                                       cluster1_1_1_meanScaleOUTPUT, cluster1_1_1_stdScaleOUTPUT,
                                       cluster1_1_1_maxOUTPUT, cluster1_1_1_minOUTPUT,
                                       grandChild1_1_1_maxVal, grandChild1_1_1_minVal, // grandChild
                                       inputPCAANNReg1_1_1, modelANNReg1_1_1, outputStandardizeANNReg1_1_1,
                                       listParticles[p]->m_T_gas, listParticles[p]->m_Yk_gas,
                                       numVarANN, nsp,
                                       input_childLocalStandardized, 
                                       outputStandardizeANN_Vec, mixture);
                            break;
                            
                        // Cluster 1_1_2
                        case 2:
                            if(rank==0)
                            {
                                cout << " Case1_1_2"<<  endl;
                            }
                            advanceANN_withoutPCA(cluster1_1_2_meanScaleINPUT, cluster1_1_2_stdScaleINPUT,
                                       cluster1_1_2_meanScaleOUTPUT, cluster1_1_2_stdScaleOUTPUT,
                                       cluster1_1_2_maxOUTPUT, cluster1_1_2_minOUTPUT,
                                       grandChild1_1_2_maxVal, grandChild1_1_2_minVal, // grandChild
                                       inputPCAANNReg1_1_2, modelANNReg1_1_2, outputStandardizeANNReg1_1_2,
                                       listParticles[p]->m_T_gas, listParticles[p]->m_Yk_gas,
                                       numVarANN, nsp,
                                       input_childLocalStandardized, 
                                       outputStandardizeANN_Vec, mixture);
                            break;
                            
                            // Cluster 1_1_3
                        case 3:
                            if(rank==0)
                            {
                                cout << " Case1_1_3"<<  endl;
                            }
                            advanceANN_withoutPCA(cluster1_1_3_meanScaleINPUT, cluster1_1_3_stdScaleINPUT,
                                       cluster1_1_3_meanScaleOUTPUT, cluster1_1_3_stdScaleOUTPUT,
                                       cluster1_1_3_maxOUTPUT, cluster1_1_3_minOUTPUT,
                                       grandChild1_1_3_maxVal, grandChild1_1_3_minVal, // grandChild
                                       inputPCAANNReg1_1_3, modelANNReg1_1_3, outputStandardizeANNReg1_1_3,
                                       listParticles[p]->m_T_gas, listParticles[p]->m_Yk_gas,
                                       numVarANN, nsp,
                                       input_childLocalStandardized, 
                                       outputStandardizeANN_Vec, mixture);
                            break;
                            
                        // Cluster 1_1_4
                        case 4:
                            if(rank==0)
                            {
                                cout << " Case1_1_4"<<  endl;
                            }
                            advanceANN_withoutPCA(cluster1_1_4_meanScaleINPUT, cluster1_1_4_stdScaleINPUT,
                                       cluster1_1_4_meanScaleOUTPUT, cluster1_1_4_stdScaleOUTPUT,
                                       cluster1_1_4_maxOUTPUT, cluster1_1_4_minOUTPUT,
                                       grandChild1_1_4_maxVal, grandChild1_1_4_minVal, // grandChild
                                       inputPCAANNReg1_1_4, modelANNReg1_1_4, outputStandardizeANNReg1_1_4,
                                       listParticles[p]->m_T_gas, listParticles[p]->m_Yk_gas,
                                       numVarANN, nsp,
                                       input_childLocalStandardized, 
                                       outputStandardizeANN_Vec, mixture);
                            break;
                            
                            // Cluster 1_1_5
                        case 5:
                            if(rank==0)
                            {
                                cout << " Case1_1_5"<<  endl;
                            }
                            advanceANN_withoutPCA(cluster1_1_5_meanScaleINPUT, cluster1_1_5_stdScaleINPUT,
                                       cluster1_1_5_meanScaleOUTPUT, cluster1_1_5_stdScaleOUTPUT,
                                       cluster1_1_5_maxOUTPUT, cluster1_1_5_minOUTPUT,
                                       grandChild1_1_5_maxVal, grandChild1_1_5_minVal, // grandChild
                                       inputPCAANNReg1_1_5, modelANNReg1_1_5, outputStandardizeANNReg1_1_5,
                                       listParticles[p]->m_T_gas, listParticles[p]->m_Yk_gas,
                                       numVarANN, nsp,
                                       input_childLocalStandardized, 
                                       outputStandardizeANN_Vec, mixture);
                            break;
                            
                            // Cluster 1_1_6
                        case 6:
                            if(rank==0)
                            {
                                cout << " Case1_1_6"<<  endl;
                            }
                            advanceANN_withoutPCA(cluster1_1_6_meanScaleINPUT, cluster1_1_6_stdScaleINPUT,
                                       cluster1_1_6_meanScaleOUTPUT, cluster1_1_6_stdScaleOUTPUT,
                                       cluster1_1_6_maxOUTPUT, cluster1_1_6_minOUTPUT,
                                       grandChild1_1_6_maxVal, grandChild1_1_6_minVal, // grandChild
                                       inputPCAANNReg1_1_6, modelANNReg1_1_6, outputStandardizeANNReg1_1_6,
                                       listParticles[p]->m_T_gas, listParticles[p]->m_Yk_gas,
                                       numVarANN, nsp,
                                       input_childLocalStandardized, 
                                       outputStandardizeANN_Vec, mixture);
                            break;
                            
                            // Cluster 1_1_7
                        case 7:
                            if(rank==0)
                            {
                                cout << " Case1_1_7"<<  endl;
                            }
                            advanceANN_withoutPCA(cluster1_1_7_meanScaleINPUT, cluster1_1_7_stdScaleINPUT,
                                       cluster1_1_7_meanScaleOUTPUT, cluster1_1_7_stdScaleOUTPUT,
                                       cluster1_1_7_maxOUTPUT, cluster1_1_7_minOUTPUT,
                                       grandChild1_1_7_maxVal, grandChild1_1_7_minVal, // grandChild
                                       inputPCAANNReg1_1_7, modelANNReg1_1_7, outputStandardizeANNReg1_1_7,
                                       listParticles[p]->m_T_gas, listParticles[p]->m_Yk_gas,
                                       numVarANN, nsp,
                                       input_childLocalStandardized, 
                                       outputStandardizeANN_Vec, mixture);
                            break;
                            
                            
                            // Cluster 1_1_8
                        case 8:
                            if(rank==0)
                            {
                                cout << " Case1_1_8"<<  endl;
                            }
                            advanceANN_withoutPCA(cluster1_1_8_meanScaleINPUT, cluster1_1_8_stdScaleINPUT,
                                       cluster1_1_8_meanScaleOUTPUT, cluster1_1_8_stdScaleOUTPUT,
                                       cluster1_1_8_maxOUTPUT, cluster1_1_8_minOUTPUT,
                                       grandChild1_1_8_maxVal, grandChild1_1_8_minVal, // grandChild
                                       inputPCAANNReg1_1_8, modelANNReg1_1_8, outputStandardizeANNReg1_1_8,
                                       listParticles[p]->m_T_gas, listParticles[p]->m_Yk_gas,
                                       numVarANN, nsp,
                                       input_childLocalStandardized, 
                                       outputStandardizeANN_Vec, mixture);
                            break;
                            
                            // Cluster 1_1_9
                        case 9:
                            if(rank==0)
                            {
                                cout << " Case1_1_9"<<  endl;
                            }
                            advanceANN_withoutPCA(cluster1_1_9_meanScaleINPUT, cluster1_1_9_stdScaleINPUT,
                                       cluster1_1_9_meanScaleOUTPUT, cluster1_1_9_stdScaleOUTPUT,
                                       cluster1_1_9_maxOUTPUT, cluster1_1_9_minOUTPUT,
                                       grandChild1_1_9_maxVal, grandChild1_1_9_minVal, // grandChild
                                       inputPCAANNReg1_1_9, modelANNReg1_1_9, outputStandardizeANNReg1_1_9,
                                       listParticles[p]->m_T_gas, listParticles[p]->m_Yk_gas,
                                       numVarANN, nsp,
                                       input_childLocalStandardized, 
                                       outputStandardizeANN_Vec, mixture);
                            break;
                            
                            
                            
                    } // END switch(indexGrandChildCluster) // cluster1_1_X
                    
                }// END if(indexChildCluster==1)
                else // // if(indexChildCluster==2) - childCluster1_2
                {
                    // Classifier grand child
                    // LOCAL CHILD Standardize Y ,T for Kmeans LOCAL CHILD: standardizedX = (X - meanScale) / standardDeviationScale
                    // Standardize T
                    input_childLocalStandardized[numVarANN-1] = listParticles[p]->m_T_gas - cluster1_2_meanScaleINPUT[numVarANN-1];
                    input_childLocalStandardized[numVarANN-1] /= cluster1_2_stdScaleINPUT[numVarANN-1];
                    // Standardize Y
                    for (int kANN=0; kANN<nsp-1; kANN++) // Don't take N2 (last species)
                    {
                        input_childLocalStandardized[kANN] = listParticles[p]->m_Yk_gas[kANN] - cluster1_2_meanScaleINPUT[kANN];
                        input_childLocalStandardized[kANN] /= cluster1_2_stdScaleINPUT[kANN];
                        //Cross-check if there is N2 inside, if yes, should move N2 to the end of the scheme
                        if(mixture->speciesName(kANN) == "N2") {cout << "Warning : N2 in the ANN!!!" << endl;}
                    }
                    
                    
                    // CLASSIFIER LOCAL CHILD: Caculate the distance : Eucledian (Kmeans) or Projection parition (LPCA)
                    // >> Determine indexGrandChildCluster
                    // LPCA example: cv::Mat d0 = (x-r0).t() * eigenVec0.t() * eigenVec0 * (x-r0);
                    // Eucledean: Can be calculated by std::vector
                    
                    // Declare vector x
                    xChild = cv::Mat(numVarANN, 1, CV_32F, input_childLocalStandardized.data());
                    d1_2_0 = cv::norm(xChild,r1_2_0,cv::NORM_L2SQR); // Calculate Euclidean distance from x to cluster 0
                    d1_2_1 = cv::norm(xChild,r1_2_1,cv::NORM_L2SQR);
                    d1_2_2 = cv::norm(xChild,r1_2_2,cv::NORM_L2SQR);
                    d1_2_3 = cv::norm(xChild,r1_2_3,cv::NORM_L2SQR);
                    d1_2_4 = cv::norm(xChild,r1_2_4,cv::NORM_L2SQR);
                    d1_2_5 = cv::norm(xChild,r1_2_5,cv::NORM_L2SQR);
                    d1_2_6 = cv::norm(xChild,r1_2_6,cv::NORM_L2SQR);
                    d1_2_7 = cv::norm(xChild,r1_2_7,cv::NORM_L2SQR);
                    d1_2_8 = cv::norm(xChild,r1_2_8,cv::NORM_L2SQR);
                    d1_2_9 = cv::norm(xChild,r1_2_9,cv::NORM_L2SQR);

                    
                    // Group distance to d
                    cv::Mat dChild;
                    dChild.push_back(d1_2_0);
                    dChild.push_back(d1_2_1);
                    dChild.push_back(d1_2_2);
                    dChild.push_back(d1_2_3);
                    dChild.push_back(d1_2_4);
                    dChild.push_back(d1_2_5);
                    dChild.push_back(d1_2_6);
                    dChild.push_back(d1_2_7);
                    dChild.push_back(d1_2_8);
                    dChild.push_back(d1_2_9);

                    
                    // Take position of mean (argmin) = indexParentCluster
                    cv::minMaxLoc(dChild, NULL, NULL, &minLocChild, NULL); // d, minVal, maxVal, minLoc, maxLoc
                    indexGrandChildCluster = minLocChild.y; //Point has 2 coordinates (x,y)
                    
                    if(rank==0) // Print to check
                    {
                        cout << "Particles " << p << " || parentCluster = " << indexParentCluster;
                        cout << " || childCluster = " << indexChildCluster << " || grandChildCluster = " << indexGrandChildCluster << " ";
                    }
                    
                    // Save indexCluster of this particle
                    if(activateIndexCluster)
                    {
                        if(rank==0)
                        {
                            ofstream indexClusterFile("outputs/indexCluster.dat",ios::app); //ios::app = append at the end of the file
                            indexClusterFile << t << "	";
                            indexClusterFile << p << "	";
                            indexClusterFile << indexParentCluster <<  "	";
                            indexClusterFile << indexChildCluster <<  "	";
                            indexClusterFile << indexGrandChildCluster << endl;
                            indexClusterFile.close();
                        }
                    }
                    
                    switch(indexGrandChildCluster) // cluster1_2_X
                    {
                        // Cluster 1_2_0
                        case 0:
                            if(rank==0)
                            {
                                cout << " Case1_2_0"<<  endl;
                            }
                            advanceANN_withoutPCA(cluster1_2_0_meanScaleINPUT, cluster1_2_0_stdScaleINPUT,
                                       cluster1_2_0_meanScaleOUTPUT, cluster1_2_0_stdScaleOUTPUT,
                                       cluster1_2_0_maxOUTPUT, cluster1_2_0_minOUTPUT,
                                       grandChild1_2_0_maxVal, grandChild1_2_0_minVal, // grandChild
                                       inputPCAANNReg1_2_0, modelANNReg1_2_0, outputStandardizeANNReg1_2_0,
                                       listParticles[p]->m_T_gas, listParticles[p]->m_Yk_gas,
                                       numVarANN, nsp,
                                       input_childLocalStandardized, 
                                       outputStandardizeANN_Vec, mixture);
                            break;
                            
                        // Cluster 1_2_1
                        case 1:
                            if(rank==0)
                            {
                                cout << " Case1_2_1"<<  endl;
                            }
                            advanceANN_withoutPCA(cluster1_2_1_meanScaleINPUT, cluster1_2_1_stdScaleINPUT,
                                       cluster1_2_1_meanScaleOUTPUT, cluster1_2_1_stdScaleOUTPUT,
                                       cluster1_2_1_maxOUTPUT, cluster1_2_1_minOUTPUT,
                                       grandChild1_2_1_maxVal, grandChild1_2_1_minVal, // grandChild
                                       inputPCAANNReg1_2_1, modelANNReg1_2_1, outputStandardizeANNReg1_2_1,
                                       listParticles[p]->m_T_gas, listParticles[p]->m_Yk_gas,
                                       numVarANN, nsp,
                                       input_childLocalStandardized, 
                                       outputStandardizeANN_Vec, mixture);
                            break;
                            
                            // Cluster 1_2_2
                        case 2:
                            if(rank==0)
                            {
                                cout << " Case1_2_2"<<  endl;
                            }
                            advanceANN_withoutPCA(cluster1_2_2_meanScaleINPUT, cluster1_2_2_stdScaleINPUT,
                                       cluster1_2_2_meanScaleOUTPUT, cluster1_2_2_stdScaleOUTPUT,
                                       cluster1_2_2_maxOUTPUT, cluster1_2_2_minOUTPUT,
                                       grandChild1_2_2_maxVal, grandChild1_2_2_minVal, // grandChild
                                       inputPCAANNReg1_2_2, modelANNReg1_2_2, outputStandardizeANNReg1_2_2,
                                       listParticles[p]->m_T_gas, listParticles[p]->m_Yk_gas,
                                       numVarANN, nsp,
                                       input_childLocalStandardized, 
                                       outputStandardizeANN_Vec, mixture);
                            break;
                            
                        // Cluster 1_2_3
                        case 3:
                            if(rank==0)
                            {
                                cout << " Case1_2_3"<<  endl;
                            }
                            advanceANN_withoutPCA(cluster1_2_3_meanScaleINPUT, cluster1_2_3_stdScaleINPUT,
                                       cluster1_2_3_meanScaleOUTPUT, cluster1_2_3_stdScaleOUTPUT,
                                       cluster1_2_3_maxOUTPUT, cluster1_2_3_minOUTPUT,
                                       grandChild1_2_3_maxVal, grandChild1_2_3_minVal, // grandChild
                                       inputPCAANNReg1_2_3, modelANNReg1_2_3, outputStandardizeANNReg1_2_3,
                                       listParticles[p]->m_T_gas, listParticles[p]->m_Yk_gas,
                                       numVarANN, nsp,
                                       input_childLocalStandardized, 
                                       outputStandardizeANN_Vec, mixture);
                            break;
                            
                            // Cluster 1_2_4
                        case 4:
                            if(rank==0)
                            {
                                cout << " Case1_2_4"<<  endl;
                            }
                            advanceANN_withoutPCA(cluster1_2_4_meanScaleINPUT, cluster1_2_4_stdScaleINPUT,
                                       cluster1_2_4_meanScaleOUTPUT, cluster1_2_4_stdScaleOUTPUT,
                                       cluster1_2_4_maxOUTPUT, cluster1_2_4_minOUTPUT,
                                       grandChild1_2_4_maxVal, grandChild1_2_4_minVal, // grandChild
                                       inputPCAANNReg1_2_4, modelANNReg1_2_4, outputStandardizeANNReg1_2_4,
                                       listParticles[p]->m_T_gas, listParticles[p]->m_Yk_gas,
                                       numVarANN, nsp,
                                       input_childLocalStandardized, 
                                       outputStandardizeANN_Vec, mixture);
                            break;
                            
                            // Cluster 1_2_5
                        case 5:
                            if(rank==0)
                            {
                                cout << " Case1_2_5"<<  endl;
                            }
                            advanceANN_withoutPCA(cluster1_2_5_meanScaleINPUT, cluster1_2_5_stdScaleINPUT,
                                       cluster1_2_5_meanScaleOUTPUT, cluster1_2_5_stdScaleOUTPUT,
                                       cluster1_2_5_maxOUTPUT, cluster1_2_5_minOUTPUT,
                                       grandChild1_2_5_maxVal, grandChild1_2_5_minVal, // grandChild
                                       inputPCAANNReg1_2_5, modelANNReg1_2_5, outputStandardizeANNReg1_2_5,
                                       listParticles[p]->m_T_gas, listParticles[p]->m_Yk_gas,
                                       numVarANN, nsp,
                                       input_childLocalStandardized, 
                                       outputStandardizeANN_Vec, mixture);
                            break;
                            
                            // Cluster 1_2_6
                        case 6:
                            if(rank==0)
                            {
                                cout << " Case1_2_6"<<  endl;
                            }
                            advanceANN_withoutPCA(cluster1_2_6_meanScaleINPUT, cluster1_2_6_stdScaleINPUT,
                                       cluster1_2_6_meanScaleOUTPUT, cluster1_2_6_stdScaleOUTPUT,
                                       cluster1_2_6_maxOUTPUT, cluster1_2_6_minOUTPUT,
                                       grandChild1_2_6_maxVal, grandChild1_2_6_minVal, // grandChild
                                       inputPCAANNReg1_2_6, modelANNReg1_2_6, outputStandardizeANNReg1_2_6,
                                       listParticles[p]->m_T_gas, listParticles[p]->m_Yk_gas,
                                       numVarANN, nsp,
                                       input_childLocalStandardized, 
                                       outputStandardizeANN_Vec, mixture);
                            break;
                            
                            
                            // Cluster 1_2_7
                        case 7:
                            if(rank==0)
                            {
                                cout << " Case1_2_7"<<  endl;
                            }
                            advanceANN_withoutPCA(cluster1_2_7_meanScaleINPUT, cluster1_2_7_stdScaleINPUT,
                                       cluster1_2_7_meanScaleOUTPUT, cluster1_2_7_stdScaleOUTPUT,
                                       cluster1_2_7_maxOUTPUT, cluster1_2_7_minOUTPUT,
                                       grandChild1_2_7_maxVal, grandChild1_2_7_minVal, // grandChild
                                       inputPCAANNReg1_2_7, modelANNReg1_2_7, outputStandardizeANNReg1_2_7,
                                       listParticles[p]->m_T_gas, listParticles[p]->m_Yk_gas,
                                       numVarANN, nsp,
                                       input_childLocalStandardized, 
                                       outputStandardizeANN_Vec, mixture);
                            break;
                            
                            // Cluster 1_2_8
                        case 8:
                            if(rank==0)
                            {
                                cout << " Case1_2_8"<<  endl;
                            }
                            advanceANN_withoutPCA(cluster1_2_8_meanScaleINPUT, cluster1_2_8_stdScaleINPUT,
                                       cluster1_2_8_meanScaleOUTPUT, cluster1_2_8_stdScaleOUTPUT,
                                       cluster1_2_8_maxOUTPUT, cluster1_2_8_minOUTPUT,
                                       grandChild1_2_8_maxVal, grandChild1_2_8_minVal, // grandChild
                                       inputPCAANNReg1_2_8, modelANNReg1_2_8, outputStandardizeANNReg1_2_8,
                                       listParticles[p]->m_T_gas, listParticles[p]->m_Yk_gas,
                                       numVarANN, nsp,
                                       input_childLocalStandardized, 
                                       outputStandardizeANN_Vec, mixture);
                            break;
                            
                            // Cluster 1_2_9
                        case 9:
                            if(rank==0)
                            {
                                cout << " Case1_2_9"<<  endl;
                            }
                            advanceANN_withoutPCA(cluster1_2_9_meanScaleINPUT, cluster1_2_9_stdScaleINPUT,
                                       cluster1_2_9_meanScaleOUTPUT, cluster1_2_9_stdScaleOUTPUT,
                                       cluster1_2_9_maxOUTPUT, cluster1_2_9_minOUTPUT,
                                       grandChild1_2_9_maxVal, grandChild1_2_9_minVal, // grandChild
                                       inputPCAANNReg1_2_9, modelANNReg1_2_9, outputStandardizeANNReg1_2_9,
                                       listParticles[p]->m_T_gas, listParticles[p]->m_Yk_gas,
                                       numVarANN, nsp,
                                       input_childLocalStandardized, 
                                       outputStandardizeANN_Vec, mixture);
                            break;
                            
                            
                            
                    } // END switch(indexGrandChildCluster) // cluster1_2_X
                    
                } // END if(indexChildCluster==2)
                
                // Release OPENCV MAT
                dParent.release();
                
            } // END parentCluster1
            
            // Release OPENCV MAT
            d.release();
            
            
        } // END for (int p=0; p<nTot; p++)
        
 
        
        if(rank==0) cout << "Number of particles modified by EMST check = " << modifEMSTParticle << endl;
        
    }// END flagANN
	
	if(rank==0)	//HT@2020.08.20 check after test
        {
                cout << "Particle 700 before - T= " << listParticles[700]->m_T_gas << " |Enthalpy = " <<  listParticles[700]->m_H_gas << endl;
        }


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
    r0.release();
    r1.release();
    
    r0_0.release();
    r0_1.release();
    r0_2.release();
    r0_3.release();
    r0_4.release();
    r0_5.release();
    r0_6.release();
    r0_7.release();
    r0_8.release();
    r0_9.release();
    r0_10.release();
    r0_11.release();
    r0_12.release();
    r0_13.release();
    r0_14.release();
    r0_15.release();
    r0_16.release();
    r0_17.release();
    r0_18.release();
    r0_19.release();
    
    r1_0.release();
    r1_1.release();
    r1_2.release();

    r1_0_0.release();
    r1_0_1.release();
    r1_0_2.release();
    r1_0_3.release();
    r1_0_4.release();
    r1_0_5.release();
    r1_0_6.release();
    r1_0_7.release();
    r1_0_8.release();
    r1_0_9.release();
    r1_0_10.release();
    r1_0_11.release();
    r1_0_12.release();
    r1_0_13.release();
    r1_0_14.release();
    r1_0_15.release();
    r1_0_16.release();
    r1_0_17.release();
    r1_0_18.release();
    r1_0_19.release();

    r1_1_0.release();
    r1_1_1.release();
    r1_1_2.release();
    r1_1_3.release();
    r1_1_4.release();
    r1_1_5.release();
    r1_1_6.release();
    r1_1_7.release();
    r1_1_8.release();
    r1_1_9.release();

    r1_2_0.release();
    r1_2_1.release();
    r1_2_2.release();
    r1_2_3.release();
    r1_2_4.release();
    r1_2_5.release();
    r1_2_6.release();
    r1_2_7.release();
    r1_2_8.release();
    r1_2_9.release();
    // HuuTri@20211006 : Commented - ANN without PCA
    /*input_childLocalStandardized_Mat.release();
    inputPCA.release();
    */// HuuTri@20211006 : Commented - ANN without PCA

    x.release();
    xParent.release();
   
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
}


computeMultipleInlet::~computeMultipleInlet() //Destructeur
{}

