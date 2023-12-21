#include "../main/conditions.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iterator>

void read_string_array(string key, string value[20], string line)
{
	int i =0;
	stringstream s_line;
	s_line << line;
	string k;
	s_line >> k;
	while (k == key)
	{
		s_line >> k;
		while (k == " ") s_line >> k;
		s_line >> value[i];
		i++;
	}
}

void read_double(string key, double* value, string line)
{
	stringstream s_line;
	s_line << line;
	string k;
	s_line >> k;
	if (k == key)
	{
		s_line >> k;
		while (k == " ") s_line >> k;
		s_line >> *value;
	}
}

void read_string(string key, string* value, string line)
{
	stringstream s_line;
	s_line << line;
	string k;
	s_line >> k;
	if (k == key)
	{
		s_line >> k;
		while (k == " ") s_line >> k;
		s_line >> *value;
	}
}

void read_compo(string key, string* value, string line)
{
	stringstream s_line;
	s_line << line;
	string k;
	s_line >> k;
	if (k == key)
	{
		s_line >> k;
		while (k == " ") s_line >> k;
		getline(s_line, *value);
	}
}

void read_vector(string key, vector<string> *value, string line)
{
	ostringstream s_line;
	stringstream s;
	s << line;
	string k;
	s >> k;
	if (k == key)
	{
		s >> k;
		while (k == " ") s >> k;
		s >> k;
		(*value).push_back(k);
	}
}

void read_integer(string key, int* value, string line)
{
	stringstream s_line;
	s_line << line;
	string k;
	s_line >> k;
	if (k == key)
	{
		s_line >> k;
		while (k == " ") s_line >> k;
		s_line >> *value;
	}
	
}

void read_bool(string key, bool* value, string line)
{
	stringstream s_line;
	s_line << line;
	string k;
	s_line >> k;
	if (k == key)
	{
		s_line >> k;
		while (k == " ") s_line >> k;
		s_line >> *value;
	}
	
}

void read_vector_flame(string key, vector<PremixedFlames*> *listFlames, PremixedFlames* Pmx, string line)
{
	stringstream s_line;
	s_line << line;
	string k;
	s_line >> k;
	if (k == key)
	{
		(*listFlames).push_back(Pmx);
	}
}

void read_vector_inlet(string key, vector<MultipleInlet*> *listInlets, MultipleInlet* MI, string line)
{
	stringstream s_line;
	s_line << line;
	string k;
	s_line >> k;
	if (k == key)
	{
		(*listInlets).push_back(MI);
	}
}

void read_vector_inletGB(string key, vector<MultipleInlet*> *listInlets, Characteristics_MultipleInlet* MGB, string line)
{
	stringstream s_line;
	s_line << line;
	string k;
	s_line >> k;
	if (k == key)
	{
		(*listInlets).push_back(MGB);
	}
}

void conditions(int &debuglevel,
	vector <string> &speciesToPlot, //Trajectories the user wants to plot
	string &mech_ref, //Reference detailed mechanism
	string &mech, //Current step reference chemical scheme
	string &mech_desc, //Description of the reference chemical scheme
	string &configuration, //Studied combustion regime - OPTIONS: "MultipleInlet"; "PremixedFlames";
	vector<MultipleInlet*> &listInlets, //List of flames with their characteristics for the multiple inlet problem
	vector<PremixedFlames*> &listFlames, //List of premixed flames with their characteristics
	vector<string> &listTargets, //List of target species
	string &step, //Step to perform - OPTIONS: "DRGEP_Species"; "DRGEP_Reactions"; "ComputeTrajectories"; "computeQSSCriteria"; "getQSSfile"; "getQSSfileFORTRAN"; "Optimisation"; "Lumping";
	bool &new_mixing, //Define if a new mixing of the particles is defined or if the mixing used with the previous step is kept
	bool &plot_T,
	bool &plot_U,
	vector<QSSscenario*> &listQSSscenarios, //List of scenarios to test with different QSS species 
	OptimScenario* &listOptimScenarios, //List of scenarios to optimise
	vector<string> &trajectory_ref) //Name of a reference trajectory defined by the user for the plots 
{
	ifstream input("input_file.ini");
	
	//Flame parameters
	double Tf = 0.0;
	double To = 0.0;
	double P = 0.0;
	double Phi = 0.0;
	string Yf = "";
	string Xf = "";
	string Yo = "";
	string Xo = "";
	string Initial_flame = "";
	string Final_flame = "";
	
	//Inlet parameters
	double T=0.0;
	double MassFlowRate = 0.0;
	string Yk = "";
	string Xk = "";
	int EvaporationModel = 0;
	double DropletDiameter = 0.0;
	double EvaporationTime = 0.0;
	double liquidDensity = 0.0;
	double EvaporationLatentHeat = 0.0;
	MultipleInlet* MI; 
	Characteristics_MultipleInlet* MGB; 
	double MixingTime = 0.0;
	double TimeStep = 0.0;
	int IterationNumber = 0;
	int EMST = 0;
	
	int nbFlame=0;
	int nbInlet=0;
	int PopSize = 0;
	int MaxAllowableGenerations = 0;
	int NbElitism = 0;
	double CrossoverRate = 0.0;
	double MutationRate = 0.0;
	vector<string> array1;
	PremixedFlames* Pmx;
	
	double AllowedVariation_A = 0.0;
	double AllowedVariation_b = 0.0;
	double AllowedVariation_E = 0.0;
	
	for (std::string line; getline(input,line);)
	{
		read_bool("NewMixing", &new_mixing, line);
		read_double("MixingTime", &MixingTime, line);
		read_double("TimeStep", &TimeStep, line);
		read_integer("IterationNumber", &IterationNumber, line);
		read_integer("EMST", &EMST, line);
		
		read_integer("debuglevel",&debuglevel, line);
		
		read_string("configuration", &configuration, line);
		read_string("step", &step, line);
		
		read_string("mech", &mech, line);
		
		read_string("mech_desc", &mech_desc, line);
		
		read_vector("speciesToPlot", &speciesToPlot, line);
		
		read_bool("plot_T",&plot_T, line);
		read_bool("plot_U",&plot_U, line);
		
		read_string("mech_ref", &mech_ref, line);
		
		if ((mech_ref != "None")&&(configuration != "PremixedFlames"))
		{
			read_vector("trajectory_ref", &trajectory_ref, line);
		}
		else
		{
			trajectory_ref.push_back("");
		}
		
		if (configuration == "MultipleInlet") {
			read_double("T", &T, line);
			read_double("MassFlowRate", &MassFlowRate, line);
			read_double("Pressure", &P, line);
			read_compo("Yk", &Yk, line);
			read_compo("Xk", &Xk, line);
			read_integer("EvaporationModel", &EvaporationModel, line);
			read_double("DropletDiameter", &DropletDiameter, line);
			read_double("EvaporationTime", &EvaporationTime, line);
			read_double("LiquidDensity", &liquidDensity, line);
			read_double("EvaporationLatentHeat", &EvaporationLatentHeat, line);
		} else {
			read_string("Initial_flame", &Initial_flame, line);
			read_string("Final_flame", &Final_flame, line);
			read_compo("Yf", &Yf, line);
			read_compo("Xf", &Xf, line);
			read_compo("Yo", &Yo, line);
			read_compo("Xo", &Xo, line);
			read_double("Equivalence_ratio", &Phi, line);
			read_double("Tf", &Tf, line);
			read_double("To", &To, line);
			read_double("Pressure", &P, line);
		}
		
		//If a new Premixed flame is found
		if (line == "//End")
		{
			nbFlame++;
			Pmx = new PremixedFlames(Tf, To, P, Phi, Yf, Xf, Yo, Xo, Initial_flame, Final_flame);
			read_vector_flame("//End", &listFlames, Pmx, line);
		}
			
		//If a new Inlet is found
		if (line == "//EndInlet")
		{
			nbInlet++;
			bool isEvap = (EvaporationModel > 0)?true:false;
			MI = new MultipleInlet(T, P, MassFlowRate, Xk, Yk, isEvap, DropletDiameter, EvaporationTime, liquidDensity, EvaporationLatentHeat);
			read_vector_inlet("//EndInlet", &listInlets, MI, line);
		}
		
		//Inlet of Burned Gases
		if (line == "//EndInletGB")
		{
			nbInlet++;
			bool isEvap = (EvaporationModel > 0)?true:false;
			bool isEMST = (EMST > 0)?true:false;
			MGB = new Characteristics_MultipleInlet(T, P, MassFlowRate, Xk, Yk, isEvap, DropletDiameter, EvaporationTime, liquidDensity, EvaporationLatentHeat, MixingTime, TimeStep, IterationNumber, true, isEMST);
			read_vector_inletGB("//EndInletGB", &listInlets, MGB, line);
		}
		read_vector("listTargets", &listTargets, line);
		
		read_vector("QSS", &array1, line);
		
		read_integer("PopSize",&PopSize, line);
		read_integer("MaxAllowableGenerations",&MaxAllowableGenerations, line);
		read_integer("NbElitism",&NbElitism, line);
		
		read_double("CrossoverRate", &CrossoverRate, line);
		read_double("MutationRate", &MutationRate, line);
		
		read_double("AllowedVariation_A", &AllowedVariation_A, line);
		read_double("AllowedVariation_b", &AllowedVariation_b, line);
		read_double("AllowedVariation_E", &AllowedVariation_E, line);
	}
	input.close();
	
	//------QSS step------//
	vector<string> vec(array1);
	listQSSscenarios.push_back(new QSSscenario(vec));
	
	//------Optimisation------//
	listOptimScenarios = new OptimScenario(PopSize, MaxAllowableGenerations, NbElitism, CrossoverRate, MutationRate, AllowedVariation_A, AllowedVariation_b, AllowedVariation_E);
}


