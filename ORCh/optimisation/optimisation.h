#ifndef Optim_H
#define Optim_H

#include <mpi.h>
#include <iostream>
#include "../read_write/write_QSS.h"
#include "../read_write/write.h"
#include "OptimScenario.h"
#include "../tools/fitness_criteria.h"
#include "../tools/outputs.h"


using namespace std;

//#define POP_SIZE 24
//#define MAX_ALLOWABLE_GENERATIONS 500

//#define CROSSOVER_RATE 0.75
#define MUTATION_RATE 0.02

#define NB_ELITISM 1


//Define the required accuracy for the studied variables
#define accuracy_A 0.0001
#define accuracy_Ea_R 0.0001
#define accuracy_b 0.0001




//--------------------
class Chromosome
{
   public:
   //constructeur
   Chromosome(string bits, double fitness);

   //destructeur
   virtual ~Chromosome();

   public:
   string m_bits;
   float m_fitness;
};


//--------------------
class Optim
{
   public:
   //constructeur
   Optim();

   virtual void Optimise(string mech_ref, string mech, string mech_desc, vector<bool> QSS_Species, OptimScenario* listOptimScenarios, int nbInlets_nbFlames, vector<string> listTargets, string configuration, vector<string> trajectory_ref);

   virtual void Rank(vector<Chromosome*>& Sort, vector<Chromosome*> Population, int PopSize);

   virtual void Linear_Ranking(vector<Chromosome*>& Sort, int PopSize);

   virtual string Roulette(int TotalFitness, vector<Chromosome*> Population, int PopSize);
   
   virtual double Get_random();

   virtual string GetRandomBits(int length) ;

   virtual double BinToDec(string bits);

   virtual string DecToBin(int number);

   virtual void Crossover(string &offspring1, string &offspring2, double CrossoverRate);

   virtual void Mutate(string &bits, double MutationRate);

   virtual void Elitism(int cPop, vector<Chromosome*>& Temporary, vector<Chromosome*> Sort, int PopSize);

   virtual bool replace(std::string& str, const std::string& from, const std::string& to);

   virtual void Translate_Bits(vector<double>& A, vector<double>& b, vector<double>& E, 
      vector<int> NBits_A, vector<int> NBits_b, vector<int> NBits_Ea_R, string Bits,
      vector<vector<double> > min_max_A, 
      vector<vector<double> > min_max_b,
      vector<vector<double> > min_max_Ea_R);

   virtual void Min_Max(OptimScenario* listOptimScenarios, vector<Reaction_ORCh*>& listReactions, 
      vector<vector<double> >& min_max_A, 
      vector<vector<double> >& min_max_b, 
      vector<vector<double> >& min_max_Ea_R);

   virtual void CreateAnalyticDirectory(string configuration, string fileName, int nbInlets_nbFlames, string step, vector<string> trajectory_ref);


   //destructeur
   virtual ~Optim();

   private:

};



#endif




