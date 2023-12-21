#ifndef OptimScenario_H
#define OptimScenario_H

#include <iostream>
#include <vector>

using namespace std;


//--------------------
class OptimScenario
{
   public:
   //constructeur
   OptimScenario(int PopSize, int MaxAllowableGenerations, int NbElitism, double CrossoverRate, double MutationRate, double AllowedVariation_A, 
                 double AllowedVariation_b, double AllowedVariation_E); 

   //destructeur
   virtual ~OptimScenario();

   public:
   int m_PopSize;
   int m_MaxAllowableGenerations;
   int m_NbElitism;

   double m_CrossoverRate;
   double m_MutationRate;

   double m_AllowedVariation_A;
   double m_AllowedVariation_b;
   double m_AllowedVariation_E;

};


#endif

