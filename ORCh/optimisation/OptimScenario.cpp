#include "OptimScenario.h"

//---OptimScenario---

OptimScenario::OptimScenario(int PopSize, int MaxAllowableGenerations, int NbElitism, 
                             double CrossoverRate, double MutationRate, double AllowedVariation_A, double AllowedVariation_b, double AllowedVariation_E
                        ) //Constructeur
   :m_PopSize(PopSize), m_MaxAllowableGenerations(MaxAllowableGenerations), m_NbElitism(NbElitism), m_CrossoverRate(CrossoverRate), m_MutationRate(MutationRate),
    m_AllowedVariation_A(AllowedVariation_A), m_AllowedVariation_b(AllowedVariation_b), m_AllowedVariation_E(AllowedVariation_E)
{}

OptimScenario::~OptimScenario() //Destructeur
{}


