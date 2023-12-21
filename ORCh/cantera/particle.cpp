#include "particle.h"

//---Particle---

Particle::Particle(vector<double> Yk_gas, double T_gas, double H_gas, double Z_gas, double density_gas, double N_droplets, double droplets_diameter, double P_gas_liquid, vector<double> Yk_liquid, double density_liquid, double EvaporationLatentHeat) //Constructeur
   :m_Yk_gas(Yk_gas), m_T_gas(T_gas), m_H_gas(H_gas), m_Z_gas(Z_gas), m_density_gas(density_gas), m_N_droplets(N_droplets), m_droplets_diameter(droplets_diameter), m_P_gas_liquid(P_gas_liquid), m_Yk_liquid(Yk_liquid),
    m_density_liquid(density_liquid), m_EvaporationLatentHeat(EvaporationLatentHeat)
{}

Particle::~Particle() //Destructeur
{}


