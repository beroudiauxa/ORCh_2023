#ifndef Particle_H
#define Particle_H

#include <iostream>
#include <vector>

using namespace std;


//--------------------
class Particle
{
   public:
   //constructeur
   Particle(vector<double> Yk_gas, double T_gas, double H_gas, double Z_gas, double density_gas, double N_droplets, double droplets_diameter, double P_gas_liquid, vector<double> Yk_liquid, double density_liquid, 
            double EvaporationLatentHeat);

   //destructeur
   virtual ~Particle();

   public:

   vector<double> m_Yk_gas;
   double m_T_gas;
   double m_H_gas;
   double m_Z_gas;
   double m_density_gas;

   double m_N_droplets;
   double m_droplets_diameter;

   double m_P_gas_liquid;

   vector<double> m_Yk_liquid;
   double m_density_liquid;
   double m_EvaporationLatentHeat;

};


#endif

