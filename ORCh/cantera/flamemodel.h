#ifndef FlameModel_H
#define FlameModel_H

#include <iostream>
#include <vector>

using namespace std;

//--------------------
class AutoIgnition
{
   public:
   //constructeur
   AutoIgnition(double T_fuel, double T_oxidizer, double Pressure, double delta_t, double Phi, double max_t, string Y_fuel, string X_fuel ,string Y_oxidizer, string X_oxidizer);

   //destructeur
   virtual ~AutoIgnition();

   public:
   double m_T_fuel;
   double m_T_oxidizer;
   double m_Pressure;
   double m_delta_t;
   double m_max_t;
   double m_Phi;
   string m_Y_fuel;
   string m_X_fuel;
   string m_Y_oxidizer;
   string m_X_oxidizer;

};


//--------------------
class Premixed
{
   public:
   //constructeur
   Premixed();

   //destructeur
   virtual ~Premixed();

   public:

};


//--------------------
class MultipleInlet
{
   public:
   //constructeur
   MultipleInlet(double Temperature, double Pressure, double flowRate, string X_Species, string Y_Species, bool liquid, double DropletsDiameter, double Tau_vj, double density_liquid, double EvaporationLatentHeat);

   //destructeur
   virtual ~MultipleInlet();

   public:
   double m_Temperature;
   double m_Pressure;
   double m_flowRate;
   string m_X_Species;
   string m_Y_Species;
   bool m_liquid;
   double m_DropletsDiameter; //[m]
   double m_Tau_vj; //[s]
   double m_density_liquid; //[kg/m3]
   double m_EvaporationLatentHeat; //[J/kg]

};


//--------------------
class Characteristics_MultipleInlet : public MultipleInlet
{
   public:
   //constructeur
   Characteristics_MultipleInlet(double Temperature, double Pressure, double flowRate, string X_Species, string Y_Species, bool liquid, double DropletsDiameter, double Tau_vj, double density_liquid, double EvaporationLatentHeat, 
                                 double tau_t, double delta_t, int NbIterations, bool BurnedGases, bool EMST);

   //destructeur
   virtual ~Characteristics_MultipleInlet();

   public:
   double m_tau_t;
   double m_delta_t;
   int m_NbIterations;
   bool m_BurnedGases;
   bool m_EMST;
};

//--------------------
class PremixedFlames
{
   public:
   //constructeur
   PremixedFlames(double T_fuel, double T_oxidizer, double Pressure, double Phi, string Yf_str, string Xf_str, string Yo_str, string Xo_str, string Initial_Flame, string Final_Flame);

   //destructeur
   virtual ~PremixedFlames();

   public:
   double m_T_fuel;
   double m_T_oxidizer;
   double m_Pressure;
   double m_Phi;
   string m_Yf_str;
   string m_Xf_str;
   string m_Yo_str;
   string m_Xo_str;
   string m_Initial_Flame;
   string m_Final_Flame;

};

#endif

