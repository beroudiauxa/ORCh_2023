#ifndef Write_QSS_H
#define Write_QSS_H

#include "read.h"

using namespace std;


//--------------------
class Write_QSS
{
   public:
   //constructeur
   Write_QSS();

   virtual int Get_species_number(string species, vector<Species_ORCh*> listSpecies) const;

   virtual void Write_QSS_file(string Dimension, string mech, string write_mech, vector<bool> &QSS_Species, bool optim, vector<double> A_Arrh, vector<double> b_Arrh, vector<double> E_Arrh) const;

   virtual void Check_Non_Linearity(string mech, vector<vector<double> > QSS_Criteria, int rank=0) const;

   virtual string print_sign(double input) const;


   //destructeur
   virtual ~Write_QSS();


   private:




};

#endif
















