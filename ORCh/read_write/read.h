#ifndef Read_H
#define Read_H

//using namespace std;

#include <sys/stat.h>
#include "species.h"
#include "reaction.h"


//--------------------
class Read
{
   public:
   //constructeur
   Read();

   //Read the list of chemical species present whithin 
   virtual void Read_species(string mech, vector<Species_ORCh*> &listSpecies) const;

   virtual void Read_reactions(string mech, vector<Reaction_ORCh*> &listReactions) const;

   virtual void Input_file_into_string(string mech, string& mech_string) const;

   virtual void translate_type(string initial_string, string type, string &associated_text) const;

   virtual void find_XML_key(string& initial_string, string& substracted_part, string keyword) const;

   virtual void find_XML_key(string& initial_string, string& substracted_part, string keyword,
                             string& type) const;

   virtual void find_XML_key(string& initial_string, string& substracted_part, string keyword,
                             string& type, bool &found) const;

   virtual void Species_into_Reactants_Products(string initial_string, vector<string>& Species_List,
                                                vector<double>& Coeffs) const;

   virtual void Atom_into_Species(string initial_string, vector<string>& Atom_List,
                                                vector<int>& Coeffs) const;

   virtual void Species_into_Efficiencies(string initial_string, vector<string>& Species_List,
                                          vector<double>& Coeffs) const;

   virtual void getTroeCoefficients(string initial_string, vector<double> &Coefficients) const;


   virtual void get_Arrhenius_coefficients(string Arrhenius_description, 
                                      double &A, double &b, double &E, string &Etype) const;


   //destructeur
   virtual ~Read();









   private:



};











#endif
