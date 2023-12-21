#ifndef Write_H
#define Write_H

#include "read.h"

using namespace std;



//--------------------
class Write
{
   public:
   //constructeur
   Write();

   virtual void Write_xml_file(string mech, string mech_desc, string write_mech, 
                               vector<bool> Species_to_add, vector<bool> Reactions_to_add, bool optim, vector<double> A_Arrh, vector<double> b_Arrh, vector<double> E_Arrh,
                               vector<Species_ORCh*> listSpecies_lumping, vector<Reaction_ORCh*> listReactions_lumping) const; 

   virtual void Write_reaction(ofstream& write, Reaction_ORCh* reac, vector<bool> Species_to_add, 
                               vector<Species_ORCh*> listSpecies, int j, vector<Reaction_ORCh*> listReactions, vector<bool> Reactions_fwd_to_add, 
                           vector<bool> Reactions_rev_to_add) const;


   virtual void Write_xml_for_Analytic_Applications(string mech, string mech_desc, 
                    string write_mech, vector<bool> Species_to_add) const;

   //destructeur
   virtual ~Write();


   private:




};











#endif

