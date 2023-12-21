#ifndef LUMPING_H
#define LUMPING_H

#include <mpi.h>
#include <iostream>

#include "../read_write/read.h"
#include "lumping_species.h"
#include <Cantera.h>
#include <user.h>

using namespace Cantera;
using namespace Cantera_CXX;
using namespace User;


using namespace std;

//---------------------
class lumping
{
   public:
   //constructeur
   lumping();

   virtual void lumping_mech(string mech, string trajectory_ref, string mech_desc) const;
   virtual void lumping_groups(vector<Lumping*> &listLumpingGroups, vector<Species_ORCh*> listSpecies_ref) const;
   virtual void lumping_species(vector<Lumping*> listLumpingGroups, vector<Species_ORCh*> listSpecies_ref, vector<Species_ORCh*> &listSpecies_lumping, vector<Reaction_ORCh*> &listReactions_lumping, vector<bool> &Species_to_add) const;
   virtual void lumping_equation(vector<Reaction_ORCh*> &listReactions_lumping) const;
   virtual void lumping_optimise(vector<vector<double> > min_max_A, vector<vector<double> > min_max_b, vector<vector<double> > min_max_E,
                               vector<Reaction_ORCh*> &listReactions_lumping, vector<int> AssociatedReactions, int coeff, vector<double> t_ref, vector<double> T_ref, vector<double> conc_lumped, vector<double> ReactionRate, double max_ReactionRate, bool singleReaction, int jt) const;

   virtual void getReactionRate(vector<Lumping*> listLumpingGroups, vector<Reaction_ORCh*> listReactions_ref, vector<Species_ORCh*> listSpecies_ref, vector<double> &ReactionRate, double &max_ReactionRate, vector<vector<double> > conc_ref, vector<double> T_ref, vector<double> &conc_lumped, vector<int> AssociatedReactions, vector<string> AssociatedSpecies, int &coeff, bool singleReaction, int jt) const;


   //destructeur 
   virtual ~lumping();

   private:


};

#endif

