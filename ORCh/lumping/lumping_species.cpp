#include "lumping_species.h"

//---Lumping---

Lumping::Lumping(string Name, int nb_C, int nb_H, int nb_O, int nb_N, vector<string> SpeciesList, vector<double> alphaSpecies, vector<int> SpeciesNumber, vector<vector<int> > ReactionsGroups, vector<vector<string> > ReactionsGroupsDescription) //Constructeur
   :m_Name(Name), m_C(nb_C), m_H(nb_H), m_O(nb_O), m_N(nb_N), m_SpeciesList(SpeciesList), m_alphaSpecies(alphaSpecies), m_SpeciesNumber(SpeciesNumber), m_ReactionsGroups(ReactionsGroups), m_ReactionsGroupsDescription(ReactionsGroupsDescription)
{}

Lumping::~Lumping() //Destructeur
{}




