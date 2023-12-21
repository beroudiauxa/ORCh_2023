#ifndef Lumping_H
#define Lumping_H

#include <iostream>
#include <vector>

using namespace std;


//--------------------
class Lumping
{
   public:
   //constructeur
   Lumping(string Name, int nb_C, int nb_H, int nb_O, int nb_N, vector<string> SpeciesList, vector<double> alphaSpecies, vector<int> SpeciesNumber, vector<vector<int> > ReactionsGroups, vector<vector<string> > ReactionsGroupsDescription); 

   //destructeur
   virtual ~Lumping();

   public:
   string m_Name;

   int m_C;
   int m_H;
   int m_O;
   int m_N;

   vector<string> m_SpeciesList;
   vector<double> m_alphaSpecies;
   vector<int> m_SpeciesNumber;
   vector<vector<int> > m_ReactionsGroups;
   vector<vector<string> > m_ReactionsGroupsDescription;

};

#endif

