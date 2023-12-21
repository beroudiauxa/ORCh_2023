#ifndef Species_H
#define Species_H

#include <iostream>
#include <vector>

using namespace std;


//--------------------
class Species_ORCh
{
   public:
   //constructeur
   Species_ORCh(string Name, string Description, 
           vector<double> NASACoeffs_lowT, double NASA_lowT_min, double NASA_lowT_max,
           vector<double> NASACoeffs_highT, double NASA_highT_min, double NASA_highT_max,
           int nb_C, int nb_H, int nb_O, int nb_N, int nb_Ar,
           string geometry, double LJ_welldepth, double LJ_diameter, 
           double dipoleMoment, double polarizability, double rot_relax,
           bool QSS); 

   //destructeur
   virtual ~Species_ORCh();

   public:
   string m_Name;
   string m_Description; //This is a copy of the species description from the xml file

   vector<double> m_NASACoeffs_lowT;
   double m_NASA_lowT_min;
   double m_NASA_lowT_max;
   vector<double> m_NASACoeffs_highT;
   double m_NASA_highT_min;
   double m_NASA_highT_max;

   int m_C;
   int m_H;
   int m_O;
   int m_N;
   int m_Ar;

   string m_geometry;
   double m_LJ_welldepth;
   double m_LJ_diameter;
   double m_dipoleMoment;
   double m_polarizability;
   double m_rot_relax;

   bool m_QSS;

};


#endif

