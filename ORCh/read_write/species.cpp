#include "species.h"

//---Species---

Species_ORCh::Species_ORCh(string Name, string Description,
           vector<double> NASACoeffs_lowT, double NASA_lowT_min, double NASA_lowT_max,
           vector<double> NASACoeffs_highT, double NASA_highT_min, double NASA_highT_max,
           int nb_C, int nb_H, int nb_O, int nb_N, int nb_Ar,
           string geometry, double LJ_welldepth, double LJ_diameter, 
           double dipoleMoment, double polarizability, double rot_relax,
           bool QSS) //Constructeur
   :m_Name(Name), m_Description(Description),
    m_NASACoeffs_lowT(NASACoeffs_lowT), 
    m_NASA_lowT_min(NASA_lowT_min), m_NASA_lowT_max(NASA_lowT_max),
    m_NASACoeffs_highT(NASACoeffs_highT), 
    m_NASA_highT_min(NASA_highT_min), m_NASA_highT_max(NASA_highT_max),
    m_C(nb_C), m_H(nb_H), m_O(nb_O), m_N(nb_N), m_Ar(nb_Ar),
    m_geometry(geometry), m_LJ_welldepth(LJ_welldepth), m_LJ_diameter(LJ_diameter), 
    m_dipoleMoment(dipoleMoment), m_polarizability(polarizability), m_rot_relax(rot_relax),
    m_QSS(QSS)
{}

Species_ORCh::~Species_ORCh() //Destructeur
{}


