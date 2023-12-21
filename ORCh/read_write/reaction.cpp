
#include "reaction.h"


//---Reaction---

Reaction_ORCh::Reaction_ORCh(double A, double b, double E, string ID, string equation,
                   string Etype, bool reversible, bool duplicate,
                   vector<string> ReactantSpecies, vector<string> ProductSpecies,
                   vector<double> ReactantStoichCoeffs, vector<double> ProductStoichCoeffs
                   ) //Constructeur
   :m_A(A), m_b(b), m_E(E), m_ID(ID), m_equation(equation), 
    m_Etype(Etype), m_reversible(reversible), m_duplicate(duplicate),
    m_ReactantSpecies(ReactantSpecies), m_ProductSpecies(ProductSpecies),
    m_ReactantStoichCoeffs(ReactantStoichCoeffs), m_ProductStoichCoeffs(ProductStoichCoeffs)
{}

Reaction_ORCh::~Reaction_ORCh() //Destructeur
{}





//---Simple---

Simple::Simple(double A, double b, double E, string ID, string equation,
               string Etype, bool reversible, bool duplicate,
               vector<string> ReactantSpecies, vector<string> ProductSpecies,
               vector<double> ReactantStoichCoeffs, vector<double> ProductStoichCoeffs
               ) //Constructeur
   :Reaction_ORCh(A, b, E, ID, equation, Etype, reversible, duplicate, ReactantSpecies, ProductSpecies,
             ReactantStoichCoeffs, ProductStoichCoeffs)
{}

Simple::~Simple() //Destructeur
{}





//---ThreeBody---

ThreeBody::ThreeBody(double A, double b, double E, string ID, string equation,
               string Etype, bool reversible, bool duplicate,
               vector<string> ReactantSpecies, vector<string> ProductSpecies,
               vector<double> ReactantStoichCoeffs, vector<double> ProductStoichCoeffs,
               /*ThreeBody*/
               vector<double> TBconc, vector<string> NameTBconc
               ) //Constructeur
   :Reaction_ORCh(A, b, E, ID, equation, Etype, reversible, duplicate, ReactantSpecies, ProductSpecies,
    ReactantStoichCoeffs, ProductStoichCoeffs), m_TBconc(TBconc), m_NameTBconc(NameTBconc)
{}

ThreeBody::~ThreeBody() //Destructeur
{}





//---Falloff---

FalloffR::FalloffR(double A, double b, double E, string ID, string equation,
               string Etype, bool reversible, bool duplicate,
               vector<string> ReactantSpecies, vector<string> ProductSpecies,
               vector<double> ReactantStoichCoeffs, vector<double> ProductStoichCoeffs,
               /*ThreeBody*/
               vector<double> TBconc, vector<string> NameTBconc, 
               /*Falloff*/
               double A_low, double b_low, double E_low, string Etype_low
               ) //Constructeur
   :ThreeBody(A, b, E, ID, equation, Etype, reversible, duplicate, ReactantSpecies, ProductSpecies,
    ReactantStoichCoeffs, ProductStoichCoeffs, TBconc, NameTBconc), 
    m_A_low(A_low), m_b_low(b_low), m_E_low(E_low), m_Etype_low(Etype_low)
{}

FalloffR::~FalloffR() //Destructeur
{}





//---Lindemann---

Lindemann::Lindemann(double A, double b, double E, string ID, string equation,
               string Etype, bool reversible, bool duplicate,
               vector<string> ReactantSpecies, vector<string> ProductSpecies,
               vector<double> ReactantStoichCoeffs, vector<double> ProductStoichCoeffs,
               /*ThreeBody*/
               vector<double> TBconc, vector<string> NameTBconc, 
               /*Falloff*/
               double A_low, double b_low, double E_low, string Etype_low
               ) //Constructeur
   :FalloffR(A, b, E, ID, equation, Etype, reversible, duplicate, ReactantSpecies, ProductSpecies,
    ReactantStoichCoeffs, ProductStoichCoeffs, TBconc, NameTBconc,
    A_low, b_low, E_low, Etype_low) 
{}

Lindemann::~Lindemann() //Destructeur
{}





//---Troe---

Troe::Troe(double A, double b, double E, string ID, string equation,
               string Etype, bool reversible, bool duplicate,
               vector<string> ReactantSpecies, vector<string> ProductSpecies,
               vector<double> ReactantStoichCoeffs, vector<double> ProductStoichCoeffs,
               /*ThreeBody*/
               vector<double> TBconc, vector<string> NameTBconc, 
               /*Falloff*/
               double A_low, double b_low, double E_low, string Etype_low,
               /*Troe*/
               vector<double> TroeCoeffs
               ) //Constructeur
   :FalloffR(A, b, E, ID, equation, Etype, reversible, duplicate, ReactantSpecies, ProductSpecies,
    ReactantStoichCoeffs, ProductStoichCoeffs, TBconc, NameTBconc,
    A_low, b_low, E_low, Etype_low), m_TroeCoeffs(TroeCoeffs) 
{}

Troe::~Troe() //Destructeur
{}






