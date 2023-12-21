#ifndef Reaction_H
#define Reaction_H

#include <iostream>
#include <vector>

using namespace std;


//--------------------
class Reaction_ORCh
{
   public:
   //constructeur
   Reaction_ORCh(double A, double b, double E, string ID, string equation, 
            string Etype, bool reversible, bool duplicate,
            vector<string> ReactantSpecies, vector<string> ProductSpecies,
            vector<double> ReactantStoichCoeffs, vector<double> ProductStoichCoeffs); 

   //destructeur
   virtual ~Reaction_ORCh();

   public:
   double m_A;
   double m_b;
   double m_E;

   string m_ID;
   string m_equation;
   string m_Etype;

   bool m_reversible;
   bool m_duplicate;

   vector<string> m_ReactantSpecies;
   vector<string> m_ProductSpecies;
   vector<double> m_ReactantStoichCoeffs;
   vector<double> m_ProductStoichCoeffs;

};



class Simple : public Reaction_ORCh
{
   public:
   //constructeur
   Simple(double A, double b, double E, string ID, string equation,
          string Etype, bool reversible, bool duplicate,
          vector<string> ReactantSpecies, vector<string> ProductSpecies,
          vector<double> ReactantStoichCoeffs, vector<double> ProductStoichCoeffs); 

   //destructeur
   virtual ~Simple();
   
   private:

};


class ThreeBody : public Reaction_ORCh
{
   public:
   //constructeur
   ThreeBody(double A, double b, double E, string ID, string equation,
             string Etype, bool reversible, bool duplicate,
             vector<string> ReactantSpecies, vector<string> ProductSpecies,
             vector<double> ReactantStoichCoeffs, vector<double> ProductStoichCoeffs,
             /*From ThreeBody*/ 
             vector<double> TBconc, vector<string> NameTBconc); 

   //destructeur
   virtual ~ThreeBody();
   
   public:
   vector<double> m_TBconc;
   vector<string> m_NameTBconc;


};



//---Falloff is declared as a specific ThreeBody reaction---

class FalloffR : public ThreeBody 
{
   public:
   //constructeur
   FalloffR(double A, double b, double E, string ID, string equation,
           string Etype, bool reversible, bool duplicate,
           vector<string> ReactantSpecies, vector<string> ProductSpecies,
           vector<double> ReactantStoichCoeffs, vector<double> ProductStoichCoeffs,
           /*From ThreeBody*/
           vector<double> TBconc, vector<string> NameTBconc,
           /*From Falloff*/
           double A_low, double b_low, double E_low, string m_Etype_low
           );

   //destructeur
   virtual ~FalloffR();

   public:
   double m_A_low;
   double m_b_low;
   double m_E_low;
 
   string m_Etype_low;
};


class Lindemann : public FalloffR 
{
   public:
   //constructeur
   Lindemann(double A, double b, double E, string ID, string equation,
           string Etype, bool reversible, bool duplicate,
           vector<string> ReactantSpecies, vector<string> ProductSpecies,
           vector<double> ReactantStoichCoeffs, vector<double> ProductStoichCoeffs,
           /*From ThreeBody*/
           vector<double> TBconc, vector<string> NameTBconc,
           /*From Falloff*/
           double A_low, double b_low, double E_low, string m_Etype_low
           );

   //destructeur
   virtual ~Lindemann();

};


class Troe : public FalloffR 
{
   public:
   //constructeur
   Troe(double A, double b, double E, string ID, string equation,
           string Etype, bool reversible, bool duplicate,
           vector<string> ReactantSpecies, vector<string> ProductSpecies,
           vector<double> ReactantStoichCoeffs, vector<double> ProductStoichCoeffs,
           /*From ThreeBody*/
           vector<double> TBconc, vector<string> NameTBconc,
           /*From Falloff*/
           double A_low, double b_low, double E_low, string m_Etype_low,
           /*From Troe*/
           vector<double> TroeCoeffs
           );

   //destructeur
   virtual ~Troe();

   public:
   vector<double> m_TroeCoeffs;

};


#endif


