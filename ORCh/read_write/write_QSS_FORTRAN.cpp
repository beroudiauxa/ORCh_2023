#include <fstream>
#include <sstream>
#include "math.h"
#include <stdlib.h>
#include "write_QSS_FORTRAN.h"

#include "../tools/tools.h"

//---Write_QSS_FORTRAN---
Write_QSS_FORTRAN::Write_QSS_FORTRAN() //Constructeur
{}



//=================================================================
// <SUBROUTINE NAME= print_sign >" << endl;
// <DESCRIPTION> Routine to return the string of a number with its positive or negative sign </DESCRIPTION>" << endl;
string print_sign(double input)
{
   ostringstream output;
   output.precision(8);

   if (input > 0.0)
      output << "+ " << fabs(input);
   else
      output << "- " << fabs(input);

   string string_output = output.str();
   return string_output;
}
//=================================================================





int Get_species_number(string species, vector<Species_ORCh*> listSpecies) 
{
   int species_number = -1;

   for (unsigned int k=0; k<listSpecies.size(); k++)
   {
      if (listSpecies[k]->m_Name == species)
         species_number = k;
   }

   return species_number;
}

void Write_QSS_FORTRAN::Write_QSS_file_in_FORTRAN(string mech, string write_mech, vector<bool> &QSS_Species) const
{
   vector<Species_ORCh*> listSpecies;
   vector<Reaction_ORCh*> listReactions;

   Read *r = new Read();
   r->Read_species(mech, listSpecies);
   r->Read_reactions(mech, listReactions);

   int nsp = listSpecies.size();
   int nreac = listReactions.size();


   vector<int> Index_QSS_Species_number;
   int count_nQSS = 0;
   for (int k=0; k<nsp; k++)
   {
      if (QSS_Species[k] == false)
      {
         Index_QSS_Species_number.push_back(count_nQSS);
         count_nQSS += 1;
      }
      else
         Index_QSS_Species_number.push_back(-1);
   }

   for (int k=0; k<nsp; k++)
   {
      if (QSS_Species[k])
         listSpecies[k]->m_QSS = true;
   }

   int nQSS_Species = 0;
   for (int k=0; k<nsp; k++)
   {
      if (QSS_Species[k])
      {
         nQSS_Species += 1;
      }
   }



   ofstream write(write_mech.c_str(), ofstream::out);
   write.precision(7);


   write << endl;
   write << "      !=================================================================" << endl;
   write << "      !# <SUBROUTINE NAME=...>" << endl;
   write << "      !#  <DESCRIPTION> Routine to compute the reaction rates ... </DESCRIPTION>" << endl;
   write << "      !#  <VAR NAME=... TYPE=... IO=In/Out> ... </VAR>" << endl;
   write << "      !" << endl;
   write << "      subroutine source_terms(nnode,neqs,rhoinv,tempe,w_spec,wmol,source_spec)" << endl;
   write << endl;
   write << "      implicit none" << endl;
   write << endl;
   write << "         ! ----------------------" << endl;
   write << "         " << endl;
   write << "         integer, parameter :: pr=kind(1.0d0)" << endl;
   write << "         " << endl;
   write << "         ! IN/OUT" << endl;
   write << "         integer :: neqs,nnode" << endl;
   write << "         real(pr) :: rhoinv" << endl;
   write << "         real(pr) :: tempe" << endl;
   write << "         real(pr) :: w_spec(1:neqs)" << endl;
   write << "         real(pr) :: wmol(1:neqs)" << endl;
   write << "         real(pr) :: source_spec(1:neqs,1:nnode)" << endl;
   write << "         " << endl;
   write << "         ! ----------------------" << endl;
   write << endl;

   int nSimple = 0;
   int nFalloff = 0;
   int nThreeBody = 0;
   int nTroe = 0;
   int nLindemann = 0;

   for (unsigned int j=0; j<listReactions.size(); j++)
   {
      if (dynamic_cast <Simple *> (listReactions[j]))
         nSimple += 1;

      if (dynamic_cast <ThreeBody *> (listReactions[j]))
         nThreeBody += 1;

      if (dynamic_cast <FalloffR *> (listReactions[j]))
         nFalloff += 1;

      if (dynamic_cast <Troe *> (listReactions[j]))
         nTroe += 1;

      if (dynamic_cast <Lindemann *> (listReactions[j]))
         nLindemann += 1;
   }

   write << "         real(pr) :: C(1:neqs)" << endl;
   write << "         real(pr) :: CTB(1:" << nThreeBody << ")" << endl;
   write << "         real(pr) :: ROP(1:" << listReactions.size() << ")" << endl;
   write << "         real(pr) :: ROP_lowP(1:" << nFalloff << ")" << endl;
   //write << "         double precision :: Qr(1:" << listReactions.size() << ")" << endl;
   //write << "         double precision :: SMH(1:" << nsp << ") ! Computed from NASA coefficients" << endl;
   //write << "         double precision :: EG(1:" << nsp << ") ! Species equilibrium" << endl;
   //write << "         double precision :: EQ(1:"<< listReactions.size() << ") ! Reactions equilibrium" << endl;
   //write << "         double precision :: ROP(1:"<< listReactions.size() << ") ! Reactions rate of progress" << endl;


   write << "         integer :: i, j, ispec, ns, nr" << endl;
   //write << "         double precision :: RU, PATM, SMALL" << endl;
   write << "         real(pr) :: SMALL" << endl;
   write << "         real(pr) :: T, lnT, TI" << endl;
   write << "         real(pr) :: rhoi" << endl;
   //write << "         double precision :: Tn1, Tn2, Tn3, Tn4, Tn5" << endl;
   //write << "         double precision :: PFAC, PFAC2, PFAC3" << endl;
   write << "         real(pr) :: Ctot" << endl;
   write << "         real(pr) :: Pr_value, Pcor, PrLog, Fcent, FcLog" << endl;
   write << "         real(pr) :: Xn, CprLog, FLog, Fc" << endl;

   write << "         real(pr) :: X(" << nQSS_Species << ")" << endl;

   write << "         real(pr) :: ";
   for (int i=0; i<nQSS_Species; i++)
   {
      write << "RHS" << i+1; 
      if (i<nQSS_Species-1) 
         write << ", ";
      else
         write << endl;
   }

   write << "         real(pr) :: ";
   for (int i=0; i<nQSS_Species; i++)
   {
         write << "A" << i+1 << "(" << nQSS_Species << ")" ; 
         if (i<nQSS_Species-1) 
            write << ", ";
         else
            write << endl;
   }
   write << "         real(pr) :: Den,oorhoi,CTBtmp"<< endl;
   write << "         real(pr) :: coefficient" << endl;
   write << "         real(pr) :: sum_for_substitution" << endl;


   //write << endl;
   //write << "         ! chemistry properties" << endl;
   //write << "         nspecies = " << nsp << endl;
   //write << "         nreactions = " << nreac << endl;
   //write << endl;

   //write << "         ! mixture properties" << endl;
   //write << "         RU = 8.31447215E3" << endl;
   //write << "         PATM = 1.01325E5" << endl;
   write << "         SMALL = 1.0E-15" << endl;
   write << endl;




   write << endl;
   write << "         T = tempe" << endl;
   write << "         rhoi = rhoinv" << endl;
   write << "         oorhoi = 1.0_pr/rhoi" << endl;
   write << "         do ispec=1,neqs" << endl;
   write << "            C(ispec) = max( w_spec(ispec), 0.0_pr )" << endl;
   write << "            C(ispec) = min( C(ispec), oorhoi )" << endl;
   write << "            C(ispec) = C(ispec) / (wmol(ispec)*1.0e3_pr)" << endl;
   write << "         end do" << endl;

   write << endl;
   write << "         lnT = log(T)" << endl;
   write << "         TI = 1/T" << endl;
   write << endl;


   double R = 8.31447215;


   for (unsigned int j=0; j<listReactions.size(); j++)
   {
      write << "         ROP(" << j+1 << ") = ";

      //If non-zero activation energy
      if (listReactions[j]->m_E != 0.0)
      {
         write << "exp(" << log(listReactions[j]->m_A) << "_pr";
         if (listReactions[j]->m_E > 0.0)
            write << " - " << listReactions[j]->m_E*4.184/R << "_pr*TI";
         else
            write << " + " << fabs(listReactions[j]->m_E*4.184/R) << "_pr*TI";

         if (listReactions[j]->m_b != 0.0)
         {
            if (listReactions[j]->m_b > 0.0)
               write << " + " << listReactions[j]->m_b << "_pr*lnT";
            else
               write << " - " << fabs(listReactions[j]->m_b) << "_pr*lnT";
         }
         write << ")" << endl;
      }
      else
      {
         if (listReactions[j]->m_b != 0.0)
         {
            write << "exp(" << log(listReactions[j]->m_A) << "_pr";

            if (listReactions[j]->m_b > 0.0)
               write << " + " << listReactions[j]->m_b << "_pr*lnT)" << endl;
            else
               write << " - " << fabs(listReactions[j]->m_b) << "_pr*lnT)" << endl;
         }
         else
            write << listReactions[j]->m_A << "_pr"  << endl;
      }
   }

   write << endl;

   int count_falloff = 1;
   for (unsigned int j=0; j<listReactions.size(); j++)
   {
      if (dynamic_cast <FalloffR *> (listReactions[j]))
      {
         write << "         ROP_lowP(" << count_falloff << ") = ";

         //If non-zero activation energy
         if (dynamic_cast <FalloffR *> (listReactions[j])->m_E_low != 0.0)
         {
            write << "exp(" << log(dynamic_cast <FalloffR *> (listReactions[j])->m_A_low);
            if (dynamic_cast <FalloffR *> (listReactions[j])->m_E_low > 0.0)
               write << " - " << dynamic_cast <FalloffR *> (listReactions[j])->m_E_low*4.184/R << "*TI";
            else
               write << " + " << fabs(dynamic_cast <FalloffR *> (listReactions[j])->m_E_low*4.184/R) << "*TI";

            if (dynamic_cast <FalloffR *> (listReactions[j])->m_b_low != 0.0)
            {
               if (dynamic_cast <FalloffR *> (listReactions[j])->m_b_low > 0.0)
                  write << " + " << dynamic_cast <FalloffR *> (listReactions[j])->m_b_low << "*lnT";
               else
                  write << " - " << fabs(dynamic_cast <FalloffR *> (listReactions[j])->m_b_low) << "*lnT";
            }
            write << ")" << endl;
         }
         else
         {
            if (dynamic_cast <FalloffR *> (listReactions[j])->m_b_low != 0.0)
            {
               write << "exp(" << log(dynamic_cast <FalloffR *> (listReactions[j])->m_A_low);

               if (dynamic_cast <FalloffR *> (listReactions[j])->m_b_low > 0.0)
                  write << " + " << dynamic_cast <FalloffR *> (listReactions[j])->m_b_low << "*lnT)" << endl;
               else
                  write << " - " << fabs(dynamic_cast <FalloffR *> (listReactions[j])->m_b_low) << "*lnT)" << endl;
            }
            else
               write << dynamic_cast <FalloffR *> (listReactions[j])->m_A_low << endl;
         }

         count_falloff += 1;
      }

   }


/*
   //Estimation of the equilibrium through the Gibbs free energy calculation
   write << endl;
   write << "         ! thermal data " << endl << endl;
   write << "         Tn1 = lnT - 1.0 " << endl;
   write << "         Tn2 = T " << endl;
   write << "         Tn3 = Tn2 * T " << endl;
   write << "         Tn4 = Tn3 * T " << endl;
   write << "         Tn5 = Tn4 * T " << endl;

   //Check for the different Temperatures Tmax in the low T region
   vector<double> T;
   bool Check_new_T = true;

   T.push_back(listSpecies[0]->m_NASA_lowT_max);
   int nVarious_T = 1;

   for (int k=0; k<nsp; k++)
   {
      Check_new_T = true;
      for (int n=0; n<nVarious_T; n++)
      {
         if (listSpecies[k]->m_NASA_lowT_max == T[n])
         {
            Check_new_T = false;
         }
      }

      if (Check_new_T == true)
      {
         T.push_back(listSpecies[k]->m_NASA_lowT_max);
         nVarious_T += 1;
      }
   }


   write.precision(5);
   write.setf( std::ios::fixed, std:: ios::floatfield );
   write << endl << endl;
   write << "         ! low temperature Gibbs free energy" << endl;
   //Low temperature Gibbs values
   for (int n=0; n<nVarious_T; n++)
   {
      write << "         if (T < " << T[n] << ") then" << endl;

      for (int k=0; k<nsp; k++)
      {
         if (listSpecies[k]->m_NASA_lowT_max == T[n])
         {

            write << "            SMH(" << k+1 << ") = " << print_sign(listSpecies[k]->m_NASACoeffs_lowT[0]) << "*Tn1 &" << endl
                                              << "                    " << print_sign(0.5*listSpecies[k]->m_NASACoeffs_lowT[1]) << "*Tn2 &" << endl
                                              << "                    " << print_sign(listSpecies[k]->m_NASACoeffs_lowT[2]/6) << "*Tn3 &" << endl
                                              << "                    " << print_sign(listSpecies[k]->m_NASACoeffs_lowT[3]/12) << "*Tn4 &" << endl
                                              << "                    " << print_sign(0.05*listSpecies[k]->m_NASACoeffs_lowT[4]) << "*Tn5 &" << endl
                                              << "                    " << print_sign(-listSpecies[k]->m_NASACoeffs_lowT[5]) << "*TI &" << endl
                                              << "                    " << print_sign(listSpecies[k]->m_NASACoeffs_lowT[6]) << endl;
         }
      }
      write << "         end if " << endl;
   }

   write << endl << endl;
   write << "         ! high temperature Gibbs free energy" << endl;
   //High temperature Gibbs values
   for (int n=0; n<nVarious_T; n++)
   {
      write << "         if (T > " << T[n] << ") then" << endl;

      for (int k=0; k<nsp; k++)
      {
         if (listSpecies[k]->m_NASA_lowT_max == T[n])
         {

            write << "            SMH(" << k+1 << ") = " << print_sign(listSpecies[k]->m_NASACoeffs_highT[0]) << "*Tn1 &" << endl
                                              << "                    " << print_sign(0.5*listSpecies[k]->m_NASACoeffs_highT[1]) << "*Tn2 &" << endl
                                              << "                    " << print_sign(listSpecies[k]->m_NASACoeffs_highT[2]/6) << "*Tn3 &" << endl
                                              << "                    " << print_sign(listSpecies[k]->m_NASACoeffs_highT[3]/12) << "*Tn4 &" << endl
                                              << "                    " << print_sign(0.05*listSpecies[k]->m_NASACoeffs_highT[4]) << "*Tn5 &" << endl
                                              << "                    " << print_sign(-listSpecies[k]->m_NASACoeffs_highT[5]) << "*TI &" << endl
                                              << "                    " << print_sign(listSpecies[k]->m_NASACoeffs_highT[6]) << endl;
         }
      }
      write << "         end if" << endl;
   }


   write << endl << endl;

   write << "         do ns=1," << nsp << endl;
   write << "            EG(ns) = exp(SMH(ns)) " << endl;
   write << "         end do " << endl;

   write << endl;
   write << "         PFAC = PATM / (RU*T)" << endl;
   write << "         PFAC2 = PFAC*PFAC" << endl;
   write << "         PFAC3 = PFAC2*PFAC" << endl;

   bool Add_multiplication;

   //Equilibrium on the SimpleReactions
   write << endl << endl;
   write << "         ! equilibrium constants " << endl;

   for (unsigned int j=0; j<listReactions.size(); j++)
   {
      //If the reaction is not reversible, it is not necessary to estimate the equilibrium
      if (listReactions[j]->m_reversible)
      {
         write << "         EQ(" << j+1 <<") = ";

         Add_multiplication = false;
         for (unsigned int nj=0; nj<listReactions[j]->m_ProductSpecies.size(); nj++)
         {
            if (Get_species_number(listReactions[j]->m_ProductSpecies[nj], listSpecies) != -1)
            {
               if (Add_multiplication)
               {
                  write << "*";
               }
               write << "EG(" << Get_species_number(listReactions[j]->m_ProductSpecies[nj], listSpecies)+1 << ")";

               if (listReactions[j]->m_ProductStoichCoeffs[nj] == 2)
               {
                   write << "*EG(" << Get_species_number(listReactions[j]->m_ProductSpecies[nj], listSpecies)+1 << ")";
               }

               Add_multiplication = true;
            }
         }
         for (unsigned int nj=0; nj<listReactions[j]->m_ReactantSpecies.size(); nj++)
         {
            if (Get_species_number(listReactions[j]->m_ReactantSpecies[nj], listSpecies) != -1)
            {
               write << "/EG(" << Get_species_number(listReactions[j]->m_ReactantSpecies[nj], listSpecies)+1 << ")";

               if (listReactions[j]->m_ReactantStoichCoeffs[nj] == 2)
               {
                   write << "/EG(" << Get_species_number(listReactions[j]->m_ReactantSpecies[nj], listSpecies)+1 << ")";
               }
            }
         }

         double Delta_n = 0.0;
         for (unsigned int k=0; k<listReactions[j]->m_ReactantSpecies.size(); k++)
         {
            Delta_n -= listReactions[j]->m_ReactantStoichCoeffs[k];
         }
         for (unsigned int k=0; k<listReactions[j]->m_ProductSpecies.size(); k++)
         {
            Delta_n += listReactions[j]->m_ProductStoichCoeffs[k];
         }

         if (Delta_n == 1)
         {
            write << "*PFAC";
         }
         else if (Delta_n == -1)
         {
            write << "/PFAC";
         }
         else if (Delta_n != 0)
         {
            cout << "ERROR - Delta n is different from 1 or -1" << "  value: " << Delta_n << endl;
         }
         write << endl;
      }
   }



   write << endl << endl << endl;
   for (unsigned int j=0; j<listReactions.size(); j++)
   {
      if (listReactions[j]->m_reversible)
      {
         write << "         Qr(" << j+1 << ") = Qf(" << j+1 << ") / MAX(EQ(" << j+1 << "), SMALL)" << endl;
      }
      else
         write << "         Qr(" << j+1 << ") = 0.0" << endl;
   }
*/

   write << endl << endl << endl;
   write << "         ! three body concentration " << endl << endl;
   write << "         Ctot = sum(C(1:neqs))" << endl;


   write.precision(2);
   write << endl << endl;
   write << "         ! three body efficiencies" << endl;

   int count_threebody = 1;

   for (unsigned int j=0; j<listReactions.size(); j++)
   {
      if (dynamic_cast <ThreeBody *> (listReactions[j]))
      {
         write << "         CTB(" << count_threebody << ") = Ctot";
         for (unsigned int k=0; k<(dynamic_cast <ThreeBody *> (listReactions[j]))->m_TBconc.size(); k++)
         {

            if (Index_QSS_Species_number
                  [Get_species_number((dynamic_cast <ThreeBody *> (listReactions[j]))->
                                                  m_NameTBconc[k], listSpecies)] != -1)
            {
               if ((dynamic_cast <ThreeBody *> (listReactions[j]))->m_TBconc[k]-1 > 0.0)
                  write << " + " << (dynamic_cast <ThreeBody *> (listReactions[j]))->m_TBconc[k]-1;
               else
                  write << " - " << fabs((dynamic_cast <ThreeBody *> (listReactions[j]))->m_TBconc[k]-1);
               write << "*C(";
               write << Index_QSS_Species_number[Get_species_number((dynamic_cast <ThreeBody *> (listReactions[j]))->m_NameTBconc[k], listSpecies)]+1;
               write << ")";
               if (k%2 == 0 && k<(dynamic_cast <ThreeBody *> (listReactions[j]))->m_TBconc.size()-1)
               {
                  write << " & " << endl << "            ";
               }
            }

         }
         write << endl;
         count_threebody += 1;
      }
   }

   write << endl << endl;


   count_falloff = 1;
   count_threebody = 1;

   for (unsigned int j=0; j<listReactions.size(); j++)
   {
      if (dynamic_cast <ThreeBody *> (listReactions[j]))
      {
         if (dynamic_cast <Lindemann *> (listReactions[j]))
         {
            write << endl;
            write << "         ! lindemann type" << endl;
            write << "         Pr_value = ROP_lowP(" << count_falloff << ") * CTB(" << count_threebody <<")"
                              << " / ROP(" << j+1 << ")" << endl;
            write << "         Pcor = Pr_value / (1.0 + Pr_value)" << endl;
            write << "         ROP(" << j+1 << ") = ROP(" << j+1 << ") * Pcor" << endl;
//            if (listReactions[j]->m_reversible)
//               write << "         Qr(" << j+1 << ") = Qr(" << j+1 << ") * Pcor" << endl;

            count_falloff += 1;
         }

         else if (dynamic_cast <Troe *> (listReactions[j]))
         {
            write << endl;
            write << "         ! troe type" << endl;
            write << "         Pr_value = ROP_lowP(" << count_falloff << ") * CTB(" << count_threebody << ")"
                              << " / ROP(" << j+1 << ")" << endl;
            write << "         Pcor = Pr_value / (1.0 + Pr_value)" << endl;
            write << "         PrLog = log10(MAX(Pr_value, SMALL))" << endl;
            write << "         Fcent = " << (1-(dynamic_cast <Troe *> (listReactions[j]))->m_TroeCoeffs[0])
                         << "*exp(-T/" << (dynamic_cast <Troe *> (listReactions[j]))->m_TroeCoeffs[1] << ")+"
                         << (dynamic_cast <Troe *> (listReactions[j]))->m_TroeCoeffs[0];
            write << "*exp(-T/" << (dynamic_cast <Troe *> (listReactions[j]))->m_TroeCoeffs[2] << ") &" << endl;
            write << "                 +exp(" << -(dynamic_cast <Troe *> (listReactions[j]))->m_TroeCoeffs[3] << "/T)" << endl;
            write << "         FcLog = log10(MAX(Fcent, SMALL))" << endl;
            write << "         Xn = 0.75 - 1.27*FcLog" << endl;
            write << "         CprLog = PrLog - (0.4 + 0.67*FcLog)" << endl;
            write << "         FLog = FcLog/(1.0 + (CprLog/(Xn-0.14*CprLog))**2)" << endl;
            write << "         Fc = 10**FLog" << endl;
            write << "         Pcor = Pcor*Fc" << endl;
            write << "         ROP(" << j+1 << ") = ROP(" << j+1 << ") * Pcor" << endl;
//            if (listReactions[j]->m_reversible)
//               write << "         Qr(" << j+1 << ") = Qr(" << j+1 << ") * Pcor" << endl;

            count_falloff += 1;
         }
         count_threebody += 1;
      }
   }

   write << endl;

//Multiply by the concentration of the known species (not the one of the QSS species)

   count_threebody = 1;

   for (unsigned int j=0; j<listReactions.size(); j++)
   {
      write << "         ROP(" << j+1 << ") = ROP(" << j+1 << ")";

      for (unsigned int k=0; k<listReactions[j]->m_ReactantSpecies.size(); k++)
      {
         if (QSS_Species[Get_species_number(listReactions[j]->m_ReactantSpecies[k], listSpecies)] != true)
         {
              write << "*C(" << Index_QSS_Species_number[Get_species_number(listReactions[j]->m_ReactantSpecies[k], listSpecies)]+1 << ")";

              if (listReactions[j]->m_ReactantStoichCoeffs[k] != 1.0)
              {
                if (listReactions[j]->m_ReactantStoichCoeffs[k] == 2.0)
                {
                   write << "*C(" << Index_QSS_Species_number[Get_species_number(listReactions[j]->m_ReactantSpecies[k], listSpecies)]+1 << ")";
                }
                else
                {
                   cout << "Value of the StoichCoeff " << listReactions[j]->m_ReactantStoichCoeffs[k] << endl;
                   cout << "Error" << endl;
                }
              }
          } //end if QSS species
      }

      if ((dynamic_cast <ThreeBody *> (listReactions[j])))
      {
         if (!(dynamic_cast <FalloffR *> (listReactions[j])))
         {
            write << " *CTB(" << count_threebody << ")";
         }
      }
      write << endl;

      if ((dynamic_cast <ThreeBody *> (listReactions[j])))
      {
         count_threebody += 1;
      }
   }

   write << endl << endl;

   count_threebody = 0;

/*
   for (unsigned int j=0; j<listReactions.size(); j++)
   {
      if (listReactions[j]->m_reversible)
      {
         write << "         Qr(" << j+1 << ") = Qr(" << j+1 << ")";

         for (unsigned int k=0; k<listReactions[j]->m_ProductSpecies.size(); k++)
         {
            if (QSS_Species[Get_species_number(listReactions[j]->m_ProductSpecies[k], listSpecies)] != true)
            {
                 write << "*C(" << Index_QSS_Species_number[Get_species_number(listReactions[j]->m_ProductSpecies[k], listSpecies)]+1 << ")";

                 if (listReactions[j]->m_ProductStoichCoeffs[k] != 1.0)
                 {
                   if (listReactions[j]->m_ProductStoichCoeffs[k] == 2.0)
                   {
                      write << "*C(" << Index_QSS_Species_number[Get_species_number(listReactions[j]->m_ProductSpecies[k], listSpecies)]+1 << ")";
                   }
                   else
                   {
                      cout << "Value of the StoichCoeff " << listReactions[j]->m_ProductStoichCoeffs[k] << endl;
                      cout << "Error" << endl;
                   }
                 }
             } //end if QSS species
         }

         if ((dynamic_cast <ThreeBody *> (listReactions[j])))
         {
            if ((dynamic_cast <FalloffR *> (listReactions[j])) == false)
            {
               write << " *CTB(" << count_threebody << ")";
            }
         }
         write << endl;
      }
      if ((dynamic_cast <ThreeBody *> (listReactions[j])))
      {
         count_threebody += 1;
      }
   }
*/











   //write << endl;
   //write << "         do nr=1," << listReactions.size() << endl;
   //write << "            if (Qf(nr) < 1E-15) then" << endl;
   //write << "               Qf(nr) = 1E-15" << endl;
   //write << "            end if" << endl;
   //write << "            if (Qr(nr) < 1E-15) then" << endl;
   //write << "               Qr(nr) = 1E-15;" << endl;
   //write << "            end if" << endl;
   //write << "         end do" << endl;

   vector<bool> Reaction_fwd_to_zero (listReactions.size(), false);
   vector<bool> Reaction_rev_to_zero (listReactions.size(), false);

   for (unsigned int j=0; j<listReactions.size(); j++)
   {
      for (unsigned int k=0; k<listReactions[j]->m_ReactantSpecies.size(); k++)
      {
         if (QSS_Species[Get_species_number(listReactions[j]->m_ReactantSpecies[k], listSpecies)])
         {
            if (listReactions[j]->m_ReactantStoichCoeffs[k] == 2)
            {
               Reaction_fwd_to_zero[j] = true;
            }

            for (unsigned int ksec=0; ksec<listReactions[j]->m_ReactantSpecies.size(); ksec++)
            {
               if (ksec != k)
               {
                  if (QSS_Species[Get_species_number(listReactions[j]->m_ReactantSpecies[ksec], listSpecies)])
                  {
                      Reaction_fwd_to_zero[j] = true;
                  }
               }
            } //end for ksecond
         }
      }

      for (unsigned int k=0; k<listReactions[j]->m_ProductSpecies.size(); k++)
      {
         if (QSS_Species[Get_species_number(listReactions[j]->m_ProductSpecies[k], listSpecies)])
         {
            if (listReactions[j]->m_ProductStoichCoeffs[k] == 2)
            {
               Reaction_rev_to_zero[j] = true;
            }

            for (unsigned int ksec=0; ksec<listReactions[j]->m_ProductSpecies.size(); ksec++)
            {
               if (ksec != k)
               {
                  if (QSS_Species[Get_species_number(listReactions[j]->m_ProductSpecies[ksec], listSpecies)])
                  {
                      Reaction_rev_to_zero[j] = true;
                  }
               }
            } //end for ksecond
         }
      }
   }

   for (unsigned int j=0; j<listReactions.size(); j++)
   {
      if (Reaction_fwd_to_zero[j])
         write << "         ROP(" << j+1 << ") = 0.0" << endl;
 //     if (Reaction_rev_to_zero[j])
 //        write << "         Qr(" << j+1 << ") = 0.0" << endl;
   }


 
   //Solve for the QSS species
   write << endl << endl << endl;
   write << "         ! solve for the QSS species " << endl;



if (nQSS_Species > 0)
{
   vector<vector<bool> > QSS_matrix (nQSS_Species, vector<bool> (nQSS_Species, false));

   for (int k=0; k<nQSS_Species; k++)
   {
         for (int kb=0; kb<nQSS_Species; kb++)
         {
            if (k == kb)
               QSS_matrix[k][kb] = true;
            else
               QSS_matrix[k][kb] = false;
         }
   }

       write << "         X(1:" << nQSS_Species << ") = 0.0_pr" << endl;
    write << endl;


   vector <int> QSSindex (nsp, -1);

   int i = 0;
   for (int k=0; k<nsp; k++)
   {
      if (QSS_Species[k])
      {
         QSSindex[k] = i;
         i += 1;
      }
   }


   string Den;
   int count_den;
   int count_RHS;
   bool look_on_the_other_side;
   bool write_into_RHS;

   vector<string> QSS_SpeciesCoeffs (nQSS_Species+1, "");


















   for (int k=0; k<nsp; k++)
   {
      for (int ka=0; ka<nQSS_Species+1; ka++)
         QSS_SpeciesCoeffs[ka] = "";

      Den = "";
      if (QSS_Species[k])
      {
         count_den = 1;
         count_RHS = 1;



         //WE SHOULD ADD A COUNT_DEN and COUNT_RHS TO AVOID AVING    "&"    AT THE END OF THE LINE











         for (unsigned int j=0; j<listReactions.size(); j++)
         {
            look_on_the_other_side = false;

            for (unsigned int ka=0; ka<listReactions[j]->m_ReactantSpecies.size(); ka++)
            {
               if (k == Get_species_number(listReactions[j]->m_ReactantSpecies[ka], listSpecies))
               {
                  stringstream sstm;
                  sstm << " + ROP(" << j+1 << ")";
                  if (count_den == 5)
                  {
                     sstm << " &" << endl << "          ";
                     count_den = 0;
                  }
                  Den.append(sstm.str());
                  count_den += 1;
                  look_on_the_other_side = true;
               }
            }
            
            if (look_on_the_other_side)
            {
               write_into_RHS = true;
               //Look for the other QSS species involved (which can only be present on the other side of the reaction)

               for (unsigned int ka=0; ka<listReactions[j]->m_ProductSpecies.size(); ka++)
               {
                  if (QSS_Species[Get_species_number(listReactions[j]->m_ProductSpecies[ka], listSpecies)])
                  {
                     //stringstream sstm;
                     //sstm << " + Qr(" << j+1 << ")";

                     //QSS_SpeciesCoeffs[QSSindex[Get_species_number(listReactions[j]->m_ProductSpecies[ka], listSpecies)]]
                     //                 .append(sstm.str());

                     QSS_matrix[QSSindex[k]]
                               [QSSindex[Get_species_number(listReactions[j]->m_ProductSpecies[ka], listSpecies)]] = true;
                     write_into_RHS = false;

                  }
               }

               //if (write_into_RHS == true)
               //{
               //       stringstream sstm;
               //       sstm << " + Qr(" << j+1 << ")";
               //       if (count_RHS == 5)
               //       {
               //          sstm << " &" << endl << "          ";
               //          count_RHS = 0;
               //       }
               //       QSS_SpeciesCoeffs[nQSS_Species].append(sstm.str());
               //       count_RHS += 1;
               //}
               
            } //end if (look_on_the_other_side)

      //---Repeat this for a QSS species into the products

            look_on_the_other_side = false;
            for (unsigned int ka=0; ka<listReactions[j]->m_ProductSpecies.size(); ka++)
            {
               if (k == Get_species_number(listReactions[j]->m_ProductSpecies[ka], listSpecies))
               {
               //   stringstream sstm;
               //   sstm << " + Qr(" << j+1 << ")";
               //   if (count_den == 5)
               //   {
               //      sstm << " &" << endl << "          ";
               //      count_den = 0;
               //   }
               //   Den.append(sstm.str());
               //   count_den += 1;
                  look_on_the_other_side = true;
               }
            }

            if (look_on_the_other_side)
            {
               write_into_RHS = true;
               //Look for the other QSS species involved (which can only be present on the other side of the reaction)

               for (unsigned int ka=0; ka<listReactions[j]->m_ReactantSpecies.size(); ka++)
               {
                  if (QSS_Species[Get_species_number(listReactions[j]->m_ReactantSpecies[ka], listSpecies)])
                  {
                     stringstream sstm;
                     sstm << " + ROP (" << j+1 << ")";

                     QSS_SpeciesCoeffs[QSSindex[Get_species_number(listReactions[j]->m_ReactantSpecies[ka], listSpecies)]]
                                      .append(sstm.str());

                     QSS_matrix[QSSindex[k]]
                               [QSSindex[Get_species_number(listReactions[j]->m_ReactantSpecies[ka], listSpecies)]] = true;
                     write_into_RHS = false;
                  }
               }
               if (write_into_RHS == true)
               {
                      stringstream sstm;
                      sstm << " + ROP(" << j+1 << ")";
                      if (count_RHS == 5)
                      {
                         sstm << " &" << endl << "          ";
                         count_RHS = 0;
                      }
                      QSS_SpeciesCoeffs[nQSS_Species].append(sstm.str());
                      count_RHS += 1;
               }
            } //end if (look_on_the_other_side)
         } //end for (unsigned int j=0; j<listReactions.size(); j++)





         size_t remove_first_plus = Den.find("+");
         Den.erase(0, remove_first_plus+1);

         write << endl << endl << "         ! QSS species " << listSpecies[k]->m_Name << endl;

         if (Den != "")
         {
            write << "         Den = " << Den << endl;
            write << "         Den = MAX(Den, 1.0)" << endl;
         }
         else
            //should be removed to something else
            write << "         Den = 0.0" << endl;

            write << "         A"<< QSSindex[k]+1 <<"(1:" << nQSS_Species << ") = 0.0_pr"  << endl;
         for (int ka=0; ka<nQSS_Species; ka++)
         {
            if (QSS_SpeciesCoeffs[ka] != "")
            {
               remove_first_plus = QSS_SpeciesCoeffs[ka].find("+");
               QSS_SpeciesCoeffs[ka].erase(0, remove_first_plus+1);
               write << "         A" << QSSindex[k]+1 << "(" << ka+1 << ") = ";
               write << "(" << QSS_SpeciesCoeffs[ka] << ") /Den" << endl;
            }
     //       else
     //          write << "         A" << QSSindex[k]+1 << "_" << ka+1 << " = 0.0" << endl;

         }
         write << endl;


         remove_first_plus = QSS_SpeciesCoeffs[nQSS_Species].find("+");
         QSS_SpeciesCoeffs[nQSS_Species].erase(0, remove_first_plus+1);

         if (QSS_SpeciesCoeffs[nQSS_Species] != "")
         {
            write << "         RHS" << QSSindex[k]+1 << " = ";
            write << "-(" << QSS_SpeciesCoeffs[nQSS_Species] << ") /Den" << endl;
         }
         else
         {
            write << "         RHS" << QSSindex[k]+1 << " = 0.0";
         }

      } //end if QSS_Species[k]
   } //end for k=0; k<nsp; k++

   write << endl << endl;

   for (int k=0; k<nQSS_Species; k++)
   {
      write << "         A" << k+1 << "(" << k+1 << ") = -1" << endl;
   }

   write << endl << endl << endl;
   write << "         ! solve the matrix" << endl;

   write << endl << endl;

   vector<double> b (nQSS_Species);
   vector<vector<double> > A (nQSS_Species,vector<double> (nQSS_Species));

   for (int i=0; i<nQSS_Species; i++)
   {
      b[i] = 2.0*(i+1); //Give ordinary values to the b vector
      for (int j=0; j<nQSS_Species; j++)
      {
         if (QSS_matrix[i][j] == true)
         {
            double rand = 0.0;
            give_random_number_from_0_to_1(rand);
            A[i][j] = 100*rand;
         }
      }
   }

   vector<double> x (nQSS_Species);

   //Elimination
   write << endl << "         ! //-----Elimination-----" << endl;

    for(int i=0; i<nQSS_Species-1; i++)
    {
        for(int j=i+1; j<nQSS_Species; j++)
        {
            double Coeff = -(A[j][i]/A[i][i]);
            bool storeCoeff = true;

            A[j][i] = 0.0;
            for(int k = i+1; k < nQSS_Species; k++)
            {
                A[j][k] = A[j][k] + Coeff*A[i][k];

                if (Coeff != 0.0)
                {
                   if (A[i][k] != 0.0)
                   {
                      if (storeCoeff)
                      {
                         write << "         coefficient = -(A" << j+1 << "_" << i+1 << "/A" << i+1 << "_" << i+1 << ")" << endl;
                         storeCoeff = false;
                      }
                      write << "         A" << j+1 << "_" << k+1 << " = " << "A" << j+1 << "_" << k+1 << " + coefficient*A" << i+1 << "_" << k+1 << endl;
                   }
                }
            }

            b[j] = b[j] + Coeff*b[i];

            if (storeCoeff == false)
               write << "         RHS" << j+1 << " = RHS" << j+1 << " + coefficient*RHS" << i+1 << endl << endl;
        }
    }



    //Back-substitution
    write << endl << "         ! //-----Back substitution-----//" << endl;

    x[nQSS_Species-1] = b[nQSS_Species-1]/A[nQSS_Species-1][nQSS_Species-1];
    write << "         X" << nQSS_Species << " = "
          << "RHS" << nQSS_Species << " /A" << nQSS_Species << "_" << nQSS_Species << endl;

    for(int i = nQSS_Species-2; i>=0; i--)
    {
        stringstream sstm;
        sstm << "         sum_for_substitution = ";
        bool Coeff = false;

        double sum = 0;
        for(int j = i; j < nQSS_Species; j++)
        {
            sum += x[j]*A[i][j];

           if ((x[j]*A[i][j]) != 0.0)
           {
              sstm << "+ X" << j+1 << "*A" << i+1 << "_" << j+1 << "";
              Coeff = true;
           }

        }

        if (Coeff)
        {
          write << endl << sstm.str() << endl;
          write << "         X" << i+1 << " = (RHS" << i+1 << "-sum_for_substitution)/A" << i+1 << "_" << i+1 << endl;
        }
        else
        {
          write << "         x" << i+1 << " = RHS" << i+1 << "/A" << i+1 << "_" << i+1 << endl;
        }

        x[i] = (b[i]-sum)/A[i][i];
    }


   write << endl << endl << endl << endl;
   write << "         ! multiply the forward rates by the concentration of the QSS species " << endl;
   write << endl << endl;
   for (unsigned int j=0; j<listReactions.size(); j++)
   {
       for (int k=0; k<nsp; k++)
       {
          if (QSS_Species[k])
          {
             for (unsigned int ka=0; ka<listReactions[j]->m_ReactantSpecies.size(); ka++)
             {
                if (k == Get_species_number(listReactions[j]->m_ReactantSpecies[ka], listSpecies))
                {
                   //Within the forward part of this Simple reaction there is a QSS species
                   write << "         ROP(" << j+1 << ") = ROP("<< j+1 <<  ") * X" << QSSindex[k]+1 << endl;
                }
             }

            //for (unsigned int ka=0; ka<listReactions[j]->m_ProductSpecies.size(); ka++)
            //{
            //    if (k == Get_species_number(listReactions[j]->m_ProductSpecies[ka], listSpecies))
            //    {
            //       //Within the forward part of this Simple reaction there is a QSS species
            //       write << "         Qr(" << j+1 << ") = Qr("<< j+1 << ") * X" << QSSindex[k]+1 << endl;
            //    }
            // }
          }
       }
   }

} //End if (nQSS_Species > 0) 
else
{
   write << endl << "         ! //---0 QSS expressions---//" << endl;
}

  //write << endl;
  //for (unsigned int j=0; j<listReactions.size(); j++)
  //{
  //   if (listReactions[j]->m_reversible)
  //      write << "         ROP(" << j+1 << ") = ROPQf(" << j+1 << ") - Qr(" << j+1 << ")" << endl;
  //   else
  //      write << "         ROP(" << j+1 << ") = Qf(" << j+1 << ")" << endl;

  //}



  write << endl << endl << endl;
  for (int k=0; k<nsp; k++)
  {
     if (QSS_Species[k] == false)
     {
        write << "         ! species " << Index_QSS_Species_number[k]+1 << endl;
        write << "         source_spec(" << Index_QSS_Species_number[k]+1 << ") =";

        int count_total_element = 0;
        for (unsigned int j=0; j<listReactions.size(); j++)
        {
           for (unsigned int ka=0; ka<listReactions[j]->m_ReactantSpecies.size(); ka++)
           {
              if (k == Get_species_number(listReactions[j]->m_ReactantSpecies[ka], listSpecies))
              {
                 if (listReactions[j]->m_ReactantStoichCoeffs[ka] == 2)
                    count_total_element += 1;
                 else if (listReactions[j]->m_ReactantStoichCoeffs[ka] == 1)
                    count_total_element += 1;
              }
           }
          for (unsigned int ka=0; ka<listReactions[j]->m_ProductSpecies.size(); ka++)
           {
              if (k == Get_species_number(listReactions[j]->m_ProductSpecies[ka], listSpecies))
              {
                 if (listReactions[j]->m_ProductStoichCoeffs[ka] == 2)
                    count_total_element += 1;
                 else if (listReactions[j]->m_ProductStoichCoeffs[ka] == 1)
                    count_total_element += 1;
              }
           }
        }


        int count_element = 0;

        for (unsigned int j=0; j<listReactions.size(); j++)
        {
           for (unsigned int ka=0; ka<listReactions[j]->m_ReactantSpecies.size(); ka++)
           {
              if (k == Get_species_number(listReactions[j]->m_ReactantSpecies[ka], listSpecies))
              {
                 if (listReactions[j]->m_ReactantStoichCoeffs[ka] == 2)
                 {
                    write << "-2*ROP(" << j+1 << ") ";
                    count_element += 1;
                 }
                 else if (listReactions[j]->m_ReactantStoichCoeffs[ka] == 1)
                 {
                    write << "-ROP(" << j+1 << ") ";
                    count_element += 1;
                 }
                 else
                    cout << "ERROR stoichiometric value different from 1 or 2 " << endl;

                 if (count_element > 2 && (count_element != count_total_element))
                 {
                    count_element = 0;               
                    write << " &" << endl << "                   ";
                 }

              }
           }

           for (unsigned int ka=0; ka<listReactions[j]->m_ProductSpecies.size(); ka++)
           {
              if (k == Get_species_number(listReactions[j]->m_ProductSpecies[ka], listSpecies))
              {
                 if (listReactions[j]->m_ProductStoichCoeffs[ka] == 2)
                 {
                    write << "+2*ROP(" << j+1 << ") ";
                    count_element += 1;
                 }
                 else if (listReactions[j]->m_ProductStoichCoeffs[ka] == 1)
                 {
                    write << "+ROP(" << j+1 << ") ";
                    count_element += 1;
                 }
                 else
                    cout << "ERROR stoichiometric value different from 1 or 2 " << endl;

                 if (count_element > 2 && (count_element != count_total_element))
                 {
                    count_element = 0;
                    write << " &" << endl << "                   ";
                 }
              }
           }
        }

        write << endl << endl;
     }
  }











   write << endl;
   write << endl;
   write << "         do ns=1,neqs" << endl;
   write << "            source_spec(ns) = source_spec(ns)*wmol(ns)*1e+3_pr" << endl;
   write << "         end do" << endl;
   write << "" << endl;
   write << endl;
   write << "      end subroutine source_terms" << endl;
   write << "      !# </SUBROUTINE>" << endl;
   write << "      !=================================================================" << endl;




   write.close();































}

Write_QSS_FORTRAN::~Write_QSS_FORTRAN() //Destructeur
{}

