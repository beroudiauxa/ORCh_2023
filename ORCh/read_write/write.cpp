#include <fstream>

#include "write.h"
#include <math.h>

//---Write---

Write::Write() //Constructeur
{}



double round_to_digits(double value, int digits)
{

   double fractpart, intpart;
   fractpart = modf(value*pow(10, digits), &intpart);

   if (fractpart > 0.5)
      return (intpart+1)/pow(10, digits);
   else
      return intpart/pow(10, digits);

}



void Write::Write_xml_file(string mech, string mech_desc, string write_mech,
                           vector<bool> Species_to_add, vector<bool> Reactions_to_add, bool optim, vector<double> A_Arrh, vector<double> b_Arrh, vector<double> E_Arrh, 
                           vector<Species_ORCh*> listSpecies_lumping, vector<Reaction_ORCh*> listReactions_lumping) const
{




   vector<Species_ORCh*> listSpecies;
   vector<Reaction_ORCh*> listReactions;

   Read *r = new Read();
   r->Read_species(mech, listSpecies);
   r->Read_reactions(mech, listReactions);



   vector<bool> Reactions_fwd_to_add (listReactions.size(), false);
   vector<bool> Reactions_rev_to_add (listReactions.size(), false);


   for (int j=0; j<listReactions.size(); j++)
   {
      if (Reactions_to_add[j])
      {
         Reactions_fwd_to_add[j] = true;
         Reactions_rev_to_add[j] = true;
      }
   }

   ofstream write_file(write_mech.c_str(), ofstream::out);


   if (optim)
   {
      for (unsigned int j=0; j<listReactions.size(); j++)
      {
         listReactions[j]->m_A = A_Arrh[j];
         listReactions[j]->m_b = b_Arrh[j];
         listReactions[j]->m_E = E_Arrh[j];
      }

   }


   if (listSpecies_lumping.size() > 0)
   {
      for (unsigned int k=0; k<listSpecies.size(); k++)
      {
         listSpecies[k]->m_Name = listSpecies_lumping[k]->m_Name;
         for (int p=0; p<7; p++)
         {
            listSpecies[k]->m_NASACoeffs_highT[p] = listSpecies_lumping[k]->m_NASACoeffs_highT[p];
            listSpecies[k]->m_NASACoeffs_lowT[p] = listSpecies_lumping[k]->m_NASACoeffs_lowT[p];
         }
         listSpecies[k]->m_dipoleMoment = listSpecies_lumping[k]->m_dipoleMoment;
         listSpecies[k]->m_polarizability = listSpecies_lumping[k]->m_polarizability;

         listSpecies[k]->m_rot_relax = listSpecies_lumping[k]->m_rot_relax;
         listSpecies[k]->m_LJ_diameter = listSpecies_lumping[k]->m_LJ_diameter;
         listSpecies[k]->m_LJ_welldepth = listSpecies_lumping[k]->m_LJ_welldepth;

      }



      for (unsigned int j=0; j<listReactions.size(); j++)
      {
         listReactions[j]->m_A = listReactions_lumping[j]->m_A;
         if (dynamic_cast <FalloffR *> (listReactions[j]))
         {
            (dynamic_cast <FalloffR *> (listReactions[j]))->m_A_low = (dynamic_cast <FalloffR *> (listReactions_lumping[j]))->m_A_low;
         }

         listReactions[j]->m_equation = listReactions_lumping[j]->m_equation;
         for (int k=0; k<listReactions[j]->m_ReactantSpecies.size(); k++)
         {
            listReactions[j]->m_ReactantSpecies[k] = listReactions_lumping[j]->m_ReactantSpecies[k];
         }

         for (int k=0; k<listReactions[j]->m_ProductSpecies.size(); k++)
         {
            listReactions[j]->m_ProductSpecies[k] = listReactions_lumping[j]->m_ProductSpecies[k];
         }

         if (listReactions_lumping[j]->m_duplicate == true)
            listReactions[j]->m_duplicate = true;


         if (listReactions_lumping[j]->m_equation == "")
         {
            Reactions_fwd_to_add[j] = false;
            Reactions_rev_to_add[j] = false;
         }
      }
   }





   //Add the phase 
   write_file << "<?xml version=\"1.0\"?>" << endl;
   write_file << "<ctml>" << endl;
   write_file << "  <validate reactions=\"yes\" species=\"yes\"/>" << endl << endl;

   write_file << "  <!-- phase " << mech_desc << "     -->" << endl;
   write_file << "  <phase dim=\"3\" id=\"" << mech_desc << "\">" << endl;

   write_file << "    <elementArray datasrc=\"elements.xml\">O  H  C";

   for (unsigned int k=0; k<listSpecies.size(); k++)
   { 
      if (Species_to_add[k])
      {
         if (listSpecies[k]->m_Name == "N2")
            write_file << "  N";

         if (listSpecies[k]->m_Name == "AR")
            write_file << "  Ar";
      }
   }
   write_file << "</elementArray>" << endl;
   write_file << "    <speciesArray datasrc=\"#species_data\">" << endl;
   write_file << "    ";

   int countSpecies = 0;
   for (unsigned int k=0; k<listSpecies.size(); k++)
   {
      if (Species_to_add[k])
      {
         write_file << listSpecies[k]->m_Name << "  ";

         if ((countSpecies+1)%10 == 0)
            write_file << endl << "    ";

         countSpecies += 1;
      }
   }

   write_file << "</speciesArray>" << endl;
   write_file << "    <reactionArray datasrc=\"#reaction_data\"/>" << endl;
   write_file << "    <state>" << endl;
   write_file << "      <temperature units=\"K\">300.0</temperature>" << endl;
   write_file << "      <pressure units=\"Pa\">101325.0</pressure>" << endl;
   write_file << "    </state>" << endl;
   write_file << "    <thermo model=\"IdealGas\"/>" << endl;
   write_file << "    <kinetics model=\"GasKinetics\"/>" << endl;
   write_file << "    <transport model=\"Mix\"/>" << endl;
   write_file << "  </phase>" << endl << endl;

   //Add the chemical species involved
   write_file << "  <!-- species definitions     -->" << endl;
   write_file << "  <speciesData id=\"species_data\">" << endl << endl;



   //Write the species description
   //
   //1st species Name
   //2nd species NASA coefficients...

   for (unsigned int k=0; k<listSpecies.size(); k++)
   {
      if (Species_to_add[k])
      {
         write_file << "   <!--species " << listSpecies[k]->m_Name << "   -->" << endl;
         write_file << "   <species name=\"" << listSpecies[k]->m_Name << "\">" << endl;
         write_file << "     <atomArray>";

         if (listSpecies[k]->m_C != 0)
            write_file << "C:" << listSpecies[k]->m_C << " ";
         if (listSpecies[k]->m_H != 0)
            write_file << "H:" << listSpecies[k]->m_H << " ";
         if (listSpecies[k]->m_O != 0)
            write_file << "O:" << listSpecies[k]->m_O << " ";
         if (listSpecies[k]->m_N != 0)
            write_file << "N:" << listSpecies[k]->m_N << " ";
         if (listSpecies[k]->m_Ar != 0)
            write_file << "Ar:" << listSpecies[k]->m_Ar << " ";

         write_file << "</atomArray>" << endl;
         write_file << "     <thermo>" << endl;
         write_file << "       <NASA Tmax=\"" << listSpecies[k]->m_NASA_lowT_max << "\" Tmin=\"" << listSpecies[k]->m_NASA_lowT_min << "\" P0=\"100000.0\">" << endl;
         write_file << "         <floatArray name=\"coeffs\" size=\"7\">" << endl << "           ";
         write_file.precision(9);
         write_file.fixed;
         for (int c=0; c<7; c++)
         {
            if (c > 0)
               write_file << ",  ";
            if (c == 4)
               write_file << endl << "           ";
            write_file << listSpecies[k]->m_NASACoeffs_lowT[c];
         }
         write_file << "</floatArray>" << endl;
         write_file << "       </NASA>" << endl;
         write_file << "       <NASA Tmax=\"" << listSpecies[k]->m_NASA_highT_max << "\" Tmin=\"" << listSpecies[k]->m_NASA_highT_min << "\" P0=\"100000.0\">" << endl;
         write_file << "         <floatArray name=\"coeffs\" size=\"7\">" << endl << "           ";
         for (int c=0; c<7; c++)
         {
            if (c > 0)
               write_file << ",  ";
            if (c == 4)
               write_file << endl << "           ";
            write_file << listSpecies[k]->m_NASACoeffs_highT[c];
         }
         write_file.precision(5);
         write_file.scientific;
         write_file << "</floatArray>" << endl;
         write_file << "       </NASA>" << endl;
         write_file << "     </thermo>" << endl;
         write_file << "     <transport model=\"gas_transport\">" << endl;
         write_file << "       <string title=\"geometry\">" << listSpecies[k]->m_geometry << "</string>" << endl;
         write_file << "       <LJ_welldepth units=\"K\">" << listSpecies[k]->m_LJ_welldepth << "</LJ_welldepth>" << endl;
         write_file << "       <LJ_diameter units=\"A\">" << listSpecies[k]->m_LJ_diameter << "</LJ_diameter>" << endl;
         write_file << "       <dipoleMoment units=\"Debye\">" << listSpecies[k]->m_dipoleMoment << "</dipoleMoment>" << endl;
         write_file << "       <polarizability units=\"A3\">" << listSpecies[k]->m_polarizability << "</polarizability>" << endl;
         write_file << "       <rotRelax>" << listSpecies[k]->m_rot_relax << "</rotRelax>" << endl;
         write_file << "     </transport>" << endl;
         write_file << "   </species>" << endl << endl;
      }



         //write_file << "    " << listSpecies[k]->m_Description << endl << endl;
   }
   write_file << "  </speciesData>" << endl << endl;















   //Add the reactions
   write_file << "  <reactionData id=\"reaction_data\">" << endl << endl;
   int countReactions = 0;



   for (unsigned int j=0; j<listReactions.size(); j++)
   {
      if (Reactions_rev_to_add[j] == false)
         listReactions[j]->m_reversible = false;

      if (Reactions_fwd_to_add[j])
      {
         bool Print = true;
         for (unsigned int k=0; k<listSpecies.size(); k++)
         {
            for (unsigned int kr=0; kr<listReactions[j]->m_ReactantSpecies.size(); kr++)
            {
               if ((listSpecies[k]->m_Name == listReactions[j]->m_ReactantSpecies[kr])
                    && Species_to_add[k] == false)
               {
                  Print = false;
               }
            }

            for (unsigned int kp=0; kp<listReactions[j]->m_ProductSpecies.size(); kp++)
            {
               if ((listSpecies[k]->m_Name == listReactions[j]->m_ProductSpecies[kp])
                    && Species_to_add[k] == false)
               {
                  Print = false;
               }
            }
         }

         if (Print)
         {
            countReactions += 1;
            Write_reaction(write_file, listReactions[j], Species_to_add, listSpecies, countReactions, listReactions, Reactions_fwd_to_add, Reactions_rev_to_add);
         }
      }
   }

   write_file << "  </reactionData>" << endl << endl;
   write_file << "</ctml>" << endl << endl;





}










void Write::Write_xml_for_Analytic_Applications(string mech, string mech_desc, 
                    string write_mech, vector<bool> Species_to_add) const
{

   cout << "Analytic scheme 2/2 : " <<  write_mech << endl;
   cout << endl << endl;

   vector<Species_ORCh*> listSpecies;
   vector<Reaction_ORCh*> listReactions;

   Read *r = new Read();
   r->Read_species(mech, listSpecies);
   r->Read_reactions(mech, listReactions);

   ofstream write_file(write_mech.c_str(), ofstream::out);

   //Add the phase 
   write_file << "<?xml version=\"1.0\"?>" << endl;
   write_file << "<ctml>" << endl;
   write_file << "  <validate reactions=\"yes\" species=\"yes\"/>" << endl << endl;

   write_file << "  <!-- phase " << mech_desc << "     -->" << endl;
   write_file << "  <phase dim=\"3\" id=\"" << mech_desc << "\">" << endl;

   write_file << "    <elementArray datasrc=\"elements.xml\">O  H  C";

   for (unsigned int k=0; k<listSpecies.size(); k++)
   {
      if (Species_to_add[k])
      {
         if (listSpecies[k]->m_Name == "N2")
            write_file << "  N";

         if (listSpecies[k]->m_Name == "AR")
            write_file << "  Ar";
      }
   }
   write_file << "</elementArray>" << endl;
   write_file << "    <speciesArray datasrc=\"#species_data\">" << endl;
   write_file << "    ";

   int countSpecies = 0;
   for (unsigned int k=0; k<listSpecies.size(); k++)
   {
      if (Species_to_add[k])
      {
         write_file << listSpecies[k]->m_Name << "  ";

         if ((countSpecies+1)%10 == 0)
            write_file << endl << "    ";

         countSpecies += 1;
      }
   }

   write_file << "</speciesArray>" << endl;
   write_file << "    <reactionArray datasrc=\"#reaction_data\"/>" << endl;
   write_file << "    <state>" << endl;
   write_file << "      <temperature units=\"K\">300.0</temperature>" << endl;
   write_file << "      <pressure units=\"Pa\">101325.0</pressure>" << endl;
   write_file << "    </state>" << endl;
   write_file << "    <thermo model=\"IdealGas\"/>" << endl;
   write_file << "    <kinetics model=\"GasKinetics\"/>" << endl;
   write_file << "    <transport model=\"Mix\"/>" << endl;
   write_file << "  </phase>" << endl << endl;

   //Add the chemical species involved
   write_file << "  <!-- species definitions     -->" << endl;
   write_file << "  <speciesData id=\"species_data\">" << endl << endl;

   for (unsigned int k=0; k<listSpecies.size(); k++)
   {
      if (Species_to_add[k])
         write_file << "    " << listSpecies[k]->m_Description << endl << endl;
   }
   write_file << "  </speciesData>" << endl << endl;




   //Add the reactions
   write_file << "  <reactionData id=\"reaction_data\">" << endl << endl;
   write_file << "  </reactionData>" << endl;
   write_file << "</ctml>" << endl << endl;

   write_file.close();


}



















void Write::Write_reaction(ofstream& write, Reaction_ORCh* reac, vector<bool> Species_to_add, 
                           vector<Species_ORCh*> listSpecies, int j, vector<Reaction_ORCh*> listReactions, vector<bool> Reactions_fwd_to_add, 
                           vector<bool> Reactions_rev_to_add) const
{
   write << "    <!-- reaction " << j << " -->" << endl;
   write << "    <reaction";


   bool duplicate_already_written = false;

   if (reac->m_duplicate)
   {
      for (unsigned int j=0; j<listReactions.size(); j++)
      {
         if (duplicate_already_written == false)
         {
            if (Reactions_fwd_to_add[j] || Reactions_rev_to_add[j])
            {
               if (reac->m_equation == listReactions[j]->m_equation)
               {
                  //cout << reac->m_equation << endl;
                  //cout << listReactions[j]->m_equation << endl;
                  write << " duplicate=\"yes\"";
                  duplicate_already_written = true;
               }
            }
         }
      }
   }

   write << " reversible=\"";
   if (reac->m_reversible)
      write << "yes";
   else
      write << "no";

   if (dynamic_cast <ThreeBody *> (reac))
   {
      if (dynamic_cast <FalloffR *> (reac))
      {
         write << "\" type=\"falloff\"";
      }
      else
      {
         write << "\" type=\"threeBody\"";
      }
   }
   else
   {
      write << "\"";
   }


    
   if (fabs(reac->m_A) > 1)
   {
      reac->m_A = round_to_digits(double(reac->m_A), 3);
   }

   reac->m_b = round_to_digits(double(reac->m_b), 2);
   reac->m_E = round_to_digits(double(reac->m_E), 0);



   write << " id=\"" << j << "\">" << endl;
   write << "      <equation>" << reac->m_equation << "</equation>" << endl;
   write << "      <rateCoeff>" << endl;
   write << "        <Arrhenius>" << endl;
   write.precision(3);
   write << "           <A>" << reac->m_A << "</A>" << endl;
   write.precision(6);
   write << "           <b>" << reac->m_b << "</b>" << endl;
   write << "           <E" << reac->m_Etype << ">" << reac->m_E << "</E>" << endl;
   write << "        </Arrhenius>" << endl;

   if (dynamic_cast <FalloffR *> (reac))
   {
      write << "        <Arrhenius name=\"k0\">" << endl;
      write << "           <A>" << (dynamic_cast <FalloffR *> (reac))->m_A_low << "</A>" << endl;
      write << "           <b>" << (dynamic_cast <FalloffR *> (reac))->m_b_low << "</b>" << endl;
      write << "           <E" << (dynamic_cast <FalloffR *> (reac))->m_Etype_low << ">" <<
                                  (dynamic_cast <FalloffR *> (reac))->m_E_low << "</E>" << endl;
      write << "        </Arrhenius>" << endl;
   }


   if (dynamic_cast <ThreeBody *> (reac))
   {
      write << "        <efficiencies default=\"1.0\">";

      for (unsigned int k=0; k<(dynamic_cast <ThreeBody *> (reac))->m_TBconc.size(); k++)
      {
         for (unsigned int ka=0; ka<Species_to_add.size(); ka++)
         {
            if (((dynamic_cast <ThreeBody *> (reac))->m_NameTBconc[k] == listSpecies[ka]->m_Name)  && Species_to_add[ka])
            {
               write << (dynamic_cast <ThreeBody *> (reac))->m_NameTBconc[k] << ":";
               write << (dynamic_cast <ThreeBody *> (reac))->m_TBconc[k];
               if (k<(dynamic_cast <ThreeBody *> (reac))->m_TBconc.size()-1)
                  write << "  ";
            }
         }
      }
      write << "</efficiencies>" << endl;
   }

   if (dynamic_cast <Troe *> (reac))
   {
      write << "        <falloff type=\"Troe\">";
      for (unsigned int k=0; k<(dynamic_cast <Troe *> (reac))->m_TroeCoeffs.size(); k++)
      {
         write << (dynamic_cast <Troe *> (reac))->m_TroeCoeffs[k];
         if (k<(dynamic_cast <Troe *> (reac))->m_TroeCoeffs.size()-1)
            write << " ";
      }
      write << "</falloff>" << endl;
   }

   if (dynamic_cast <Lindemann *> (reac))
   {
      write << "        <falloff type=\"Lindemann\"/>" << endl;
   }

   write << "      </rateCoeff>" << endl;

   write << "      <reactants>";
   for (unsigned int k=0; k<reac->m_ReactantSpecies.size(); k++)
   {
      write << reac->m_ReactantSpecies[k] << ":";
      write << reac->m_ReactantStoichCoeffs[k];
      if (k<reac->m_ReactantSpecies.size()-1)
         write << " ";
   } 
   write << "</reactants>" << endl;

   write << "      <products>";
   for (unsigned int k=0; k<reac->m_ProductSpecies.size(); k++)
   {
      write << reac->m_ProductSpecies[k] << ":";
      write << reac->m_ProductStoichCoeffs[k];
      if (k<reac->m_ProductSpecies.size()-1)
         write << " ";
   } 
   write << "</products>" << endl;

   write << "    </reaction>" << endl << endl;
}







Write::~Write() //Destructeur
{}



