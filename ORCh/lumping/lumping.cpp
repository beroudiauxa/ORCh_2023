#include "lumping.h"
#include "../read_write/write.h"
#include "../optimisation/optimisation.h"

#include <sstream>


//---lumping---

lumping::lumping() //Constructeur
{}

//----------------------------------------Lumping----------------------------------------//

void lumping::lumping_mech(string mech, string trajectory_ref, string mech_desc) const
{
   vector<Species_ORCh*> listSpecies_ref;
   vector<Reaction_ORCh*> listReactions_ref;

   vector<Species_ORCh*> listSpecies_lumping;
   vector<Reaction_ORCh*> listReactions_lumping;

   Read *r = new Read();
   r->Read_species(mech, listSpecies_ref);
   r->Read_reactions(mech, listReactions_ref);

   r->Read_species(mech, listSpecies_lumping);
   r->Read_reactions(mech, listReactions_lumping);

   int nsp_ref = listSpecies_ref.size();
   int nreac_ref = listReactions_ref.size();


   //------------------------------------------------------------
   //1st find every group of isomers
   //------------------------------------------------------------
   
   vector<Lumping*> listLumpingGroups;
   lumping_groups(listLumpingGroups, listSpecies_ref);


   cout << "Into lumping function : ------------------------------------------------------------------------------------------------------------------------" << endl;

   cout << endl;







   //------------------------------------------------------------
   //Compute the concentration of each species from AI calculation
   //------------------------------------------------------------
/*

   ifstream in_ref("./outputs/Ignite1_luche.dat", ios::in);
   string line;
   int nbLines_ref = 0;
   while(getline(in_ref, line))
      nbLines_ref++;
   in_ref.close();

   vector<vector<double> > Yk_ref (nsp_ref, vector<double> (nbLines_ref, 0.0));
   vector<double> T_ref (nbLines_ref);
   vector<double> t_ref (nbLines_ref);
   vector<double> density_ref (nbLines_ref);
   vector<vector<double> > conc_ref (nsp_ref, vector<double> (nbLines_ref, 0.0));
   vector<double> molecularWeight (nsp_ref);

   for (int k=0; k<nsp_ref; k++)
   {
      molecularWeight[k] = listSpecies_ref[k]->m_C*12+listSpecies_ref[k]->m_H+listSpecies_ref[k]->m_C*16+listSpecies_ref[k]->m_C*14+listSpecies_ref[k]->m_Ar*40;
   }

   ifstream in_data("./outputs/Ignite1_luche.dat", ios::in);

   for (int i=0; i<nbLines_ref; i++)
   {
      in_data >> t_ref[i] >> T_ref[i];
      for (int k=0; k<nsp_ref; k++)
         in_data >> Yk_ref[k][i];
     in_data >> density_ref[i];

      for (int k=0; k<nsp_ref; k++)
         conc_ref[k][i] = density_ref[i]*Yk_ref[k][i]/molecularWeight[k];
   }
   in_data.close();
*/




   //------------------------------------------------------------
   //Compute the concentration of each species from MultipleInlet calculation
   //------------------------------------------------------------

   string path_ref = "./outputs/";
   stringstream s_traj;
   s_traj << trajectory_ref;
   path_ref.append(s_traj.str()).append("_0.dat");   


   cout << "lumping computed from the air trajectory : " << path_ref << endl << endl;

   ifstream in_ref(path_ref.c_str(), ios::in);
   string line;
   int nbLines_ref = 0;
   while(getline(in_ref, line))
      nbLines_ref++;
   in_ref.close();
   vector<vector<double> > Yk_ref (nsp_ref, vector<double> (nbLines_ref, 0.0));
   vector<double> T_ref (nbLines_ref);
   vector<double> t_ref (nbLines_ref);
   vector<double> density_ref (nbLines_ref);
   vector<vector<double> > conc_ref (nsp_ref, vector<double> (nbLines_ref, 0.0));
   vector<double> molecularWeight (nsp_ref);

   for (int k=0; k<nsp_ref; k++)
   {
      molecularWeight[k] = listSpecies_ref[k]->m_C*12+listSpecies_ref[k]->m_H+listSpecies_ref[k]->m_C*16+listSpecies_ref[k]->m_C*14+listSpecies_ref[k]->m_Ar*40;
   }

   ifstream in_data(path_ref.c_str(), ios::in);
  
   getline(in_data, line);//Skip the first commented line of the trajectory file

   for (int i=0; i<nbLines_ref; i++)
   {
       
      in_data >> t_ref[i];
      in_data >> T_ref[i];
      for (int k=0; k<nsp_ref; k++)
         in_data >> Yk_ref[k][i];
      //in_data >> density_ref[i];
      density_ref[i] = 1; //It is not printed within the Trajectory_ref0_luche.dat, this has to be improved 

   
      for (int k=0; k<nsp_ref; k++)
        conc_ref[k][i] = density_ref[i]*Yk_ref[k][i]/molecularWeight[k];
//         cout << conc_ref[k][i] << endl;
      
    }
    
   in_data.close();
   //for (int i=0; i<nbLines_ref; i++)
   //{
   //   cout << "Temperature at i " << T_ref[i] << endl;
   //}
   //getchar();




   //------------------------------------------------------------
   //Compute the alpha for each isomer
   //------------------------------------------------------------
   for (unsigned int group=0; group<listLumpingGroups.size(); group++)
   {
      double sum_conc = 0.0;
      vector<double> sum_conc_k (nsp_ref, 0.0);
      for (int i=0; i<nbLines_ref; i++)
      {
         for (unsigned int k=0; k<listLumpingGroups[group]->m_SpeciesList.size(); k++)
         {
            sum_conc += conc_ref[listLumpingGroups[group]->m_SpeciesNumber[k]][i];
            sum_conc_k[k] += conc_ref[listLumpingGroups[group]->m_SpeciesNumber[k]][i];
         }
      }
      for (unsigned int k=0; k<listLumpingGroups[group]->m_SpeciesList.size(); k++)
      {
         listLumpingGroups[group]->m_alphaSpecies[k] = sum_conc_k[k]/sum_conc;
      }
   }




   for (unsigned int group=0; group<listLumpingGroups.size(); group++)
   {
      cout << "Leading species " << listLumpingGroups[group]->m_Name << " with " << listLumpingGroups[group]->m_C << " carbone atoms " << listLumpingGroups[group]->m_H << " hydrogen and " << listLumpingGroups[group]->m_O << " oxygen" << endl;
      for (unsigned int k=0; k<listLumpingGroups[group]->m_SpeciesList.size(); k++)
      {
         cout << "      Associated to " << listLumpingGroups[group]->m_SpeciesList[k] << "  " << listLumpingGroups[group]->m_alphaSpecies[k] << endl;
      }
      cout << endl;
   }

 //  getchar();




   vector<bool> Species_to_add (nsp_ref, true);
   lumping_species(listLumpingGroups, listSpecies_ref, listSpecies_lumping, listReactions_lumping, Species_to_add);
   lumping_equation(listReactions_lumping);



   vector<bool> LumpedReaction (nreac_ref, false);
   for (int j=0; j<nreac_ref; j++)
   {
      for (unsigned int group=0; group<listLumpingGroups.size(); group++)
      {
         for (unsigned int kr=0; kr<listReactions_lumping[j]->m_ReactantSpecies.size(); kr++)
         {
            if (listLumpingGroups[group]->m_Name == listReactions_lumping[j]->m_ReactantSpecies[kr])
               LumpedReaction[j] = true;
         }
   
         for (unsigned int kp=0; kp<listReactions_lumping[j]->m_ProductSpecies.size(); kp++)
         {
            if (listLumpingGroups[group]->m_Name == listReactions_lumping[j]->m_ProductSpecies[kp])
               LumpedReaction[j] = true;
         }
      }
   }










   vector<bool> LinkedReaction (nreac_ref, false);
   //Remove Falloff reactions from the list of reaction to associate
   for (int j=0; j<nreac_ref; j++)
   {
      if (dynamic_cast <ThreeBody *> (listReactions_lumping[j]))
      {
         LinkedReaction[j] = true;
      }
      if (LumpedReaction[j] == false)
      {
         LinkedReaction[j] = true;
      }
   }

   for (int j=0; j<nreac_ref; j++)
   {
      vector<int> AssociatedReactions;
      if (LumpedReaction[j])
      {
         AssociatedReactions.push_back(j);

         vector<int> StoichCoeffsCompReactant (nsp_ref, 0);
         vector<int> StoichCoeffsCompProduct (nsp_ref, 0);
         for (int k=0; k<nsp_ref; k++)
         {
            for (unsigned int kr=0; kr<listReactions_lumping[j]->m_ReactantSpecies.size(); kr++)
            {
               if (listSpecies_lumping[k]->m_Name == listReactions_lumping[j]->m_ReactantSpecies[kr])
                  StoichCoeffsCompReactant[k] = listReactions_lumping[j]->m_ReactantStoichCoeffs[kr];
            }

            for (unsigned int kp=0; kp<listReactions_lumping[j]->m_ProductSpecies.size(); kp++)
            {
               if (listSpecies_lumping[k]->m_Name == listReactions_lumping[j]->m_ProductSpecies[kp])
                  StoichCoeffsCompProduct[k] = listReactions_lumping[j]->m_ProductStoichCoeffs[kp];
            }
         }
         int diff = 0;
         for (int k=0; k<nsp_ref; k++)
         {
            diff += fabs(StoichCoeffsCompReactant[k]-StoichCoeffsCompProduct[k]);
         }
         if (diff == 0)
         {
           listReactions_lumping[j]->m_equation = "";
         }
         else
         {

            for (int jb=j+1; jb<nreac_ref;  jb++)
            {
               if (LinkedReaction[jb] == false)
               {
                  //Add duplicate to similar reactions
                  vector<int> StoichCoeffsCompReactant1 (nsp_ref, 0);
                  vector<int> StoichCoeffsCompProduct1 (nsp_ref, 0);
                  for (int k=0; k<nsp_ref; k++)
                  {
                     for (unsigned int kr=0; kr<listReactions_lumping[jb]->m_ReactantSpecies.size(); kr++)
                     {
                        if (listSpecies_lumping[k]->m_Name == listReactions_lumping[jb]->m_ReactantSpecies[kr])
                           StoichCoeffsCompReactant1[k] = listReactions_lumping[jb]->m_ReactantStoichCoeffs[kr];
                     }

                     for (unsigned int kp=0; kp<listReactions_lumping[jb]->m_ProductSpecies.size(); kp++)
                     {
                        if (listSpecies_lumping[k]->m_Name == listReactions_lumping[jb]->m_ProductSpecies[kp])
                           StoichCoeffsCompProduct1[k] = listReactions_lumping[jb]->m_ProductStoichCoeffs[kp];
                     }
                  }

                  int diff1 = 0;
                  for (int k=0; k<nsp_ref; k++)
                  {
                     diff1 += fabs(StoichCoeffsCompReactant[k]-StoichCoeffsCompReactant1[k]);
                     diff1 += fabs(StoichCoeffsCompProduct[k]-StoichCoeffsCompProduct1[k]);
                  }
                  if (diff1 == 0)
                  {
                    listReactions_lumping[j]->m_duplicate = true;
                    listReactions_lumping[jb]->m_duplicate = true;
                    //listReactions_lumping[jb]->m_equation = "";
                    //listReactions_lumping[j]->m_A += listReactions_lumping[jb]->m_A;
                    
                    AssociatedReactions.push_back(jb);
                    LinkedReaction[jb] = true;
                  }
               }
            }
  

            //for (int j=0; j<AssociatedReactions.size(); j++)
            //{
            //   cout << listReactions_lumping[AssociatedReactions[j]]->m_equation << endl;
            //}
            //cout << endl;




         //AssociatedReactions sort
         if (AssociatedReactions.size() > 0)
         {
            vector<string> AssociatedSpecies (AssociatedReactions.size(), "");
            bool ReactantLumped = false;
            for (unsigned int jt=0; jt<AssociatedReactions.size(); jt++)
            {
               for (unsigned int group=0; group<listLumpingGroups.size(); group++)
               {
                  for (unsigned int kr=0; kr<listReactions_ref[AssociatedReactions[jt]]->m_ReactantSpecies.size(); kr++)
                  {
                     if (listLumpingGroups[group]->m_Name == listReactions_lumping[AssociatedReactions[jt]]->m_ReactantSpecies[kr])
                     {
                        for (unsigned int kg=0; kg<listLumpingGroups[group]->m_SpeciesList.size(); kg++)
                        {
                           if (listLumpingGroups[group]->m_SpeciesList[kg] == listReactions_ref[AssociatedReactions[jt]]->m_ReactantSpecies[kr])
                           {
                              AssociatedSpecies[jt] = listLumpingGroups[group]->m_SpeciesList[kg];
                              ReactantLumped = true;
                           }
                        }
                     }
                  }

                  for (unsigned int kp=0; kp<listReactions_ref[AssociatedReactions[jt]]->m_ProductSpecies.size(); kp++)
                  {
                  }
               }

            }




            vector<double> ReactionRate (nbLines_ref, 0.0);
            double max_ReactionRate = 0.0;
            vector<double> conc_lumped (nbLines_ref, 0.0);
            int coeff = 1;


            vector<vector<double> > min_max_A (1, vector<double> (2, 0));
            vector<vector<double> > min_max_b (1, vector<double> (2, 0));
            vector<vector<double> > min_max_E (1, vector<double> (2, 0));


            if (ReactantLumped)
            {
               for (unsigned int jt=0; jt<AssociatedReactions.size(); jt++)
               {
                  getReactionRate(listLumpingGroups, listReactions_ref, listSpecies_ref, ReactionRate, max_ReactionRate, conc_ref, T_ref, conc_lumped, AssociatedReactions, AssociatedSpecies, coeff, true, jt);

                  //Optimise a single reaction but with lumped species (concentration of the lumped species is different)
                  cout << "Optimise a single reaction with lumped species" << endl;
                  double A = listReactions_ref[AssociatedReactions[jt]]->m_A;
                  double lumpedA = listReactions_lumping[AssociatedReactions[jt]]->m_A;
                  double E = listReactions_ref[AssociatedReactions[jt]]->m_E;
                  double b = listReactions_ref[AssociatedReactions[jt]]->m_b;

                  cout << "A " << A << " lumpedA " << lumpedA << endl;
                  cout << "E " << E << endl;
                  cout << "b " << b << endl;

                     min_max_A[0][0] = log10(A-0.2*fabs(A));
                     min_max_A[0][1] = log10(A+fabs(A));
                
                  if (min_max_A[0][0] > log10(lumpedA-0.2*fabs(lumpedA)))
                     min_max_A[0][0] = log10(lumpedA-0.2*fabs(lumpedA));

                  if (min_max_A[0][1] < log10(lumpedA+fabs(lumpedA)))
                     min_max_A[0][1] = log10(lumpedA+fabs(lumpedA));

                  min_max_b[0][0] = (b-fabs(0.2*b));
                  min_max_b[0][1] = (b+fabs(0.2*b));

                  min_max_E[0][0] = (E-fabs(0.2*E))/8.314;
                  min_max_E[0][1] = (E+fabs(0.2*E))/8.314;

                  //cout << "min_A " << min_max_A[0][0] << endl;
                  //cout << "max_A " << min_max_A[0][1] << endl;
                  //cout << endl;
                  //cout << "min_b " << min_max_b[0][0] << endl;
                  //cout << "max_b " << min_max_b[0][1] << endl;
                  //cout << endl;
                  //cout << "min_E " << min_max_E[0][0] << endl;
                  //cout << "max_E " << min_max_E[0][1] << endl;
                  
                 lumping_optimise(min_max_A, min_max_b, min_max_E, listReactions_lumping, AssociatedReactions, coeff, t_ref, T_ref, conc_lumped, ReactionRate, max_ReactionRate, true, jt);
               }
            }
            else
            { 
               /*
               //DO NOTHING IN THIS LOOP
               cout << "Do Nothing, you are with a single reaction with a lumped species only on the products" << endl;
               */
            }


            if (AssociatedReactions.size() > 1)
            {
               if (ReactantLumped)
               {
                  getReactionRate(listLumpingGroups, listReactions_ref, listSpecies_ref, ReactionRate, max_ReactionRate, conc_ref, T_ref, conc_lumped, AssociatedReactions, AssociatedSpecies, coeff, false, 0);

                  double max_A = -1e25;
                  double max_b = -1e25;
                  double max_E = -1e25;
                  double min_A = 1e25;
                  double min_b = 1e25;
                  double min_E = 1e25;

                  for (unsigned int jt=0; jt<AssociatedReactions.size(); jt++)
                  {
                     //Compute Arrhenius for each reaction
                     double A = listReactions_ref[AssociatedReactions[jt]]->m_A;
                     double lumpedA = listReactions_lumping[AssociatedReactions[jt]]->m_A;
                     double E = listReactions_ref[AssociatedReactions[jt]]->m_E;
                     double b = listReactions_ref[AssociatedReactions[jt]]->m_b;

                     cout << listReactions_ref[AssociatedReactions[jt]]->m_equation << endl;
                     cout << listReactions_lumping[AssociatedReactions[jt]]->m_equation << endl;

                     if (min_A > A) min_A = A;
                     if (min_A > lumpedA) min_A = lumpedA;
                     if (min_E > E) min_E = E;
                     if (min_b > b) min_b = b;

                     if (max_A < A) max_A = A;
                     if (max_A < lumpedA) max_A = lumpedA;
                     if (max_E < E) max_E = E;
                     if (max_b < b) max_b = b;
                  }

                  cout << "Multiple reactions to associate with lumped Reactants " << endl;
                  cout << endl << endl;

                  min_max_A[0][0] = log10(AssociatedReactions.size()*min_A-0.2*fabs(AssociatedReactions.size()*min_A));
                  min_max_A[0][1] = log10(AssociatedReactions.size()*max_A+3*fabs(AssociatedReactions.size()*max_A));

                  min_max_b[0][0] = (min_b-fabs(0.8*min_b));
                  min_max_b[0][1] = (max_b+fabs(2*max_b));

                  min_max_E[0][0] = (min_E-fabs(0.8*min_E))/8.314;
                  min_max_E[0][1] = (max_E+fabs(2*max_E))/8.314;

                  cout << "min_A " << min_max_A[0][0] << endl;
                  cout << "max_A " << min_max_A[0][1] << endl;
                  cout << endl;
                  cout << "min_b " << min_max_b[0][0] << endl;
                  cout << "max_b " << min_max_b[0][1] << endl;
                  cout << endl;
                  cout << "min_E " << min_max_E[0][0] << endl;
                  cout << "max_E " << min_max_E[0][1] << endl;
                  
                  lumping_optimise(min_max_A, min_max_b, min_max_E, listReactions_lumping, AssociatedReactions, coeff, t_ref, T_ref, conc_lumped, ReactionRate, max_ReactionRate, false, 0);
                  for (unsigned int jt=0; jt<AssociatedReactions.size(); jt++)
                  {
                     cout << "lumping after  " << listReactions_lumping[AssociatedReactions[jt]]->m_equation << endl;
                     cout << "lumping duplicate  " << listReactions_lumping[AssociatedReactions[jt]]->m_duplicate << endl;
                     cout << endl << endl;
                  }
               }
               else
               {
                  cout << "Multiple reactions with lumped Products, no optimization" << endl;
                  vector<bool> LumpedR (AssociatedReactions.size(), false);
                  for (unsigned int jt=0; jt<AssociatedReactions.size(); jt++)
                  {
                     //double A = listReactions_lumping[AssociatedReactions[jt]]->m_A;
                     double E = listReactions_lumping[AssociatedReactions[jt]]->m_E;
                     double b = listReactions_lumping[AssociatedReactions[jt]]->m_b;

                     for (unsigned int jb=jt+1; jb<AssociatedReactions.size(); jb++)
                     {
                        if (LumpedR[jb] == false)
                        {
                           //double A1 = listReactions_lumping[AssociatedReactions[jb]]->m_A;
                           double E1 = listReactions_lumping[AssociatedReactions[jb]]->m_E;
                           double b1  = listReactions_lumping[AssociatedReactions[jb]]->m_b;
                           if (E1 == E && b1 == b)
                           {
                              LumpedR[jb] = true;
                              LumpedR[jt] = true;
                              listReactions_lumping[AssociatedReactions[jb]]->m_equation = "";
                              listReactions_lumping[AssociatedReactions[jt]]->m_A += listReactions_lumping[AssociatedReactions[jb]]->m_A;
                           } //end if (E1 == E && b1 == b)
                        } //end if (LumpedR[jb] == false)
                     } //end for (int jb=jt+1 jb<AssociatedReactions.size(); jb++)
                  } //end (int jt=0; jt<AssociatedReactions.size(); jt++)
                  cout << endl << endl;
                  cout << "Go on to optimisation" << endl;

                  //------------------------------------------------------//
                  //Optimise between two reactions with same lumped product
                  //------------------------------------------------------//
                  







               } //else (for ProductLumped)
            }





            }
         }
      }
   }

   Write *w = new Write();
   w->Write_xml_file(mech, mech_desc, "./outputs/mechanisms/Lumped.xml", Species_to_add, vector<bool> (nreac_ref, true), false, vector<double> (), vector<double> (), vector<double> (), listSpecies_lumping, listReactions_lumping);






}








































void lumping::getReactionRate(vector<Lumping*> listLumpingGroups, vector<Reaction_ORCh*> listReactions_ref, vector<Species_ORCh*> listSpecies_ref, vector<double> &ReactionRate, double &max_ReactionRate, vector<vector<double> > conc_ref, vector<double> T_ref, vector<double> &conc_lumped, vector<int> AssociatedReactions, vector<string> AssociatedSpecies, int &coeff, bool singleReaction, int jt_ref) const
{

   int nsp_ref = listSpecies_ref.size();
   for (unsigned int i=0; i<T_ref.size(); i++)
   {
      conc_lumped[i] = 0.0;
      ReactionRate[i] = 0.0;
   }
   max_ReactionRate = 0.0;

   int group_ref = -1;
   for (unsigned int group=0; group<listLumpingGroups.size(); group++)
   {
      for (unsigned int k=0; k<listLumpingGroups[group]->m_SpeciesList.size(); k++)
      {
         for (unsigned int jt_b=0; jt_b<AssociatedReactions.size(); jt_b++)
         {
            if (listLumpingGroups[group]->m_SpeciesList[k] == AssociatedSpecies[jt_b])
            {
               group_ref = group;
            }
         }
      }
   }

   if (group_ref != -1)
   {
      for (unsigned int i=0; i<T_ref.size(); i++)
      {
         conc_lumped[i] = 0.0;
         for (int k=0; k<nsp_ref; k++)
         {
            for (unsigned int kb=0; kb<listLumpingGroups[group_ref]->m_SpeciesList.size(); kb++)
            {
               if (listSpecies_ref[k]->m_Name == listLumpingGroups[group_ref]->m_SpeciesList[kb])
               {

                  //cout << "node " << i << " conc_ref[k][i] " << conc_ref[k][i] << endl;

                  conc_lumped[i] += conc_ref[k][i];
               }
            }
         }
      }
   }


   if (singleReaction)
   {
      double A = listReactions_ref[AssociatedReactions[jt_ref]]->m_A;
      double E = listReactions_ref[AssociatedReactions[jt_ref]]->m_E;
      double b = listReactions_ref[AssociatedReactions[jt_ref]]->m_b;

 
      for (unsigned int i=0; i<T_ref.size(); i++)
      {
         double conc_term = 0.0;
         for (int k=0; k<nsp_ref; k++)
         {
            if (listSpecies_ref[k]->m_Name == AssociatedSpecies[jt_ref])
            {
               for (unsigned int kr=0; kr<listReactions_ref[AssociatedReactions[jt_ref]]->m_ReactantSpecies.size(); kr++)
               {
                  coeff = listReactions_ref[AssociatedReactions[jt_ref]]->m_ReactantStoichCoeffs[kr];
                  conc_term = pow(conc_ref[k][i], coeff);
               }
            }
         }

         //cout << "Node " << i << " species conc lumped " << conc_lumped[i] << endl;
         //getchar();
 
         if (conc_lumped[i] != 0.0 || i == 0)
            ReactionRate[i] = A*pow(T_ref[i],b)*exp(-E/(8.314*T_ref[i]))*conc_term;
         else //Concentrations will be the same, just adapt the Arrhenius terms
            ReactionRate[i] = A*pow(T_ref[i],b)*exp(-E/(8.314*T_ref[i]));
 
         if (ReactionRate[i] > max_ReactionRate)
            max_ReactionRate = ReactionRate[i];
      }
   }
   else
   {
      for (unsigned int jt=0; jt<AssociatedReactions.size(); jt++)
      {
         //Compute Arrhenius for each reaction
         double A = listReactions_ref[AssociatedReactions[jt]]->m_A;
         double E = listReactions_ref[AssociatedReactions[jt]]->m_E;
         double b = listReactions_ref[AssociatedReactions[jt]]->m_b;

         int nsp_ref = listSpecies_ref.size();

         for (unsigned int i=0; i<T_ref.size(); i++)
         {
            double conc_term = 0.0;
            for (int k=0; k<nsp_ref; k++)
            {
               if (listSpecies_ref[k]->m_Name == AssociatedSpecies[jt])
               {
                  for (unsigned int kr=0; kr<listReactions_ref[AssociatedReactions[jt]]->m_ReactantSpecies.size(); kr++)
                  {
                     coeff = listReactions_ref[AssociatedReactions[jt]]->m_ReactantStoichCoeffs[kr];
                     conc_term = pow(conc_ref[k][i], coeff);
                  }
               }
            }

            if (conc_lumped[i] != 0.0 || i == 0)
               ReactionRate[i] += A*pow(T_ref[i],b)*exp(-E/(8.314*T_ref[i]))*conc_term;
            else //Concentrations will be the same, just adapt the Arrhenius terms
               ReactionRate[i] += A*pow(T_ref[i],b)*exp(-E/(8.314*T_ref[i]));

            if (ReactionRate[i] > max_ReactionRate)
               max_ReactionRate = ReactionRate[i];
         }
      }
   }


   //ofstream storeConc ("outputs/conc.dat");
   //for (int i=0; i<T_ref.size(); i++)
   //{
   //   storeConc << conc_lumped[i] << "  ";
   //   for (int k=0; k<nsp_ref; k++)
   //   {
   //      for (int kb=0; kb<listLumpingGroups[group_ref]->m_SpeciesList.size(); kb++)
   //      {
   //         if (listSpecies_ref[k]->m_Name == listLumpingGroups[group_ref]->m_SpeciesList[kb])
   //         {
   //            storeConc << conc_ref[k][i] << "  ";
   //         }
   //      }
   //   }
   //   storeConc << endl;

   //}
   //storeConc.close();





}


void lumping::lumping_optimise(vector<vector<double> > min_max_A, vector<vector<double> > min_max_b, vector<vector<double> > min_max_E,
                               vector<Reaction_ORCh*> &listReactions_lumping, vector<int> AssociatedReactions, int coeff, 
                               vector<double> t_ref, vector<double> T_ref, vector<double> conc_lumped, vector<double> ReactionRate, 
                               double max_ReactionRate, bool singleReaction, int jt_ref) const
{

   //Optimise the constants of the lumped reaction
   Optim *o = new Optim();
               
   int PopSize =  250;
   int MaxAllowableGenerations = 300;

   double CrossoverRate = 0.75;
   double MutationRate = 0.02;

   vector<int> nBits_A (1, 0);
   double accuracyA = 0.000001;

   vector<int> nBits_b (1, 0);
   double accuracyb = 0.0001;

   vector<int> nBits_E (1, 0);
   double accuracyE = 0.0001;

   if (min_max_A[0][0] == min_max_A[0][1])
      nBits_A[0] = 0;
   else
      nBits_A[0] = int(log((min_max_A[0][1]-min_max_A[0][0])/accuracyA)/log(2)) +1;

   if (min_max_b[0][0] == min_max_b[0][1])
      nBits_b[0] = 0;
   else
      nBits_b[0] = int(log((min_max_b[0][1]-min_max_b[0][0])/accuracyb)/log(2)) +1;

   if (min_max_E[0][0] == min_max_E[0][1])
      nBits_E[0] = 0;
   else
      nBits_E[0] = int(log((min_max_E[0][1]-min_max_E[0][0])/accuracyE)/log(2)) +1;

   cout << "nBits_A " << nBits_A[0] << endl;
   cout << "nBits_b " << nBits_b[0] << endl;
   cout << "nBits_E " << nBits_E[0] << endl;

   //seed the random number generator
   srand((int)time(NULL));

   vector<Chromosome*> Population;
   for (int k=0; k<PopSize; k++)
   {
      Population.push_back(new Chromosome("", 0.0));
      Population[k]->m_bits = "a";
   }

   int NbTotBits = nBits_A[0] + nBits_b[0] + nBits_E[0];

   for (int k=0; k<PopSize; k++)
   {
      string bits_random = o->GetRandomBits(NbTotBits);
      Population[k]->m_bits = bits_random;
   }

   double A_solution = 0.0;
   double E_solution = 0.0;
   double b_solution = 0.0;

   //we will set this flag if a solution has been found
   double best_fitness = 0.0;
   bool bFound = false;
   int count = 0;
   
   while(!bFound)
   {
      vector<double> A_lumped (1, 0.0);
      vector<double> b_lumped (1, 0.0);
      vector<double> E_lumped (1, 0.0);
      for (int i=0; i<PopSize; i++)
      {
         o->Translate_Bits(A_lumped, b_lumped, E_lumped, nBits_A, nBits_b, nBits_E, Population[i]->m_bits, min_max_A, min_max_b, min_max_E);

         //cout << "A_lumped " << A_lumped[0] << endl;
         //cout << "b_lumped " << b_lumped[0] << endl;
         //cout << "E_lumped " << E_lumped[0] << endl;
         //getchar();
         
         //cout << "T_ref.size() " << T_ref.size() << endl;
         
         double fitness_lumping = 0.0;
         vector<double> ReactionRate_Lumped (T_ref.size(), 0.0);
         for (unsigned int n=0 /*To avoid accounting for the first point.*/; n<T_ref.size(); n++)
         {
            if (conc_lumped[n] != 0.0 || n == 0)
               ReactionRate_Lumped[n] = A_lumped[0]*pow(T_ref[n],b_lumped[0])*exp(-E_lumped[0]/(8.314*T_ref[n]))*pow(conc_lumped[n], coeff);
            else 
               ReactionRate_Lumped[n] = A_lumped[0]*pow(T_ref[n],b_lumped[0])*exp(-E_lumped[0]/(8.314*T_ref[n]));

            if (ReactionRate[n] > 1e-4*max_ReactionRate)
            {
               fitness_lumping -= fabs(ReactionRate_Lumped[n]-ReactionRate[n])/ReactionRate[n];
            }
         }
         Population[i]->m_fitness = fitness_lumping;
      } //end for on PopSize

      //for (int i=0; i<PopSize; i++)
      //   cout << " Population[" << i << "]->m_fitness " << Population[i]->m_fitness << endl;
      //getchar();

      vector<Chromosome*> Temporary;
      vector<Chromosome*> Sort;
       
      for (int k=0; k<PopSize; k++)
      {
         Temporary.push_back(new Chromosome("", 0.0));
         Sort.push_back(new Chromosome("", 0.0));
      }
       
      o->Rank(Sort, Population, PopSize);
       
      double Total_fitness_linear_ranking = 0.0;

      best_fitness = Sort[PopSize-1]->m_fitness;

      o->Linear_Ranking(Sort, PopSize);
      
      for (int i=0; i<PopSize; i++)
         Total_fitness_linear_ranking += Sort[i]->m_fitness;
      
      int cPop = 0;
      o->Elitism(cPop, Temporary, Sort, PopSize);
      cPop = NB_ELITISM;

      //loop until we have created PopSize new chromosomes
      while (cPop < PopSize)
      {
         // we are going to create the new population by grabbing members of the old population
         // two at a time via roulette wheel selection.
         string offspring1 = o->Roulette(int(Total_fitness_linear_ranking), Sort, PopSize);
         string offspring2 = o->Roulette(int(Total_fitness_linear_ranking), Sort, PopSize);

         //add crossover dependent on the crossover rate
         o->Crossover(offspring1, offspring2, CrossoverRate);

         //now mutate dependent on the mutation rate
         o->Mutate(offspring1, MutationRate);
         o->Mutate(offspring2, MutationRate);

        //add these offspring to the new population. (assigning zero as their
        //fitness scores)
        Temporary[cPop]->m_bits = offspring1;
        Temporary[cPop]->m_fitness = 0.0;
        cPop += 1;

        if (cPop < PopSize)
        {
           Temporary[cPop]->m_bits = offspring2;
           Temporary[cPop]->m_fitness = 0.0;
           cPop += 1;
        }
         
      }//end loop

      //copy temp population into main population array
      for (int i=0; i<PopSize; i++)
         Population[i] = Temporary[i];

      //Store the best solution
      o->Translate_Bits(A_lumped, b_lumped, E_lumped, nBits_A, nBits_b, nBits_E, Population[0]->m_bits, min_max_A, min_max_b, min_max_E);
      A_solution = A_lumped[0];
      b_solution = b_lumped[0];
      E_solution = E_lumped[0];

      ++count;
      if ( count > MaxAllowableGenerations)
         bFound = true;
   } //end while bfound


   cout << "fitness " << fabs(best_fitness) << endl;

   if (singleReaction)
   {
      listReactions_lumping[AssociatedReactions[jt_ref]]->m_A = A_solution;
      listReactions_lumping[AssociatedReactions[jt_ref]]->m_b = b_solution;
      listReactions_lumping[AssociatedReactions[jt_ref]]->m_E = E_solution;
      cout << "A_solution " << A_solution << endl;
      cout << "b_solution " << b_solution << endl;
      cout << "E_solution " << E_solution << endl;
      cout << listReactions_lumping[AssociatedReactions[jt_ref]]->m_equation << endl;
   }

   if (singleReaction == false)
   {
      if (fabs(best_fitness) < 2) 
      {
         listReactions_lumping[AssociatedReactions[0]]->m_A = A_solution;
         listReactions_lumping[AssociatedReactions[0]]->m_b = b_solution;
         listReactions_lumping[AssociatedReactions[0]]->m_E = E_solution;

         for (unsigned int j=1; j<AssociatedReactions.size(); j++)
            listReactions_lumping[AssociatedReactions[j]]->m_equation = "";
      }
      else
      {
         for (unsigned int j=0; j<AssociatedReactions.size(); j++)
            listReactions_lumping[AssociatedReactions[j]]->m_duplicate = true;
      }
   }


   ofstream storeRates ("outputs/compareRates.dat");
   for (unsigned int i=0; i<T_ref.size(); i++)
   {
       double ReactionRate_Lumped;
   
       if (conc_lumped[i] != 0.0)
          ReactionRate_Lumped = A_solution*pow(T_ref[i],b_solution)*exp(-E_solution/(8.314*T_ref[i]))*pow(conc_lumped[i], coeff);
       else 
          ReactionRate_Lumped = A_solution*pow(T_ref[i],b_solution)*exp(-E_solution/(8.314*T_ref[i]));
   
      storeRates << t_ref[i] << "  " << ReactionRate[i] << "  " << ReactionRate_Lumped << "  " << endl;
   }
   storeRates.close();
   
//   getchar();







}







void lumping::lumping_equation(vector<Reaction_ORCh*> &listReactions_lumping) const
{

   int nreac_ref = listReactions_lumping.size();

   for (int j=0; j<nreac_ref; j++)
   {
      string equation = "";
      for (unsigned int kr=0; kr<listReactions_lumping[j]->m_ReactantSpecies.size(); kr++)
      {
         if (listReactions_lumping[j]->m_ReactantStoichCoeffs[kr] > 1)
         {
            stringstream s_coeff;
            s_coeff << listReactions_lumping[j]->m_ReactantStoichCoeffs[kr];
            equation.append(s_coeff.str()).append(" ");
         }

         equation.append(listReactions_lumping[j]->m_ReactantSpecies[kr]);
            
         if (kr < listReactions_lumping[j]->m_ReactantSpecies.size()-1)
            equation.append(" + ");
         else
         {
            if (dynamic_cast <ThreeBody *> (listReactions_lumping[j]))
            {
               if (dynamic_cast <FalloffR *> (listReactions_lumping[j]))
               {
                  equation.append (" (+ M)");
               }
               else
               {
                  equation.append (" + M");
               }
            }

            if (listReactions_lumping[j]->m_reversible)
               equation.append(" [=] ");
            else
               equation.append(" =] ");
         }
      }
      for (unsigned int kp=0; kp<listReactions_lumping[j]->m_ProductSpecies.size(); kp++)
      {
         if (listReactions_lumping[j]->m_ProductStoichCoeffs[kp] > 1)
         {
            stringstream s_coeff;
            s_coeff << listReactions_lumping[j]->m_ProductStoichCoeffs[kp];
            equation.append(s_coeff.str()).append(" ");
         }

         equation.append(listReactions_lumping[j]->m_ProductSpecies[kp]);

         if (kp < listReactions_lumping[j]->m_ProductSpecies.size()-1)
            equation.append(" + ");
         else
         {
            if (dynamic_cast <ThreeBody *> (listReactions_lumping[j]))
            {
               if (dynamic_cast <FalloffR *> (listReactions_lumping[j]))
               {
                  equation.append (" (+ M)");
               }
               else
               {
                  equation.append (" + M");
               }
            }
         }
      }
      listReactions_lumping[j]->m_equation = equation;
   }


}





void lumping::lumping_groups(vector<Lumping*> &listLumpingGroups, vector<Species_ORCh*> listSpecies_ref) const
{

   int nsp_ref = listSpecies_ref.size();
   vector<bool> Selected_for_lumping (nsp_ref, false);

   for (int k=0; k<nsp_ref; k++)
   {
      if (Selected_for_lumping[k] == false)
      {
         bool new_group = true;
         for (int kl=k+1; kl<nsp_ref; kl++)
         {
            if (Selected_for_lumping[kl] == false)
            {
               if ((listSpecies_ref[k]->m_C == listSpecies_ref[kl]->m_C) &&
                   (listSpecies_ref[k]->m_H == listSpecies_ref[kl]->m_H) &&
                   (listSpecies_ref[k]->m_O == listSpecies_ref[kl]->m_O) &&
                   (listSpecies_ref[k]->m_N == listSpecies_ref[kl]->m_N))
               {
                  Selected_for_lumping[k] = true;
                  Selected_for_lumping[kl] = true;

                  if (new_group)
                  {
                     new_group = false;
                     string name = "";
                     if (listSpecies_ref[k]->m_C != 0)
                     {
                        name.append("C");
                        if (listSpecies_ref[k]->m_C > 1)
                        {
                           stringstream s_C;
                           s_C << listSpecies_ref[k]->m_C;
                           name.append(s_C.str());
                        }
                     }
                     if (listSpecies_ref[k]->m_H != 0)
                     {
                        name.append("H");
                        if (listSpecies_ref[k]->m_H > 1)
                        {
                           stringstream s_H;
                           s_H << listSpecies_ref[k]->m_H;
                           name.append(s_H.str());
                        }
                     }
                     if (listSpecies_ref[k]->m_O != 0)
                     {
                        name.append("O");
                        if (listSpecies_ref[k]->m_O > 1)
                        {
                           stringstream s_O;
                           s_O << listSpecies_ref[k]->m_O;
                           name.append(s_O.str());
                        }
                     }
                     if (listSpecies_ref[k]->m_N != 0)
                     {
                        name.append("N");
                        if (listSpecies_ref[k]->m_N > 1)
                        {
                           stringstream s_N;
                           s_N << listSpecies_ref[k]->m_N;
                           name.append(s_N.str());
                        }
                     }
                     name.append("(L)");
                     listLumpingGroups.push_back(new Lumping (name, listSpecies_ref[k]->m_C, listSpecies_ref[k]->m_H, listSpecies_ref[k]->m_O, listSpecies_ref[k]->m_N, 
                                                 vector<string> (1, listSpecies_ref[k]->m_Name), vector<double> (1, 0.0), vector<int> (1, k), vector<vector<int> > (), vector<vector<string> > ()));
                  }

                  int number = listLumpingGroups.size();
                  listLumpingGroups[number-1]->m_SpeciesList.push_back(listSpecies_ref[kl]->m_Name);
                  listLumpingGroups[number-1]->m_alphaSpecies.push_back(0.0);
                  listLumpingGroups[number-1]->m_SpeciesNumber.push_back(kl);
               }
            }
         }
      }
   }



}






void lumping::lumping_species(vector<Lumping*> listLumpingGroups, vector<Species_ORCh*> listSpecies_ref, vector<Species_ORCh*> &listSpecies_lumping, vector<Reaction_ORCh*> &listReactions_lumping, vector<bool> &Species_to_add) const
{


   int nreac_ref = listReactions_lumping.size();

   for (unsigned int group=0; group<listLumpingGroups.size(); group++)
   {
      if (group == 1)
      {
            for (int j=0; j<nreac_ref; j++)
            {
               for (unsigned int k=0; k<listLumpingGroups[group]->m_SpeciesList.size(); k++)
               {
                  for (unsigned int jk=0; jk<listReactions_lumping[j]->m_ReactantSpecies.size(); jk++)
                  {
                     if (listReactions_lumping[j]->m_ReactantSpecies[jk] == listLumpingGroups[group]->m_SpeciesList[k])
                     {
                        listReactions_lumping[j]->m_ReactantSpecies[jk] = listLumpingGroups[group]->m_Name;
                        listReactions_lumping[j]->m_A *= pow(listLumpingGroups[group]->m_alphaSpecies[k], listReactions_lumping[j]->m_ReactantStoichCoeffs[jk]);

                        if (dynamic_cast <FalloffR *> (listReactions_lumping[j]))
                           (dynamic_cast <FalloffR *> (listReactions_lumping[j]))->m_A_low *= pow(listLumpingGroups[group]->m_alphaSpecies[k], listReactions_lumping[j]->m_ReactantStoichCoeffs[jk]);
                     }
                  }
   
                  for (unsigned int jk=0; jk<listReactions_lumping[j]->m_ProductSpecies.size(); jk++)
                  {
                     if (listReactions_lumping[j]->m_ProductSpecies[jk] == listLumpingGroups[group]->m_SpeciesList[k])
                        listReactions_lumping[j]->m_ProductSpecies[jk] = listLumpingGroups[group]->m_Name;
                  }
               }
            }


            for (unsigned int k=1; k<listLumpingGroups[group]->m_SpeciesList.size(); k++)
            {
               Species_to_add[listLumpingGroups[group]->m_SpeciesNumber[k]] = false;
            }
            listSpecies_lumping[listLumpingGroups[group]->m_SpeciesNumber[0]]->m_Name = listLumpingGroups[group]->m_Name;

            for (int p=0; p<7; p++)
            {
               double sum_NASA_high = 0.0;
               double sum_NASA_low = 0.0;
               for (unsigned int k=0; k<listLumpingGroups[group]->m_SpeciesList.size(); k++)
               {
                  sum_NASA_high += listSpecies_lumping[listLumpingGroups[group]->m_SpeciesNumber[k]]->m_NASACoeffs_highT[p]*listLumpingGroups[group]->m_alphaSpecies[k];
                  sum_NASA_low += listSpecies_lumping[listLumpingGroups[group]->m_SpeciesNumber[k]]->m_NASACoeffs_lowT[p]*listLumpingGroups[group]->m_alphaSpecies[k];
               }
               listSpecies_lumping[listLumpingGroups[group]->m_SpeciesNumber[0]]->m_NASACoeffs_highT[p] = sum_NASA_high;
               listSpecies_lumping[listLumpingGroups[group]->m_SpeciesNumber[0]]->m_NASACoeffs_lowT[p] = sum_NASA_low;
            }
            double sum_DipoleMoment = 0.0;
            double sum_polarizability = 0.0;
            double sum_rotRelax = 0.0;
            double sum_LJ_welldepth = 0.0;
            double sum_LJ_diameter = 0.0;
            for (unsigned int k=0; k<listLumpingGroups[group]->m_SpeciesList.size(); k++)
            {
               sum_DipoleMoment += listSpecies_ref[listLumpingGroups[group]->m_SpeciesNumber[k]]->m_dipoleMoment*listLumpingGroups[group]->m_alphaSpecies[k];
               sum_polarizability += listSpecies_ref[listLumpingGroups[group]->m_SpeciesNumber[k]]->m_polarizability*listLumpingGroups[group]->m_alphaSpecies[k];
               sum_rotRelax += listSpecies_ref[listLumpingGroups[group]->m_SpeciesNumber[k]]->m_rot_relax*listLumpingGroups[group]->m_alphaSpecies[k];
               sum_LJ_diameter += listSpecies_ref[listLumpingGroups[group]->m_SpeciesNumber[k]]->m_LJ_diameter*listLumpingGroups[group]->m_alphaSpecies[k];
               sum_LJ_welldepth += listSpecies_ref[listLumpingGroups[group]->m_SpeciesNumber[k]]->m_LJ_welldepth*listLumpingGroups[group]->m_alphaSpecies[k];
            }
            listSpecies_lumping[listLumpingGroups[group]->m_SpeciesNumber[0]]->m_dipoleMoment = sum_DipoleMoment;
            listSpecies_lumping[listLumpingGroups[group]->m_SpeciesNumber[0]]->m_polarizability = sum_polarizability;
            listSpecies_lumping[listLumpingGroups[group]->m_SpeciesNumber[0]]->m_rot_relax = sum_rotRelax;
            listSpecies_lumping[listLumpingGroups[group]->m_SpeciesNumber[0]]->m_LJ_diameter = sum_LJ_diameter;
            listSpecies_lumping[listLumpingGroups[group]->m_SpeciesNumber[0]]->m_LJ_welldepth = sum_LJ_welldepth;
      }
   } //group 
}










lumping::~lumping() //Destructeur
{}

