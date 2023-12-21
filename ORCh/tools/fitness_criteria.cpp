#include <fstream>
#include <sstream>
#include "fitness_criteria.h"





void fit_function_0D(string mech_ref, string mech, double& fitness, int nbInlets, vector<string> listTargets, vector<string> trajectory_ref, int nbToKeep, string step, int rank)
{

   vector<Species_ORCh*> listSpecies_ref;
   vector<Reaction_ORCh*> listReactions_ref;

   Read *r = new Read();
   r->Read_species(mech_ref, listSpecies_ref);
   r->Read_reactions(mech_ref, listReactions_ref);

   int nsp_init = listSpecies_ref.size();
   int nsp;

   fitness = 0.0;

   for (int n=0; n<nbInlets; n++)
   {
   
      vector<Species_ORCh*> listSpecies;
      
      stringstream s_step;
      s_step << step;

      stringstream s_nInlet;
      s_nInlet << n;

      if (trajectory_ref[0] == "" && step != "ComputeTrajectories")
      {
        trajectory_ref[0] = "Ref_";
        trajectory_ref[0].append(s_step.str());
        stringstream s_nsp_init;
        s_nsp_init << nsp_init;
        trajectory_ref[0].append(s_nsp_init.str());
      }
      if (trajectory_ref[0] == "" && step == "ComputeTrajectories")
      {
        trajectory_ref[0] = "Trajectory";
      }

 
       
      
      string getRef = "";
      if (step == "DRGEP_Reactions" || step == "DRGEP_Species" || step == "Optim" || step == "ComputeTrajectories")   
         getRef = "outputs/Stochastic/";
      stringstream s_trajectory_ref;
      s_trajectory_ref << trajectory_ref[0];
      getRef.append(s_trajectory_ref.str()).append("_").append(s_nInlet.str()).append(".dat");


      ifstream in_ref(getRef.c_str(), ios::in);
      string line;
      int nbLines_ref = 0;
      while(getline(in_ref, line))
         nbLines_ref++;
      in_ref.close();


     //Get the trajectory of the reduced scheme
      string pathFile;
      stringstream s_nbToKeep;
      s_nbToKeep << nbToKeep; 
      stringstream s_rank;
      s_rank << rank;

      if (step == "DRGEP_Species" || step == "DRGEP_Reactions")
      {
         pathFile = "outputs/Stochastic/Reduced_";
         pathFile.append(s_step.str()).append(s_nbToKeep.str()).append("_").append(s_nInlet.str()).append(".dat");
      }
      else if (step == "QSS")
      {
        pathFile = "Reduced_";
	pathFile.append(s_step.str()).append(s_nbToKeep.str()).append("_").append(s_nInlet.str()).append(".dat");
      }
      else if (step == "Optim")
      {
         pathFile = "analytic_schemes/Ref";
         pathFile.append(s_rank.str()).append("/Reduced_QSS");
	 pathFile.append(s_nbToKeep.str()).append("_").append(s_nInlet.str()).append(".dat");
      }
      else if (step == "ComputeTrajectories") 
      {  
         pathFile = "outputs/Stochastic/Trajectory_"; 
         pathFile.append(s_nInlet.str()).append(".dat"); 
      }

     //Get the reduced scheme
        string pathMech; 

      if (step == "DRGEP_Species")
      {
         pathMech = "outputs/mechanisms/drgepSpecies";
         pathMech.append(s_nbToKeep.str()).append(".xml");
      }
      else if (step == "DRGEP_Reactions")
      {
         pathMech = "outputs/mechanisms/drgepReactions";
         pathMech.append(s_nbToKeep.str()).append(".xml");
      }
      else if (step == "QSS")
        pathMech = "scheme.xml";
      else if (step == "Optim")
      {
         pathMech = "analytic_schemes/Ref";
         pathMech.append(s_rank.str()).append("/scheme.xml");
      }
      else if (step == "ComputeTrajectories")
         pathMech = mech;
      
     

      Read *r = new Read();
      r->Read_species(pathMech.c_str(), listSpecies);

      
      int nsp_ref = listSpecies_ref.size();
      nsp = listSpecies.size();


      //Compter le nombre de ligne 
         ifstream in_check(pathFile.c_str(), ios::in);
                 int nbLines_check = 0;
                        while(getline(in_check, line))
                           nbLines_check++;
                           in_check.close();
      
      vector<vector<double> > Yk_ref (nbLines_ref, vector<double> (nsp_ref, 0.0));
      vector<vector<double> > Yk (nbLines_check, vector<double> (nsp, -1.0));
     
      vector<double> T_ref (nbLines_ref, 0.0);
      vector<double> T (nbLines_check, 0.0);
      

       
  
      ifstream ref(getRef.c_str(), ios::in);
      getline(ref, line);//Skip the first commented line of the trajectory file
        
      double a;
      for (int i=0; i<nbLines_ref-1; i++)
      {
         ref >> a;
         ref >> T_ref[i];
         for (int k=0; k<nsp_ref; k++)
         {
            ref >> Yk_ref[i][k];
         }            

      }
      ref.close();


      ifstream data(pathFile.c_str(), ios::in);
      getline(data, line); //Skip the first commented line of the trajectory file
      for (int i=0; i<nbLines_check-1; i++)
      {
         data >> a;
         data >> T[i];
         for (int k=0; k<nsp; k++)
            data >> Yk[i][k];

      }
      data.close();



      double Coeff = 0.05;
      for (int i=0; i<nbLines_check; i++)
      {
         for (int k=0; k<nsp_ref; k++)
         {
            for (int kbis=0; kbis<nsp; kbis++)
            {

               if (listSpecies_ref[k]->m_Name == listSpecies[kbis]->m_Name)
               {
                  Coeff = 0.05;
                  for (unsigned int kthird=0; kthird<listTargets.size(); kthird++)
                  {
                     if (listSpecies_ref[k]->m_Name == listTargets[kthird])
                        Coeff = 3.0;
                                                
                  }
            
                  if (Yk_ref[i][k] > 1e-06)
                  { 
                     fitness -= Coeff*fabs(Yk[i][kbis]-Yk_ref[i][k])/Yk_ref[i][k]; 
                  }
               }
            }
         }
      }


      if (nbLines_check == 0)
      {
         cout << "Error, scheme not generated " << endl;
         fitness += -100000000;
      }

      if (nbLines_ref == 0)
      {
         cout << "Error, reference trajectory not found" << endl;
         fitness += -100000000;
      }
         }

   fitness *= 1000;
   if (rank == 0)
      {
         if (step == "DRGEP_Species")
         {
             cout << "Reference scheme :  " << nsp_init << " species" << endl;
             cout << "Current scheme :  " << nsp << " species" << endl;
         }

         if (step == "DRGEP_Reactions")
         {
             cout << "Reference scheme :  " <<  listReactions_ref.size() << " reactions" << endl;
             cout << "Current scheme :  " << nbToKeep << " reactions" << endl;
         }

      }



}


void fit_function_1D(string mech_ref, string mech, double& fitness, int nbFlames, vector<string> listTargets, vector<string> trajectory_ref, int nbToKeep, string step, int rank)
{
  
   
   vector<Species_ORCh*> listSpecies_ref;
   vector<Reaction_ORCh*> listReactions_ref;

   Read *r = new Read();
   r->Read_species(mech_ref, listSpecies_ref);
   r->Read_reactions(mech_ref, listReactions_ref);

   fitness = 0.0;

   int nsp_init = listSpecies_ref.size();


   for (int n=0; n<nbFlames; n++)
   {
      stringstream s_step;
      s_step << step;
      if (trajectory_ref[0] == "" && step != "ComputeTrajectories")
      {
        trajectory_ref[0] = "Ref_";
        trajectory_ref[0].append(s_step.str());
        stringstream s_nsp_init;
        s_nsp_init << nsp_init;
        trajectory_ref[0].append(s_nsp_init.str());
      }
      if (trajectory_ref[0] == "" && step == "ComputeTrajectories")
      {
        trajectory_ref[0] = "Premixed_";
      }



      string getRef = "";
      if (step == "DRGEP_Species" || step == "DRGEP_Reactions" || step == "Optim" || step == "ComputeTrajectories")
         getRef = "outputs/Premixed/";
      stringstream s_trajectory_ref;
      s_trajectory_ref << trajectory_ref[n];
      getRef.append(s_trajectory_ref.str()).append(".dat");
     
      ifstream in_ref(getRef.c_str(), ios::in);
      string line;
      int nbLines_ref = 0;
      while(getline(in_ref, line))
         nbLines_ref++;
      in_ref.close();






     //Get the trajectory of the reduced scheme
      string pathFile;
      stringstream s_nbToKeep;
      s_nbToKeep << nbToKeep;
      stringstream s_rank;
      s_rank << rank;
      stringstream s_n;
      s_n << n;

      if (step == "DRGEP_Species" || step == "DRGEP_Reactions")
      {
         pathFile = "outputs/Premixed/Reduced_";
         pathFile.append(s_step.str()).append(s_nbToKeep.str()).append("_").append(s_n.str()).append(".dat");
      }
      else if (step == "QSS")
      {
        pathFile = "Reduced_";
        pathFile.append(s_step.str()).append(s_nbToKeep.str()).append("_").append(s_n.str()).append(".dat");
      }
      else if (step == "Optim")
      {
         pathFile = "analytic_schemes/Ref";
         pathFile.append(s_rank.str()).append("/Reduced_QSS");
         pathFile.append(s_nbToKeep.str()).append("_").append(s_n.str()).append(".dat");
      }
      else if (step == "ComputeTrajectories")
         pathFile = "outputs/Premixed/Premixed_.dat";

 
     //Get the reduced scheme 
      string pathMech;

      if (step == "DRGEP_Species")
      {
         pathMech = "outputs/mechanisms/drgepSpecies";
         pathMech.append(s_nbToKeep.str()).append(".xml");
      }
      else if (step == "DRGEP_Reactions")
      {
         pathMech = "outputs/mechanisms/drgepReactions";
         pathMech.append(s_nbToKeep.str()).append(".xml");
      }
      else if (step == "QSS")
        pathMech = "scheme.xml";
      else if (step == "Optim")
      {
         pathMech = "analytic_schemes/Ref";
         pathMech.append(s_rank.str()).append("/scheme.xml");
      }
      else if (step == "ComputeTrajectories")
         pathMech = mech;





      ifstream data(pathFile.c_str(), ios::in);
      getline(data, line);//Skip the first commented line of the trajectory file
       
      ifstream in_check(pathFile.c_str(), ios::in);
      int nbLines_check = 0;
      while(getline(in_check, line))
         nbLines_check++;
      in_check.close();

      vector<Species_ORCh*> listSpecies;
      r->Read_species(pathMech.c_str(), listSpecies);

      int nsp_ref = listSpecies_ref.size();
      int nsp = listSpecies.size();
      if (rank == 0 && step == "DRGEP_Species")
       {
          cout << "Reference scheme :  " << nsp_ref << " species" << endl;
          cout << "Current scheme :  " << nsp << " species" << endl;
       }

      if (rank == 0 && step == "DRGEP_Reactions")
       {
          cout << "Reference scheme :  " <<  listReactions_ref.size() << " reactions" << endl;
          cout << "Current scheme :  " << nbToKeep << " reactions" << endl;
       }

      vector<vector<double> > Yk_ref (nbLines_ref, vector<double> (nsp_ref, 0.0));
      vector<vector<double> > Yk (nbLines_check, vector<double> (nsp, 0.0));

      vector<double> T_ref (nbLines_ref, 0.0);
      vector<double> T (nbLines_check, 0.0);

      vector<double> U_ref (nbLines_ref, 0.0);
      vector<double> U (nbLines_check, 0.0);

      vector<double> position_ref (nbLines_ref, 0.0);
      vector<double> position (nbLines_check, 0.0);

      ifstream ref(getRef.c_str(), ios::in);
      getline(ref, line); //Skip the first commented line of the trajectory file
      
       for (int i=0; i<nbLines_ref-1; i++)
      {
         ref >> position_ref[i];
         ref >> T_ref[i];
         ref >> U_ref[i];
         for (int k=0; k<nsp_ref; k++)
         {
            ref >> Yk_ref[i][k];
         }
      }
      ref.close();


      for (int i=0; i<nbLines_check-1; i++)
      {
         data >> position[i];
         data >> T[i];
         data >> U[i];
         for (int k=0; k<nsp; k++)
         {
            data >> Yk[i][k];
         }
      }
      data.close();

      double Coeff = 0.01;
      double interpolated_Yk_ref = 0.0;

      for (int k=0; k<nsp_ref; k++)
      {
         for (int kbis=0; kbis<nsp; kbis++)
         {
            if (listSpecies_ref[k]->m_Name == listSpecies[kbis]->m_Name)
            {
               Coeff = 0.5;
               for (unsigned int kthird=0; kthird<listTargets.size(); kthird++)
               {
                  if (listSpecies_ref[k]->m_Name == listTargets[kthird])
                     Coeff = 3.0;
               }

               double fitness_numerateur = 0.0;
               double fitness_denominateur = 0.0;
               for (int i=0; i<nbLines_ref-1; i++)
               {
                  if (T_ref[i] > T_ref[0]+0.02*(T_ref[nbLines_ref-2]-T_ref[2]) && T_ref[i] < T_ref[nbLines_ref-2]-0.02*(T_ref[nbLines_ref-2]-T_ref[2]))
                  {
                     for (int ibis=0; ibis<nbLines_check; ibis++)
                     {
                        if (position[ibis]>position_ref[i] && position[ibis]<position_ref[i+1])
                        {
                           interpolated_Yk_ref = (Yk_ref[i+1][k]-Yk_ref[i][k])*(position[ibis]-position_ref[i])/(position_ref[i+1]-position_ref[i]) + Yk_ref[i][k];




          if (interpolated_Yk_ref > 1e-08)
                           {
                              fitness_numerateur += fabs(Yk[ibis][kbis]-interpolated_Yk_ref)*(position[ibis]-position[ibis-1]);
                              fitness_denominateur += interpolated_Yk_ref*fabs(position[ibis]-position[ibis-1]);
                           }
                        }
                     }
                  }
               }
               if (fitness_denominateur > 1e-15)
                  fitness -= Coeff*(fitness_numerateur/fitness_denominateur);
            }
         }
      }

      double fitness_T_numerateur = 0.0;
      double fitness_T_denominateur = 0.0;
      double fitness_U_numerateur = 0.0;
      double fitness_U_denominateur = 0.0;
      for (int i=0; i<nbLines_ref-1; i++)
      {
         if (T_ref[i] > T_ref[0]+0.02*(T_ref[nbLines_ref-2]-T_ref[2]) && T_ref[i] < T_ref[nbLines_ref-2]-0.02*(T_ref[nbLines_ref-2]-T_ref[2]))
         {
            for (int ibis=0; ibis<nbLines_check; ibis++)
            {
               if (position[ibis]>position_ref[i] && position[ibis]<position_ref[i+1])
               {
                  double interpolated_T_ref = (T_ref[i+1]-T_ref[i])*(position[ibis]-position_ref[i])/(position_ref[i+1]-position_ref[i]) + T_ref[i];
                  double interpolated_U_ref = (U_ref[i+1]-U_ref[i])*(position[ibis]-position_ref[i])/(position_ref[i+1]-position_ref[i]) + U_ref[i];
                  fitness_T_numerateur += fabs(T[ibis]-interpolated_T_ref)*(position[ibis]-position[ibis-1]);
                  fitness_T_denominateur += interpolated_T_ref*fabs(position[ibis]-position[ibis-1]);
                  fitness_U_numerateur += fabs(U[ibis]-interpolated_U_ref)*(position[ibis]-position[ibis-1]);
                  fitness_U_denominateur += interpolated_U_ref*fabs(position[ibis]-position[ibis-1]);
               }
            }
         }
      }
      if (fitness_T_denominateur > 1e-15)
         fitness -= 5*(fitness_T_numerateur/fitness_T_denominateur);



            if (fitness_U_denominateur > 1e-15)
         fitness -= 5*(fitness_U_numerateur/fitness_U_denominateur);



      if (nbLines_check == 0)
      {
         cout << "Error no such file " << endl;
         fitness += -100000000;
      }

   } //end nbFlames loop

   fitness *= 1000;
}















