#include "gnuplot.h"
#include <fstream>
#include <sstream>


//---gnuplot---

void Script_gnuplot (string step, string initial_mech, vector<string> speciesToPlot, string configuration, string mech_desc, int nbToKeep,bool plot_U, bool plot_T, vector<string> trajectory_ref, string mech_ref, int nbInlets_flames, vector<string> listTarget,string outputSchemeName, int rank)
{ 


//Finding the species and reactions lists of the initial scheme
   string choice_mech_ref;

   if (trajectory_ref[0].empty())
      choice_mech_ref = initial_mech;
   else
     choice_mech_ref = mech_ref; 


   IdealGasMix *mixture  = new IdealGasMix(choice_mech_ref,mech_desc);
   
   vector<Species_ORCh*> listSpecies_ref;
   vector<Reaction_ORCh*> listReactions_ref;
  
   Read *r = new Read();
   r->Read_species(choice_mech_ref, listSpecies_ref);
   r->Read_reactions(choice_mech_ref, listReactions_ref);

   int nsp_ref = mixture->nSpecies();  
   int nreac_ref = mixture->nReactions();


//Step chosen
   stringstream s_step;
   s_step << step;
   int nb_ref;  

   if (step == "DRGEP_Reactions")
      nb_ref = nreac_ref;
   else
      nb_ref = nsp_ref ;        


//Get reduced scheme 
  
   vector<Species_ORCh*> listSpecies;
//   string outputSchemeName;  


/*   if (step == "QSS")
      outputSchemeName = "scheme.xml";
   else if (step == "ComputeTrajectories")
{      outputSchemeName = initial_mech;
cout << "Reduced scheme for plot : " << initial_mech << endl; }
   else
   {
      outputSchemeName= "./outputs/mechanisms/drgep";
     
      if (step == "DRGEP_Species")
         outputSchemeName.append("Species");
      else if (step == "DRGEP_Reactions")
         outputSchemeName.append("Reactions");

      stringstream s_nbToKeep;
      s_nbToKeep << nbToKeep;
      outputSchemeName.append(s_nbToKeep.str()).append(".xml");
   }
*/

      Read *r2 = new Read();
      r2->Read_species(outputSchemeName, listSpecies);
 
//If the user doesnt define his own ref scheme, the default one is the initial one
    if (trajectory_ref[0]=="" && step != "ComputeTrajectories")
    {
       if (rank ==0)
          cout << endl << endl << "Default reference trajectory used." <<  endl << endl;
       trajectory_ref[0] = "Ref_";
       trajectory_ref[0].append(s_step.str());
       stringstream s_nb_ref;
       s_nb_ref << nb_ref;
       trajectory_ref[0].append(s_nb_ref.str());             
    }
    if (trajectory_ref[0]=="" && step == "ComputeTrajectories")
    {   
       if ( configuration == "MultipleInlet")
          trajectory_ref[0] = "Trajectory";
       if (configuration == "PremixedFlames")
          trajectory_ref[0] = "Premixed_";
    }
   
    vector <int> nCol; //Number of column of the ref trajectories
    vector <int> nColR; //Number of column of the reduced scheme

//Get the column numbers 
    
    int nb(0);
    
    if (configuration == "MultipleInlet")
       nb = 3;
    if (configuration == "PremixedFlames")
       nb = 4;
 

    // Get Species to Plot
    for (int i=0; i< speciesToPlot.size();i++)
    {

       for (int k=0; k< nsp_ref; k++)
       {
          if (listSpecies_ref[k]->m_Name == speciesToPlot[i])
               nCol.push_back(k + nb);
       }

          for (int k=0; k< listSpecies.size(); k++)
          {
              if ( listSpecies[k]-> m_Name  == speciesToPlot[i])
                 nColR.push_back(k + nb);
          }
            
    }   



 
   
   string gnuplotName;

   if (step != "QSS")
   {
      gnuplotName = "./outputs/";
   
      if (configuration == "PremixedFlames")
         gnuplotName.append("Premixed/makeGnu.gnu"); 

      if (configuration == "MultipleInlet")
         gnuplotName.append("Stochastic/makeGnu.gnu");
   }
   else if (step == "QSS")  
    {  gnuplotName = "makeGnu.gnu";
   }


   ofstream Gnuplot(gnuplotName.c_str());

  	   Gnuplot << "#File to print ..." << endl;
  	   Gnuplot << "" << endl;
           Gnuplot << "set encoding iso_8859_1" << endl;
 	   Gnuplot << "set terminal postscript portrait color noenhanced \"Times-Roman\" 16" << endl;
  	   Gnuplot << "" << endl;
	   Gnuplot << "set style line 1  linetype 1 linecolor rgb \"#000000\" linewidth 6 pointtype 2 pointsize 0.6" << endl;
	   Gnuplot << "set style line 3  linetype 2 linecolor rgb \"#CC0033\" linewidth 6 pointtype 8 pointsize 0.8" << endl;
	   Gnuplot << "set macros" << endl;
	   Gnuplot << "TMARGIN = \"set tmargin at screen 0.90; set bmargin at screen 0.60\"" << endl;
	   Gnuplot << "BMARGIN = \"set tmargin at screen 0.50; set bmargin at screen 0.20\"" << endl;
	   Gnuplot << "" << endl;
	   Gnuplot << "LMARGIN = \"set lmargin at screen 0.3; set rmargin at screen 0.80\"" << endl;
	   Gnuplot << "RMARGIN = \"set lmargin at screen 0.60; set rmargin at screen 0.90\"" << endl;
	   Gnuplot << "" << endl;
	   Gnuplot << "" << endl;
  	   Gnuplot << "@BMARGIN" << endl;
  	   Gnuplot << "@LMARGIN" << endl;
	   Gnuplot << "" << endl;
	   Gnuplot << "set key at graph 0.98,1.2" << endl;
	   Gnuplot << "set key font \"0,14\"" << endl;
	   Gnuplot << "set key spacing 0.6" << endl;
	   Gnuplot << "" << endl;



   if (configuration == "MultipleInlet")
   {

     //Calcul of the fitness function
     double fitness = 0;  

    fit_function_0D(choice_mech_ref, initial_mech, fitness, nbInlets_flames, listTarget, trajectory_ref, nbToKeep, step, rank);

    if (rank ==0)
    {
       cout << endl << "fitness between the reference and the current reduced scheme : " << fitness  << endl;
    }
	   Gnuplot << "set label \"fitness : " << fitness  <<"\" at graph  0.03,0.95" << endl;


      for(int inlet=0; inlet < nbInlets_flames; inlet++)
      {
           
	   Gnuplot << "set output \""<< step << nbToKeep  <<"_Inlet" << inlet << ".eps\"" << endl;

           Gnuplot << "set xlabel \"Time (ms)\" " << endl;


         if(plot_T)
         { 

           Gnuplot << "" << endl;
           Gnuplot << "" << endl;
           Gnuplot << "set ylabel \"T(K)\"" << endl;
           Gnuplot << "" << endl;
	   Gnuplot << "plot \""<< trajectory_ref[0] << "_" << inlet << ".dat\" using ($1*1000):2 with lines ls 1 title '"<< trajectory_ref[0]  <<" Inlet "<< inlet << "',\\"  << endl;
           if (step == "ComputeTrajectories")
           Gnuplot << "     \"Trajectory_" << inlet << ".dat\" using ($1*1000):2 with lines ls 3 title '"<< step  << nbToKeep << " mech Inlet "<< inlet  <<" '"<< endl;
           else
           Gnuplot << "     \"Reduced_" << step << nbToKeep << "_" << inlet << ".dat\" using ($1*1000):2 with lines ls 3 title 'Reduced "<< step  << nbToKeep << " mech Inlet "<< inlet  <<" '"<< endl;
           Gnuplot << "" << endl;

         }



         for(int i=0 ; i< nColR.size();i++)
         {
           Gnuplot << "" << endl;
           Gnuplot << "" << endl;
           Gnuplot << "set ylabel \"Y("<< speciesToPlot[i] << ")\"" << endl;
           Gnuplot << "" << endl;
           Gnuplot << "plot \""<< trajectory_ref[0] <<"_" << inlet << ".dat\" using ($1*1000):"<< nCol[i] << " with lines ls 1 title '"<< trajectory_ref[0]  <<" Inlet "<< inlet << "',\\"  << endl;
           if (step == "ComputeTrajectories")
           Gnuplot << "     \"Trajectory_" << inlet <<".dat\" using ($1*1000):"<< nColR[i] << " with lines ls 3 title '"<< step << nbToKeep << " mech Inlet "<< inlet << "'"<< endl;
           else
           Gnuplot << "     \"Reduced_" << step << nbToKeep<< "_" << inlet<< ".dat\" using ($1*1000):"<< nColR[i] << " with lines ls 3 title 'Reduced "<< step << nbToKeep << " mech Inlet "<< inlet << "'"<< endl;
           Gnuplot << "" << endl;

         }      

      } 

   }  //endif DRGEP multiple inlet  

   if (configuration == "PremixedFlames")
   {  



 
     //Calcul of the fitness function
     double fitness = 0;
    
     // Commented by Huu-Tri Nguyen - 2020.02.19
   //HT  fit_function_1D( choice_mech_ref, initial_mech, fitness, nbInlets_flames+1 ,  speciesToPlot,  trajectory_ref, nbToKeep, step, rank);

    // Fixed by Huu-Tri Nguyen - 2020.02.19 - Calculate fitness on listTarget (not speciesToPlot)
    fit_function_1D( choice_mech_ref, initial_mech, fitness, nbInlets_flames+1 ,  listTarget,  trajectory_ref, nbToKeep, step, rank);  

    if (rank ==0)
    {
       cout << endl << "fitness between the reference and the current reduced scheme : " << fitness  << endl << endl;
    }

   for (int n=0;n<nbInlets_flames+1;n++)
   {   

//if (configuration == "QSS")
//{
//   string dataSL = "dataSL.dat";
//   ofstream dat(dataSL.c_str());
//
//   dat << "#1:phi  2:Sl ref  3:Sl" << endl;
// 
//}


	   Gnuplot << "set label \"fitness :" << fitness  <<"\" at graph  0.03,0.95" << endl;

           Gnuplot << "set output \""<< step << nbToKeep << "_Flame" << n  << ".eps\"" << endl;
      

           Gnuplot << "set xrange[-10:10]" << endl;
           Gnuplot << "set xlabel \"position (mm)\" offset 0,0 " << endl;
           Gnuplot << "" << endl;
      if(plot_T)
      {
           Gnuplot << "" << endl;
           Gnuplot << "" << endl;
           Gnuplot << "set ylabel \"T(K)\"" << endl;
           Gnuplot << "" << endl;
           Gnuplot << "plot \""<< trajectory_ref[n] << ".dat\" using ($1*1000):2 with lines ls 1 title '"<< trajectory_ref[n]  << "',\\"  << endl;
           if (step == "ComputeTrajectories")
           Gnuplot << "     \"Premixed_.dat\" using ($1*1000):2 with lines ls 3 title 'Reduced " << nbToKeep << " species mech '" << endl << endl;
           else
           Gnuplot << "     \"Reduced_" << step  << nbToKeep << "_" << n  << ".dat\" using ($1*1000):2 with lines ls 3 title 'Reduced " << nbToKeep << " species mech '" << endl << endl;
      }

      if(plot_U)
      {
           Gnuplot << "" << endl;
           Gnuplot << "" << endl;
           Gnuplot << "set ylabel \"U (m/s)\"" << endl;
           Gnuplot << "" << endl;
           Gnuplot << "plot \""<< trajectory_ref[n] << ".dat\" using ($1*1000):3 with lines ls 1 title '"<< trajectory_ref[n]  <<"',\\"  << endl;
           if (step == "ComputeTrajectories")
           Gnuplot << "     \"Premixed_.dat\" using ($1*1000):3 with lines ls 3 title 'Reduced " << step  << nbToKeep << " mech '" << endl << endl;
           else
           Gnuplot << "     \"Reduced_" << step << nbToKeep << "_" << n  << ".dat\" using ($1*1000):3 with lines ls 3 title 'Reduced " << step  << nbToKeep << " mech '" << endl << endl;
      }

      for(int i=0 ; i< nColR.size();i++)
      {
    
           Gnuplot << "" << endl;
           Gnuplot << "set ylabel \"Y("<< speciesToPlot[i] << ")\"" << endl;

           Gnuplot << "plot \""<< trajectory_ref[n] << ".dat\" using ($1*1000):"<< nCol[i] << " with lines ls 1 title '"<< trajectory_ref[n]  <<"',\\"  << endl;
           if (step == "ComputeTrajectories")
           Gnuplot << "     \"Premixed_.dat\" using ($1*1000):"<< nColR[i] <<  " with lines ls 3 title 'Reduced " << step  << nbToKeep << " mech '" << endl << endl;
           else
           Gnuplot << "     \"Reduced_" << step << nbToKeep << "_" << n  << ".dat\" using ($1*1000):"<< nColR[i] <<  " with lines ls 3 title 'Reduced " << step  << nbToKeep << " mech '" << endl << endl;
      }             
   
   
   } //end loop nbFlame             

//   Gnuplot << "plot \"dataSL.dat\" u 1:2 with points ls 1 title 'ref',\\" << endl;
//   Gnuplot << "     \"dataSL.dat\" u 1:3 with lines ls 3 title 'reduced',\\" << endl;

} // endif conf premixed flames 


Gnuplot.close();

}  // close function definition


