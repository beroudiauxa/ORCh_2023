#include "outputs.h"

//----------------------------------------
//   <OutputDeterministicTrajectories> 
//----------------------------------------
//   DESCRIPTION:
//         
//   IN: 
//         
//   OUT:
//
//----------------------------------------
void OutputDeterministicTrajectories (int nInlets, int nbIterations, vector<Species_ORCh*> listSpecies, string outputName, vector<vector<vector<double> > > Ym_Trajectories, vector<vector<double> > T_Trajectories, vector<double> time_store)
{

   int nsp = listSpecies.size();
   for (int n=0; n<nInlets-1; n++)
   {
      string Name = outputName;
      stringstream s_nInlet;
      s_nInlet << n;
      Name.append(s_nInlet.str()).append(".dat");
      ofstream Trajectory(Name.c_str());
      Trajectory << "#1:time(s)  "; 
      Trajectory << "2:T(K)";
      for (int k=0; k<nsp; k++)
      {
         Trajectory << "  " << k+3 << ":Y_" << listSpecies[k]->m_Name;
      }
         Trajectory << endl;
      for (int i=0; i<nbIterations; i++)
      {
         Trajectory << time_store[i] << "  ";
         Trajectory <<  T_Trajectories[n][i];
         for (int k=0; k<nsp; k++)
         {
            Trajectory << "  " <<  Ym_Trajectories[n][i][k];
         }
         Trajectory << endl;
      }
      Trajectory.close();
   }
}


void output_datas (string outputName, vector<vector<double> > Y, vector<double> T , vector<double> U, vector<double> position, double position_max_wdot, IdealGasMix* mixture)
{

      int nsp = mixture->nSpecies();

      string Name = outputName;
      Name.append(".dat");
      ofstream Trajectory(Name.c_str());

    
      Trajectory << "#1:position(m)  "; 
      Trajectory << "2:T(K)  ";
      Trajectory << "3:U(m/s)  ";
      for (int k=0; k<nsp; k++)
      {
         Trajectory << "(" << k+4 << ":" << mixture->speciesName(k) << ") ";
      }  
      Trajectory << endl;

      for (int ino=0; ino<position.size(); ino++)
      {
         Trajectory << position[ino]-position_max_wdot << "  " << T[ino] << "  " << U[ino] << "  ";
         for (int k=0; k<nsp; k++)
         {
            Trajectory << Y[ino][k] << "  ";
         }
         Trajectory << endl;
      }

      Trajectory.close();

   }


void plotFitness(string path)
{

           path.append("make.gnu");

           ofstream Gnuplot(path.c_str());

           Gnuplot << "#File to print ..." << endl;
           Gnuplot << "" << endl;
           Gnuplot << "set encoding iso_8859_1" << endl;
           Gnuplot << "set terminal postscript portrait color noenhanced \"Times-Roman\" 12" << endl;
           Gnuplot << "" << endl;
           Gnuplot << "set style line 2000 linetype 2 linecolor rgb \"#000000\" linewidth 1.5 pointtype 5 pointsize 0.5" << endl;
           Gnuplot << "set style line 1000 linetype 2 linecolor rgb \"#00c000\" linewidth 1.5 pointtype 5 pointsize 0.5" << endl;
           Gnuplot << "" << endl;
           Gnuplot << "" << endl;
           Gnuplot << "set macros" << endl;
           Gnuplot << "TMARGIN = \"set tmargin at screen 0.90; set bmargin at screen 0.60\"" << endl;
           Gnuplot << "LMARGIN = \"set lmargin at screen 0.25; set rmargin at screen 0.80\"" << endl;
           Gnuplot << "" << endl;
           Gnuplot << "set key at graph 0.5,-0.2" << endl;
           Gnuplot << "set key font \"0,8\"" << endl;
           Gnuplot << "set key spacing 0.8" << endl;
           Gnuplot << "" << endl;
           Gnuplot << "set output \"fitnessEvolution.eps\"" << endl;
           Gnuplot << "" << endl;
           Gnuplot << "set xlabel \"generation\" " << endl;
           Gnuplot << "" << endl;
           Gnuplot << "@TMARGIN" << endl;
           Gnuplot << "@LMARGIN" << endl;
           Gnuplot << "" << endl;
           Gnuplot << "set ylabel \"fitness\"" << endl;
           Gnuplot << "plot \"fitness.dat\" using 1:2 with linespoints ls 1000 title 'mean fitness' ,\\" << endl;
           Gnuplot << "     \"fitness.dat\" using 1:3 with linespoints ls 2000 title 'best fitness'" << endl;
           Gnuplot << " "<< endl;
           Gnuplot << "" << endl;
//           Gnuplot << "@TMARGIN" << endl;
//           Gnuplot << "@LMARGIN" << endl;
//           Gnuplot << "" << endl;
//           Gnuplot << "set yrange [0.5:1]" << endl;
//           Gnuplot << "set ylabel \"variation des Aj\"" << endl;
//           Gnuplot << "" << endl;
//           Gnuplot << "plot ";
//           for (int j=0; j<listReactions.size();j++)          
//               Gnuplot << "    \"dataAj.dat\" using 1:"<< listReactions.size()*2 + j+2 <<" with linespoints ls "<< j+1 << " title 'reaction "<< listReactions[j]-> m_equation << "' ,\\" << endl;
//           
//                                 Gnuplot << " "<< endl;
//                                 Gnuplot << "" << endl;
//
           Gnuplot.close();

}
