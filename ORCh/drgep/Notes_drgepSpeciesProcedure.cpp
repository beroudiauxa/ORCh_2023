//----------------------------------------
//   <drgep_0D_species> 
//----------------------------------------
//   DESCRIPTION:
//         Run the DRGEP analysis for 0D (AutoIgnition or MultipleInlet) regimes to sort species according to their importance to the creation or destruction of the targets
//   IN: 
//         *mixture: Pointer to the composition of the mixture to get forward and reverse reaction rates
//         Targets: List of targets (species), (Targets[k] == true) if target and (Targets[k] == false) if not
//         n: Number of the inlet
//         time: Time for the integration
//   OUT:
//         R_AD_Trajectories: DRGEP matrix of interactions between species
//----------------------------------------
void drgep::drgep_0D_species(IdealGasMix *mixture, vector<bool> Targets, vector<vector<double> > &R_AD_Trajectories, int n, double time) const
{

   int nreac = mixture->nReactions();

   double* forwardRates = new double[nreac];
   double* reverseRates = new double[nreac];
   mixture->getFwdRatesOfProgress(forwardRates);
   mixture->getRevRatesOfProgress(reverseRates);

   vector<double> fwdRates (nreac, 0.0);
   vector<double> revRates (nreac, 0.0);

   for (int j=0; j<nreac; j++)
   {
      fwdRates[j] = forwardRates[j];
      revRates[j] = reverseRates[j];
   }

   //To store the data files that provide the interactions between species and get DRGEP graphs
   string outputName = "output/rAB_inlet";
   stringstream s_nInlet;
   stringstream s_time;
   s_nInlet << n;
   s_time << time;
   outputName.append(s_nInlet.str()).append("_").append(s_time.str()).append(".dat");

   drgep_species(mixture, Targets, R_AD_Trajectories, fwdRates, revRates, outputName); 
}



//----------------------------------------
//   <drgep_species> 
//----------------------------------------
//   DESCRIPTION:
//         Run the DRGEP analysis for species
//   IN: 
//         *mixture: Pointer to the composition of the mixture to get forward and reverse reaction rates
//         Targets: List of targets (species), (Targets[k] == true) if target and (Targets[k] == false) if not
//         fwdRates: Reactions forward rates (unit: ...)
//         revRates: Reactions reverse rates (unit: ...)
//         outputName: Name of the file in which the species inter-relations are written in order to plot relation graphs
//   OUT:
//         R_AD_Trajectories: DRGEP matrix of interactions between species
//----------------------------------------
void drgep::drgep_species(IdealGasMix *mixture, vector<bool> Targets, vector<vector<double> > &R_AD_Trajectories, vector<double> fwdRates, vector<double> revRates, string outputName) const
{

   int nreac = mixture->nReactions(); //Number of reactions
   int nsp = mixture->nSpecies(); //Number of chemical species
   int nel = mixture->nElements(); //Number of elements (typically: C, H, O, N, Ar)

   vector<double> Production_minus_consumption_k (nsp, 0.0);
   vector<double> Production_Atom_A (nel, 0.0);

   for (int k=0; k<nsp; k++)
   {
      double omega_k_prod = 0.0;
      double omega_k_cons = 0.0;
      for (int j=0; j<nreac; j++)
      {
         omega_k_prod += mixture->productStoichCoeff(k,j)*fwdRates[j]
                        +mixture->reactantStoichCoeff(k,j)*revRates[j];
         omega_k_cons += mixture->reactantStoichCoeff(k,j)*fwdRates[j]
                        +mixture->productStoichCoeff(k,j)*revRates[j];
      }

      Production_minus_consumption_k[k] = abs(omega_k_prod-omega_k_cons);

      for (int Atom=0; Atom<nel; Atom++)
            Production_Atom_A[Atom] += mixture->nAtoms(k,Atom)*(omega_k_prod);
   }

   vector<vector<double> > scaling_Atom_Target (nel, vector<double> (nsp));
   vector<double> scaling_Target (nsp, 0.0);

   for (int k=0; k<nsp; k++)
   {
      for (int Atom=0; Atom<nel; Atom++)
      {
         scaling_Atom_Target[Atom][k] = mixture->nAtoms(k,Atom)*Production_minus_consumption_k[k]/Production_Atom_A[Atom];

         if (scaling_Atom_Target[Atom][k] > scaling_Target[k])
            scaling_Target[k] = scaling_Atom_Target[Atom][k];
      }
   }


   vector<vector<double> > r_AB_sup (nsp, vector<double>(nsp));
   vector<vector<double> > r_AB_inf (nsp, vector<double>(nsp));
   vector<vector<double> > r_AB (nsp, vector<double>(nsp));

   //r_AB quantifies the direct inter-relation between species 
   for (int ka=0; ka<nsp; ka++)
   {
      for (int kb=0; kb<nsp; kb++)
      {
         r_AB_sup[ka][kb] = 0.0;
         r_AB_inf[ka][kb] = 0.0;
         for (int j=0; j<nreac; j++)
         {
            if (mixture->reactantStoichCoeff(kb,j) != 0.0 || mixture->productStoichCoeff(kb,j) != 0.0)
            {
               r_AB_sup[ka][kb] += abs((mixture->reactantStoichCoeff(ka,j)+mixture->productStoichCoeff(ka,j))*(fwdRates[j]-revRates[j]));
            }
            r_AB_inf[ka][kb] += abs((mixture->reactantStoichCoeff(ka,j)+mixture->productStoichCoeff(ka,j))*(fwdRates[j]-revRates[j]));
         }

         r_AB[ka][kb] = r_AB_sup[ka][kb]/r_AB_inf[ka][kb];

         if (r_AB_inf[ka][kb] == 0.0)
         {
            r_AB[ka][kb] = 0.0;
         }
      }
   }

//----------Print the species inter-relations diagrams----------//

   ofstream rAB(outputName.c_str());
   for (int ka=0; ka<nsp; ka++)
   {
      rAB << mixture->speciesName(ka) << "  ";
   }
   rAB << endl;

   for (int ka=0; ka<nsp; ka++)
   {
      for (int kb=0; kb<nsp; kb++)
      {
         rAB << r_AB[ka][kb] << "  ";
      }
      rAB << endl;
   }
   rAB.close();

//--------------------------------------------------------------//

   vector<vector<double> > r_AD_intermediate (nsp, vector<double>(nsp));
   vector<vector<double> > r_AD (nsp, vector<double>(nsp));

   int nb_interaction_level = 3;

   for (int interaction_level=2; interaction_level<2+nb_interaction_level; interaction_level++)
   {
      //First copy last solution (r_AB if it is the first DRGEP step or r_AB_DRGEP if the EP has already been applied)
      for (int ka=0; ka<nsp; ka++)
      {
         for (int kb=0; kb<nsp; kb++)
         {
            if (interaction_level == 2)
               r_AD_intermediate[ka][kb] = r_AB[ka][kb];
            else
               r_AD_intermediate[ka][kb] = r_AD[ka][kb];
         }
      }

      for (int ka=0; ka<nsp; ka++)
      {
         for (int kb=0; kb<nsp; kb++)
         {
            for (int ki=0; ki<nsp; ki++) //i : intermediate species
            {
               if (r_AD[ka][kb] < r_AD_intermediate[ka][ki]*r_AB[ki][kb])
                  r_AD[ka][kb] = r_AD_intermediate[ka][ki]*r_AB[ki][kb];
            }
         }
      }

   } //end interaction

   for (int ka=0; ka<nsp; ka++)
   {
      if (Targets[ka])
      {
         for (int kb=0; kb<nsp; kb++)
         {
            if (R_AD_Trajectories[ka][kb] < scaling_Target[ka]*r_AD[ka][kb])
            {
               R_AD_Trajectories[ka][kb] = scaling_Target[ka]*r_AD[ka][kb];
            }
         }
      }
   }
}

