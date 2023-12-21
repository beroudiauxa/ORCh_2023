#include "computeAutoIgnition.h"

void computeAutoIgnition(string mech, string outputName, string mech_desc, vector<AutoIgnition*> listIgnitions) 
{

   IdealGasMix *mixture = new IdealGasMix(mech, mech_desc);
   int nsp = mixture->nSpecies();

   double *Ym = new double[nsp];
   for (int k=0; k<nsp; k++)
      Ym[k] = 0;
   //double Hm = 0.0;

   //Composition mixture

   IdealGasMix *fuel     = new IdealGasMix(mech,mech_desc);
   IdealGasMix *oxidizer = new IdealGasMix(mech,mech_desc);
   //IdealGasMix *stoichio = new IdealGasMix(mech,mech_desc);

   //Fuel state
   if (listIgnitions[0]->m_Y_fuel != "")
   {
      fuel->setState_TPY(listIgnitions[0]->m_T_fuel, listIgnitions[0]->m_Pressure, listIgnitions[0]->m_Y_fuel);
      cout << "listIgnitions[0]->m_Yf_str " << listIgnitions[0]->m_Y_fuel << endl;
   }

   if (listIgnitions[0]->m_X_fuel != "")
   {
      fuel->setState_TPX(listIgnitions[0]->m_T_fuel, listIgnitions[0]->m_Pressure, listIgnitions[0]->m_X_fuel);
      cout << "listIgnitions[0]->m_Xf_str " << listIgnitions[0]->m_X_fuel << endl;
   }

   //double Hf  = fuel->enthalpy_mass();
   double *Yf = new double[nsp];
   double *Xf = new double[nsp];
   fuel->getMassFractions(Yf);
   fuel->getMoleFractions(Xf);

   //Oxidizer state
   if (listIgnitions[0]->m_Y_oxidizer != "")
   {
      oxidizer->setState_TPY(listIgnitions[0]->m_T_oxidizer, listIgnitions[0]->m_Pressure, listIgnitions[0]->m_Y_oxidizer);
      cout << "listIgnitions[0]->m_Yo_str " << listIgnitions[0]->m_Y_oxidizer << endl;
   }

   if (listIgnitions[0]->m_X_oxidizer != "")
   {
      oxidizer->setState_TPX(listIgnitions[0]->m_T_oxidizer, listIgnitions[0]->m_Pressure, listIgnitions[0]->m_X_oxidizer);
      cout << "listIgnitions[0]->m_Xo_str " << listIgnitions[0]->m_X_oxidizer << endl;
   }

   //double Ho  = oxidizer->enthalpy_mass();
   double *Yo = new double[nsp];
   double *Xo = new double[nsp];
   oxidizer->getMassFractions(Yo);
   oxidizer->getMoleFractions(Xo);

/*
   // get stoechiometry
   double zst;
   findStoechiometryFromCompo(fuel,oxidizer,&zst);
   
   // get fuel and oxidizer mass fraction in fuel
   double Yf_F,Yo_F;
   getYfYo(fuel,Yf_F,Yo_F);
   
   // get fuel and oxidizer mass fraction in oxidizer
   double Yf_O,Yo_O;
   getYfYo(oxidizer,Yf_O,Yo_O);
   
   // Stoichiometry state
   double *Yst = new double[nsp];
   for (int ns=0;ns<nsp;ns++)
      Yst[ns] = zst*Yf[ns]+(1.0-zst)*Yo[ns];
   stoichio->setMassFractions(Yst);
   double Hst = zst*Hf+ (1.0-zst)*Ho;
   stoichio->setState_HP(Hst,listIgnitions[0]->m_Pressure);
   
   // get fuel and oxidizer mass fraction in stoichiometric mixture
   double Yf_st,Yo_st;
   getYfYo(stoichio,Yf_st,Yo_st);
   
   printf("Yf_st: %12.5f   Yo_st: %12.5f \n", Yf_st, Yo_st);
   
   // compute useful quantities
   double s = Yo_st/Yf_st;
   double Z1_O = s*Yf_O - Yo_O;
   double Z1_F = s*Yf_F - Yo_F;
   
   double z = (Yo_O*(listIgnitions[0]->m_Phi-1.0) - Z1_O)/((Z1_F-Z1_O)-(listIgnitions[0]->m_Phi-1.0)*(Yo_F-Yo_O));
   
   // Mixture state
   for (int ns=0;ns<nsp;ns++)
      Ym[ns] = z*Yf[ns]+(1.0-z)*Yo[ns];
   Hm = z*Hf+ (1.0-z)*Ho;
   
   cout << "Mixture fraction " << z << endl;
   cout << "Phi " << listIgnitions[0]->m_Phi << endl;
   for (int k=0; k<nsp; k++)
   {
      cout << "Species " << k << "  " << Ym[k] << endl;
   }
*/









  ofstream i_delay ("outputs/i_delay.dat");

  //double t_phi_init = 0.001;

  vector<double> Phi_table (8, 0.0);
  Phi_table[0] = 0.5;
  Phi_table[1] = 0.55;
  Phi_table[2] = 0.6;
  Phi_table[3] = 0.65;
  Phi_table[4] = 0.7;
  Phi_table[5] = 0.75;
  Phi_table[6] = 0.8;
  Phi_table[7] = 0.85;
  Phi_table[8] = 0.9;
  Phi_table[9] = 0.95;
  Phi_table[10] = 1.0;

  for (int nbPhi=0; nbPhi<11; nbPhi++)
  {


   //Temporary stuff
     //double t_phi = t_phi_init*pow(1.1,nbPhi);
     double t_phi = Phi_table[nbPhi];
     double Y_O2_fuel = 0.194221;
     double Y_N2_fuel = 0.589443;
     double Y_H2O_fuel = 0.00211403;
     double Y_CH4_fuel = 0.214222;
     double T_fuel = 320;
   
     double Y_O2_coflow = 0.142234;
     double Y_N2_coflow = 0.757491;
     double Y_H2O_coflow = 0.100097;
     double Y_CH4_coflow = 0.00017827;
     double T_coflow = 1350;
   
     double Z_1 = (4*Y_CH4_fuel-t_phi*Y_O2_fuel)/(t_phi*Y_O2_coflow-4*Y_CH4_coflow)/(1+(4*Y_CH4_fuel-t_phi*Y_O2_fuel)/(t_phi*Y_O2_coflow-4*Y_CH4_coflow));
     double Z_2 = 1-Z_1;
   
     double Y_O2 = Y_O2_fuel*Z_2+Y_O2_coflow*Z_1;
     double Y_H2O = Y_H2O_fuel*Z_2+Y_H2O_coflow*Z_1;
     double Y_CH4 = Y_CH4_fuel*Z_2+Y_CH4_coflow*Z_1;
     double Y_N2 = Y_N2_fuel*Z_2+Y_N2_coflow*Z_1;
   
   
      for (int k=0; k<nsp; k++)
      {
         if (mixture->speciesName(k) == "CH4")
            Ym[k] = Y_CH4;
   
         if (mixture->speciesName(k) == "O2")
            Ym[k] = Y_O2;
   
         if (mixture->speciesName(k) == "H2O")
            Ym[k] = Y_H2O;
   
         if (mixture->speciesName(k) == "N2")
            Ym[k] = Y_N2;
      }
      cout << "Z1 " << Z_1 << endl;
      cout << "Z2 " << Z_2 << endl;
      
      double Tm = T_fuel*Z_2+T_coflow*Z_1;









   





   for (unsigned int n=0; n<listIgnitions.size(); n++)
   {
      int nbLines = listIgnitions[n]->m_max_t/listIgnitions[n]->m_delta_t;
      cout << "nbLines " << nbLines << endl;

      double *Y = new double[nsp];
      double H = 0.0;
      double T = 0.0;

      //---Particles store---
      vector<vector<double> > Y_Store(nbLines, vector<double>(nsp));
      vector<double> H_Store(nbLines, 0.0);
      vector<double> T_Store(nbLines, 0.0);
      vector<double> t_Store(nbLines, 0.0);

      IdealGasMix *IgnitionMixture = new IdealGasMix(mech, mech_desc);
      int nsp = IgnitionMixture->nSpecies();
      double t = listIgnitions[n]->m_delta_t;


/*
      if (listIgnitions[n]->m_X_Species != "")
      {
         cout << "Set the mole fraction of inlet " << n << endl;
         IgnitionMixture->setState_TPX(listIgnitions[n]->m_Temperature,
                                    listIgnitions[n]->m_Pressure,
                                    listIgnitions[n]->m_X_Species);
      }
      else if (listIgnitions[n]->m_Y_Species != "")
      {
         cout << "Set the mass fraction of inlet " << n << endl;
         IgnitionMixture->setState_TPY(listIgnitions[n]->m_Temperature,
                                    listIgnitions[n]->m_Pressure,
                                    listIgnitions[n]->m_Y_Species);
      }
*/

      //H  = IgnitionMixture->enthalpy_mass();
      //IgnitionMixture->getMassFractions(Y);
      
      IgnitionMixture->setMassFractions(Ym);
      //IgnitionMixture->setState_HP(Hm, listIgnitions[n]->m_Pressure);
      IgnitionMixture->setState_TP(Tm, listIgnitions[n]->m_Pressure);

      t = listIgnitions[n]->m_delta_t;
   
      string path = "outputs/Ignite";
      stringstream s_Phi;
      s_Phi << t_phi;
      //path.append(s_Phi.str()).append("_").append("16S.dat");
      path.append(s_Phi.str()).append("_").append("gri30.dat");
      
      
      cout << "Path " << path << endl;
      ofstream Ignite (path.c_str());

      for (int i=0; i<nbLines-1; i++)
      {
         //IgnitionMixture->setMassFractions(Y);
         //IgnitionMixture->setState_HP(H, listIgnitions[n]->m_Pressure);

         ConstPressureReactor reac;
         reac.insert(*IgnitionMixture);
         ReactorNet sim;
         sim.addReactor(reac);
      
         sim.advance(listIgnitions[n]->m_delta_t);

         H  = IgnitionMixture->enthalpy_mass();
         T = IgnitionMixture->temperature();
         IgnitionMixture->getMassFractions(Y);
      
         Ignite << t << "  " << T << "  ";
         for (int k=0; k<nsp; k++)
            Ignite << Y[k] << "  ";
         Ignite << endl;





         t = t + listIgnitions[n]->m_delta_t;
         t_Store[i+1] = t;
         T_Store[i+1] = T;

      }
      Ignite.close();

      for (int k=0; k<nsp; k++)
      {
         cout << "Species " << k << "  " << Y[k] << endl;
      }

      bool done = true;
      for (int i=0; i<nbLines; i++)
      {
         if (T_Store[i] > T_Store[1]+0.5*(T_Store[nbLines-1]-T_Store[1]) && done)
         {
            cout << "Ignition delay " << t_Store[i] << "  phi " << t_phi << " T_middle " << T_Store[i] << "  Tmin" << T_Store[0] << "  Tmax" << T_Store[nbLines-1] << endl;
            i_delay << t_phi << "  " << t_Store[i] << endl;
            done = false;
         }
      }


   }
   
}

i_delay.close();




}



