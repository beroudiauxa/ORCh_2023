#include "computePremixedFlames.h"


void computePremixedFlames(string mech, string outputName, string mech_desc, vector<PremixedFlames*> listFlames, vector<bool> Targets, string step, vector<bool> &SpeciesIntoReactants, vector<vector<double> > &R_AD_Premixed, vector<vector<double> > &max_j_on_Target, vector<vector<double> > &max_jf_on_Target, vector<vector<double> > &max_jr_on_Target, vector<vector<double> > &QSS_Criteria, bool Proceed_DRGEP_analysis)
{


   IdealGasMix* mixture  = new IdealGasMix(mech,mech_desc);
   //FreeFlame* flow = new FreeFlame(mixture);
   StFlow* flow = new StFlow(mixture);

   flow->setFreeFlow();

   Inlet1D* inlet = new Inlet1D;
   Outlet1D* outlet = new Outlet1D;

   vector<Domain1D*> domains;
   domains.push_back(inlet);
   domains.push_back(flow);
   domains.push_back(outlet);
   int flowdomain=1; // 0 is inlet, 1 is flow, 2 is outlet
   Sim1D* flame = new Sim1D(domains);

   //----------Set properties----------//
   flow->solveEnergyEqn();
   Transport* trmix = newTransportMgr("Mix", mixture);
   flow->setTransport(*trmix);
   flow->setKinetics(*mixture);
   
   cout << "PREMIXED FLAME:-------------------------------------------------------------------------------------------------" << endl;
   cout << "            Reading initial flame \"" << listFlames[0]->m_Initial_Flame<< "\" with description \"st_flame\" ----------> " << flush;
   flame->restore(listFlames[0]->m_Initial_Flame, "st_flame");
   cout << "OK" << endl;
   //cout << "----------------------------------------------------------------------------------------------------------------------------------" << endl;
   
   //----------Number of species----------//
   int nsp = mixture->nSpecies();
   int nreac = mixture->nReactions();

   double *Ym = new double[nsp];
   double Hm = 0.0;

   //----------Composition mixture----------//
   
   // Initialize Fuel and Oxidizer
   IdealGasMix *fuel     = new IdealGasMix(mech,mech_desc);
   IdealGasMix *oxidizer = new IdealGasMix(mech,mech_desc);
   IdealGasMix *stoichio = new IdealGasMix(mech,mech_desc);
   // fuel state
   //
   if (listFlames[0]->m_Yf_str != "")
   {
      fuel->setState_TPY(listFlames[0]->m_T_fuel, listFlames[0]->m_Pressure, listFlames[0]->m_Yf_str);
     // cout << "listFlames[0]->m_T_fuel " << listFlames[0]->m_T_fuel << endl;
   }

   if (listFlames[0]->m_Xf_str != "")
   {
      fuel->setState_TPX(listFlames[0]->m_T_fuel, listFlames[0]->m_Pressure, listFlames[0]->m_Xf_str);
      // cout << "listFlames[0]->m_Xf_str " << listFlames[0]->m_Xf_str << endl;
   }

   double Hf  = fuel->enthalpy_mass();
   double *Yf = new double[nsp];
   double *Xf = new double[nsp];
   fuel->getMassFractions(Yf);
   fuel->getMoleFractions(Xf);

   // oxidizer state
   //
   if (listFlames[0]->m_Yo_str != "")
   {
      oxidizer->setState_TPY(listFlames[0]->m_T_oxidizer, listFlames[0]->m_Pressure, listFlames[0]->m_Yo_str);
     // cout << "listFlames[0]->m_Yo_str " << listFlames[0]->m_Yo_str << endl;
   }

   if (listFlames[0]->m_Xo_str != "")
   {
      oxidizer->setState_TPX(listFlames[0]->m_T_oxidizer, listFlames[0]->m_Pressure, listFlames[0]->m_Xo_str);
     // cout << "listFlames[0]->m_Xo_str " << listFlames[0]->m_Xo_str << endl;
   }

   double Ho  = oxidizer->enthalpy_mass();
   double *Yo = new double[nsp];
   double *Xo = new double[nsp];
   oxidizer->getMassFractions(Yo);
   oxidizer->getMoleFractions(Xo);

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
   stoichio->setState_HP(Hst,listFlames[0]->m_Pressure);

   // get fuel and oxidizer mass fraction in stoichiometric mixture
   double Yf_st,Yo_st;
   getYfYo(stoichio,Yf_st,Yo_st);
   
   //printf("Yf_st: %12.5f   Yo_st: %12.5f \n", Yf_st, Yo_st);

   // compute useful quantities
   double s = Yo_st/Yf_st;
   double Z1_O = s*Yf_O - Yo_O;
   double Z1_F = s*Yf_F - Yo_F;
   
   double z = (Yo_O*(listFlames[0]->m_Phi-1.0) - Z1_O)/((Z1_F-Z1_O)-(listFlames[0]->m_Phi-1.0)*(Yo_F-Yo_O));

   // Mixture state
   for (int ns=0;ns<nsp;ns++)
      Ym[ns] = z*Yf[ns]+(1.0-z)*Yo[ns];
   Hm = z*Hf+ (1.0-z)*Ho;

   //cout << "Mixture fraction stoichio " << zst << endl;
   //cout << "Mixture fraction " << z << endl;
   //cout << "Z1_0 " << Z1_O << endl;
   //cout << "Z1_F " << Z1_F << endl;
   //cout << "Yf_F " << Yf_F << endl;
   //cout << "Yf_O " << Yf_O << endl;
   //cout << "Yo_F " << Yo_F << endl;
   //cout << "Yo_O " << Yo_O << endl;
   //cout << "s " << s << endl;
   cout << "               Pressure: " << listFlames[0]->m_Pressure << endl;
   cout << "               Equivalence ratio: " << listFlames[0]->m_Phi << endl;
   cout << "               Composition (mass fractions): " << endl;  
   for (int k=0; k<nsp; k++)
   {
      if (Ym[k] != 0.0)
         cout << "                                    <" << mixture->speciesName(k) << ":" << Ym[k] << ">" << endl;
   }

   mixture->setMassFractions(Ym);
   mixture->setState_HP(Hm, listFlames[0]->m_Pressure);

   for (int k=0; k<nsp; k++)
   {
      if (Ym[k] > 0.0)
      {
         SpeciesIntoReactants[k] = true;
      }
   }

   time_t start_adjusting, end_adjusting;
   time_t start_converging, end_converging;

   try
   {

      time (&start_adjusting);
      cout<<"               Adjusting flame..."<<endl;
      cout<<"               ";
      {
        double ratio=2;
        double slope=0.05;
        double curve=0.8;
        int loglevel = 1;
        double coeff_start=1.0e-3;
        flame->setRefineCriteria(flowdomain,ratio,slope,curve,-1);
        flame->set_max_isteps(30);
        flame->setJacAge(30);
        flame->newton().set_jacMax(30);
        flame->newton().setOptions(30);
        double mesh_crit=1.0e-3;
        goto_state(flame,coeff_start,loglevel);
        remesh_and_converge(flame,mesh_crit,loglevel);
      }

      time (&start_converging);
   
      cout<<endl;
      cout<<"               Converging flame..."<<endl;
      cout<<"               ";
      {
        double ratio=2;
        double slope=0.05;
        double curve=0.8;
        int loglevel = 1;
        double mesh_crit=1.0e-3;
        flame->setRefineCriteria(flowdomain,ratio,slope,curve,-1);
        flame->set_max_isteps(30);
        flame->setJacAge(20);
        flame->newton().set_jacMax(20);
        flame->newton().setOptions(20);
        remesh_and_converge(flame,mesh_crit,loglevel);
      }
      time (&end_converging);
    //  int diff_time_converging = end_converging-start_converging;
    //  if (diff_time_converging/60 > 0)
    //     cout << "-> " << int(diff_time_converging/60) << "min and " << diff_time_converging%60 << "s" << endl;
    //  else
    //     cout << diff_time_converging << "s" << endl;
   }
   catch(...)
   {
      cout << "An exception occurred. Flame not converged" << endl;
      throw -1;
   }
   cout  << "Success !"  << endl  << endl;
   flame->save(listFlames[0]->m_Final_Flame, "flame", "Saved flame");

   int npoints = flow->nPoints();

   double position_max_wdot = 0;
   double max_wdot = 0;
   vector<double> max_wdot_k (nsp, 0.0); 

//string data = "wdot.dat";
//ofstream dat(data.c_str());
//dat << "#1:wdot_CO" << endl;
  
   for (int ino=0; ino<npoints; ino++)
   {
      for (int k=0; k<nsp; k++)
      {
         if (fabs(max_wdot_k[k]) < fabs(flow->masswdot(k,ino)))
         {
            max_wdot_k[k] = flow->masswdot(k,ino);
         }

         if (max_wdot < fabs(flow->masswdot(k,ino)))
         {
            max_wdot = flow->masswdot(k,ino);
            position_max_wdot = flow->grid(ino);
         }
       


      }
   }


   int ino_start = 0;
   int ino_stop = 0;

   for (int ino=0; ino<npoints; ino++)
   {
      for (int k=0; k<nsp; k++)
      {
         if (ino_start == 0 && fabs(flow->masswdot(k,ino))>0.05*fabs(max_wdot_k[k]))
         {
            ino_start = ino;
         }
      }
   }

   for (int ino=npoints-1; ino>0; ino--)
   {
      for (int k=0; k<nsp; k++)
      {
         if (ino_stop == 0 && fabs(flow->masswdot(k,ino)>0.05*fabs(max_wdot_k[k])))
         {
            ino_stop = ino;
         }
      }
   }

   //cout << "starting node " << ino_start << "   stopping node " << ino_stop << " to be within the region with at least 5% of species wdot" << endl;
   
   vector<vector<double> > position (listFlames.size(), vector<double> (npoints, 0.0));
   vector<vector<vector<double> > > Y (listFlames.size(), vector<vector<double> > (npoints, vector<double> (nsp, 0.0)));
   vector<vector<double> > T (listFlames.size(), vector<double> (npoints, 0.0));
   vector<vector<double> > U (listFlames.size(), vector<double> (npoints, 0.0));

   for (int ino=0; ino<npoints; ino++)
   {
      position[0][ino] = flow->grid(ino);
      T[0][ino] = flame->value(flowdomain, flow->componentIndex("T"), ino);
      U[0][ino] = flame->value(flowdomain, flow->componentIndex("u"), ino);
      for (int k=0; k<nsp; k++)
      {
         string sp_name = mixture->speciesName(k);
         int sp_index = flow->componentIndex(sp_name);
         Y[0][ino][k] = flame->value(flowdomain, sp_index, ino);
      }
   }

   output_datas (outputName, Y[0], T[0], U[0], position[0], position_max_wdot, mixture);


 if (step=="ComputeTrajectories")
  {

    vector<vector<vector<double> > > CreationRate (listFlames.size(), vector<vector<double> > (npoints, vector<double> (nsp, 0.0)));
    vector<vector<vector<double> > > DestructionRate (listFlames.size(), vector<vector<double> > (npoints, vector<double> (nsp, 0.0)));
    vector<vector<vector<double> > > w_Y (listFlames.size(), vector<vector<double> > (npoints, vector<double> (nsp, 0.0)));

/*    ofstream sourceTerm("sourceTerm.dat");

    for (int ino=0; ino<npoints; ino++)
    {
       for (int k=0; k<nsp; k++)
       {
          CreationRate[0][ino][k] = flow->creationRate(k, ino);
          DestructionRate[0][ino][k] = flow->destructionRate(k, ino);
          w_Y[0][ino][k] = CreationRate[0][ino][k]-DestructionRate[0][ino][k];
          sourceTerm << w_Y[0][ino][k] << "  ";


       }


        sourceTerm << endl;

    }
    sourceTerm.close();
*/

  }





   if (step == "computeQSSCriteria")
   {
      for (int k=0; k<nsp; k++)
         QSS_Criteria[k][1] = k;

      vector<vector<vector<double> > > CreationRate (listFlames.size(), vector<vector<double> > (npoints, vector<double> (nsp, 0.0)));
      vector<vector<vector<double> > > DestructionRate (listFlames.size(), vector<vector<double> > (npoints, vector<double> (nsp, 0.0)));

      ofstream log("computeQSSCriteria.log");
      log << endl;

      for (int ino=0; ino<npoints; ino++)
      {
         for (int k=0; k<nsp; k++)
         {
            CreationRate[0][ino][k] = flow->creationRate(k, ino);
            DestructionRate[0][ino][k] = flow->destructionRate(k, ino);
         }
      }

      vector<double> inte_Net (nsp, 0.0);
      vector<double> inte_Creation (nsp, 0.0);
      vector<double> inte_Destruction (nsp, 0.0);

      for (int ino=ino_start; ino<ino_stop; ino++)
      {
         for (int k=0; k<nsp; k++)
         {
            double delta_position = position[0][ino]-position[0][ino-1];
            inte_Net[k] += abs(CreationRate[0][ino][k]-DestructionRate[0][ino][k])*delta_position;
            inte_Creation[k] += CreationRate[0][ino][k]*delta_position;
            inte_Destruction[k] += DestructionRate[0][ino][k]*delta_position;
         }
      }

      for (int k=0; k<nsp; k++)
      {
         if (QSS_Criteria[k][0] < inte_Net[k]/max(inte_Creation[k], inte_Destruction[k]))
            QSS_Criteria[k][0] = inte_Net[k]/max(inte_Creation[k], inte_Destruction[k]);
      }


      log << "mech = " << mech ;

      cout << endl  << endl << endl << "--------QSS coefficients--------" << endl;
      cout << "--------------------------------" << endl;
      log << endl  << endl << "--------QSS coefficients--------" << endl;
      log << "--------------------------------" << endl;


      for (int k=0; k<nsp; k++)
      {
         cout << "Species " << mixture->speciesName(k) << "  " << QSS_Criteria[k][0] << endl;
         log << "Species " << mixture->speciesName(k) << "  " << QSS_Criteria[k][0] << endl;
      }
      cout << "--------------------------------" << endl << endl;

   log.close();

   }




   if (Proceed_DRGEP_analysis)
   {
      if (step == "DRGEP_Species" || step == "DRGEP_Reactions")
      {
         drgep *species_relations = new drgep();

         //----------Limits should be given to compute DRGEP in the inner region of the flame (not burned gases region)----------//
         //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
         for (int ino=ino_start; ino<ino_stop; ino++)
         {
            species_relations->drgep_1D_species(mixture, flow, ino, Targets, R_AD_Premixed, position_max_wdot);
         }

         if (step == "DRGEP_Reactions")
         {
            vector<vector<double> > rj_for_k (nsp, vector<double> (nreac,0.0));
            vector<vector<double> > rjf_for_k (nsp, vector<double> (nreac,0.0));
            vector<vector<double> > rjr_for_k (nsp, vector<double> (nreac,0.0));
            for (int ino=ino_start; ino<ino_stop; ino++)
            {
               species_relations->drgep_1D_reactions(mixture, flow, ino, rj_for_k, rjf_for_k, rjr_for_k);

               //for (int j=0; j<nreac; j++)
               //   cout << "rj_for_k[kb][j] " << rj_for_k[0][j] << endl;
               //getchar();


               for (int ka=0; ka<nsp; ka++)
               {
                  for (int kb=0; kb<nsp; kb++)
                  {
                     for (int j=0; j<nreac; j++)
                     {

                        if (max_j_on_Target[ka][j] < R_AD_Premixed[ka][kb]*rj_for_k[kb][j])
                           max_j_on_Target[ka][j] = R_AD_Premixed[ka][kb]*rj_for_k[kb][j];

                        if (max_jf_on_Target[ka][j] < R_AD_Premixed[ka][kb]*rjf_for_k[kb][j])
                           max_jf_on_Target[ka][j] = R_AD_Premixed[ka][kb]*rjf_for_k[kb][j];

                        if (max_jr_on_Target[ka][j] < R_AD_Premixed[ka][kb]*rjr_for_k[kb][j])
                           max_jr_on_Target[ka][j] = R_AD_Premixed[ka][kb]*rjr_for_k[kb][j];

                     }
                  }
               }
            }
         }
      }
   }







}

