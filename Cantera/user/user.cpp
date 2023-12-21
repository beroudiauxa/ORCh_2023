// GL => OK

/*
 * This is an example of a user-defined function that can be linked into Cantera.
 */
#include "user.h"

using namespace Cantera;

namespace User {
  using namespace std;
  using namespace Cantera_CXX;


  double flame_distance(double T1, double T2, double P1, double P2, double* Y1, double* Y2, int nsp)
  {
    double distance = 0.0;
    distance = fmaxx(distance,abs(log10(T1/T2)));
    distance = fmaxx(distance,abs(log10(P1/P2)));
    for (int i=0;i<nsp;i++)
    {
      distance = fmaxx(distance,abs(abs(Y1[i])-abs(Y2[i])));
    }
    return distance;
  }

  void findStoechiometryFromCompo(IdealGasMix* Fuel, IdealGasMix* Oxidizer, double* zst)
  {
    // ==========================================================================
    // Number of species
    int nsp = Fuel->nSpecies();
    int nel = Fuel->nElements();
    ChemEquil c;
    // ==========================================================================
    // Molar mass of C atom
    int index_C = Fuel->elementIndex("C");
    double W_C;
    if (index_C<0)
    {
      W_C = 1.0e12;
    }
    else
    {
      W_C = Fuel->atomicWeight(Fuel->elementIndex("C"));
    }

    // ==========================================================================
    // Molar mass of H atom
    int index_H = Fuel->elementIndex("H");
    double W_H;
    if (index_H<0)
    {
      W_H = 1.0e12;
    }
    else
    {
      W_H = Fuel->atomicWeight(Fuel->elementIndex("H"));
    }

    // ==========================================================================
    // Molar mass of O atom
    int index_O = Fuel->elementIndex("O");
    double W_O;
    if (index_O<0)
    {
      W_O = 1.0e12;
    }
    else
    {
      W_O = Fuel->atomicWeight(Fuel->elementIndex("O"));
    }

    // ==========================================================================
    // Find mole and mass fractions of C, H and O in Fuel ( Y_C^f, Y_H^f and Y_O^f)
    double Yf_C = 0.0;
    double Yf_H = 0.0;
    double Yf_O = 0.0;
    for (int i=0; i<nsp; i++)
    {
      double Wsp = Fuel->molecularWeight(i);
      string sp_name = Fuel->speciesName(i);
      double Ysp = Fuel->massFraction(sp_name);
      for (int j=0; j<nel; j++)
      {
        int nat = Fuel->nAtoms(i,j);
        double Wel = Fuel->atomicWeight(j);
        string el_name = Fuel->elementName(j);
        double Yat = nat*Wel/Wsp;
        if (el_name=="C")
        {
          Yf_C += Ysp * Yat;
        }
        if (el_name=="H")
        {
          Yf_H += Ysp * Yat;
        }
        if (el_name=="O")
        {
          Yf_O += Ysp * Yat;
        }
      }
    }

    //cout << "Mass fraction of C in fuel  " << Yf_C << endl;
    //cout << "Mass fraction of H in fuel  " << Yf_H << endl;
    //cout << "Mass fraction of O in fuel  " << Yf_O << endl;
    //cout << "Sum mass fractions " << Yf_C+Yf_H+Yf_O << endl;

    // ==========================================================================
    // Find mole and mass fractions of C, H and O in Oxidizer ( Y_C^o, Y_H^o and Y_O^o)
    double Yo_C = 0.0;
    double Yo_H = 0.0;
    double Yo_O = 0.0;
    for (int i=0; i<nsp; i++)
    {
      double Wsp = Oxidizer->molecularWeight(i);
      string sp_name = Oxidizer->speciesName(i);
      double Ysp = Oxidizer->massFraction(sp_name);
      for (int j=0; j<nel; j++)
      {
        int nat = Oxidizer->nAtoms(i,j);
        double Wel = Oxidizer->atomicWeight(j);
        string el_name = Oxidizer->elementName(j);
        double Yat = nat*Wel/Wsp;
        if (el_name=="C")
        {
          Yo_C += Ysp * Yat;
        }
        if (el_name=="H")
        {
          Yo_H += Ysp * Yat;
        }
        if (el_name=="O")
        {
          Yo_O += Ysp * Yat;
        }
      }
    }


    // ==========================================================================
    // X_O = 2.0*X_C + 0.5*X_H for complete stoechiometric combustion
    // Y_O/W_O = 2.0*Y_C/W_C + 0.5*Y_H/W_H for complete stoechiometric combustion
    // with Y_C = z*Y_C^f + (1-z)*Y_C^o, ...

    double num = 0.0;
    num += -2.00*Yo_C/W_C;
    num += -0.50*Yo_H/W_H;
    num += +1.00*Yo_O/W_O;

    double den = 0.0;
    den += +2.00*(Yf_C-Yo_C)/W_C;
    den += +0.50*(Yf_H-Yo_H)/W_H;
    den += -1.00*(Yf_O-Yo_O)/W_O;

    *zst = num/den;
  }

  void findStoechiometryFromTmax(IdealGasMix* Fuel, IdealGasMix* Oxidizer, double* zst)
  {
    // ==========================================================================
    // Number of species
    int nspF = Fuel->nSpecies();
    int nspO = Oxidizer->nSpecies();
    if (nspF!=nspO)
    {
      cout<<"Number of species is different in Fuel and Oxidizer!!!"<<endl;
      exit(-1);
    }
    int nsp = nspF; // = nspO

    // ==========================================================================
    // Fuel
    double Tf = Fuel->temperature();
    double Pf = Fuel->pressure();
    double Hf = Fuel->enthalpy_mass();
    double *Yf=new double[nsp];
    Fuel->getMassFractions(Yf);

    // ==========================================================================
    // Oxidizer
    double To = Oxidizer->temperature();
    double Po = Oxidizer->pressure();
    double Ho = Oxidizer->enthalpy_mass();
    double *Yo=new double[nsp];
    Oxidizer->getMassFractions(Yo);

    if (abs(Pf-Po)/fmaxx(Pf,Po)>1.0e-12)
    {
      cout<<"Pressure is different in Fuel and Oxidizer!!!"<<endl;
      exit(-1);
    }
    double P = Pf; // = Po

    // ==========================================================================
    // Allocate data
    IdealGasMix* mixture  = new IdealGasMix(*Fuel);
    double *Ym = new double[nsp];

    // ==========================================================================
    // Main loop
    int nz = 10;
    bool done = false;
    double tol = 1.0e-4;
    double zl = 0.0;
    double zr = 1.0;
    double Tmax = fmaxx(Tf,To);
    ChemEquil c;
    while(!done)
    {
      int imax=-1;
      double dz = (zr-zl)/nz;
      for (int i=0; i<nz+1; i++)
      {
        double z=zl + i *dz;
        for (int ns=0;ns<nsp;ns++)
        {
          Ym[ns] = z*Yf[ns]+(1.0-z)*Yo[ns];
        }
        mixture->setMassFractions(Ym);
        double Hm = z*Hf+ (1.0-z)*Ho;
        mixture->setState_HP(Hm,P);
        bool fail_eq = false;
        try
        {
         c.equilibrate(*mixture,"HP");
        }
        catch (...)
        {
          fail_eq = true;
        }
        if (!fail_eq)
        {
          double Tm = mixture->temperature();
          if (Tm>=Tmax)
          {
            Tmax = Tm;
            imax = i;
          }
        }
      }
      if ((imax==-1)||(imax==0)||(imax==nz))
      {
        cout<<"No maximum can be found"<<endl;
        exit(-1);
      }
      double zl_new = zl + (imax-1)*dz;
      double zr_new = zl + (imax+1)*dz;
      zl = zl_new;
      zr = zr_new;
      if ((zr-zl)<tol)
      {
        done = true;
        *zst = 0.5*(zl+zr);
      }
    }

  }

  void getYfYo(IdealGasMix* mixture, double& Yf, double& Yo)
  {
    // ==========================================================================
    // Number of species
    int nsp = mixture->nSpecies();
    int nel = mixture->nElements();

    // ==========================================================================
    // Find mass fractions of C, H and O in mixture
    double Y_C = 0.0;
    double Y_H = 0.0;
    double Y_O = 0.0;
    for (int i=0; i<nsp; i++)
    {
      double Wsp = mixture->molecularWeight(i);
      string sp_name = mixture->speciesName(i);
      double Ysp = mixture->massFraction(sp_name);
      for (int j=0; j<nel; j++)
      {
        int nat = mixture->nAtoms(i,j);
        double Wel = mixture->atomicWeight(j);
        string el_name = mixture->elementName(j);
        double Yat = nat*Wel/Wsp;
        if (el_name=="C")
        {
          Y_C += Ysp * Yat;
        }
        if (el_name=="H")
        {
          Y_H += Ysp * Yat;
        }
        if (el_name=="O")
        {
          Y_O += Ysp * Yat;
        }
      }
    }

    // ==========================================================================
    // Fuel is C+H, Oxisizer is O
    Yf = Y_C + Y_H;
    Yo = Y_O;

  }

  void remesh_and_converge(Sim1D* flame, double criterion, int loglevel)
  {
    // Do not refine grid during solve
    bool refine_grid=false;

    // Pointers to various objects
    OneDim* OD = flame;
    StFlow* flow=NULL;
    Inlet1D* inlet=NULL;
    Outlet1D* outlet=NULL;
    int flowdomain=-1;

    // Number of domains in simulation
    int nd = OD->nDomains();


    // Get back inlet, flow and outlet
    for (int n = 0; n < nd; n++)
    {
      if (OD->domain(n).domainType() == cFreeFlow)
      {
        flow = static_cast<StFlow*>(&OD->domain(n));
        flowdomain = n;
      }
      if (OD->domain(n).domainType() == cInletType)
      {
        inlet = static_cast<Inlet1D*>(&OD->domain(n));
      }
      if (OD->domain(n).domainType() == cOutletType)
      {
        outlet = static_cast<Outlet1D*>(&OD->domain(n));
      }
    }

    // Get mixture
    IdealGasMix* mixture = static_cast<IdealGasMix*>(&flow->phase());

    // print info
    if (loglevel>1)
    cout<<"INITIAL MESH = ["<<setw(4)<<flow->nPoints()<<"]"<<endl;

    // main loop
    bool done=false;
    while (!done)
    {

      // number of species
      int nsp = mixture->nSpecies();

      // last values
      double *Yin=new double[nsp];
      double *Xin=new double[nsp];
      double Tin;
      double P;

      // Get flame values
      Tin = flame->value(flowdomain,flow->componentIndex("T"),0);
      P = flow->pressure();
      for (int i=0;i<nsp;i++)
      {
        string specname = mixture->speciesName(i);
        Yin[i] = flame->value(flowdomain,flow->componentIndex(specname),0);
      }

      // Set flame properties
      flow->setPressure(P);
      mixture->setState_TPY(Tin, P, Yin);
      mixture->getMoleFractions(Xin);
      inlet->setMoleFractions(Xin);
      inlet->setMdot(0.0);
      inlet->setTemperature(Tin);

      // Old number of points in the flame
      int np_old=flow->nPoints();

      // Remesh
      double mesh_distance = flame->remesh(loglevel);


      // New number of point in the flame
      int np_new=flow->nPoints();

      // Set fixed temperature
      double T_f = flame->value(flowdomain,flow->componentIndex("T"),np_new/2);
      flame->setFixedTemperature(T_f);

      // Solve
      //cout << "." << flush;
      if (loglevel>1)
         cout<<"old mesh=["<<setw(4)<<np_old<<"]"<<" new mesh = ["<<setw(4)<<np_new<<"]   distance = "<<mesh_distance<<flush;
      //OD->jacobian().resetnEvals();
      flame->solve(loglevel,refine_grid);

      // we are done!
      if (mesh_distance < criterion)
      {
        done=true;
      }
      // another step is needed
      else
      {
        done=false;
      }
    }
    if (loglevel>1)
       cout<<"Success"<<endl;
}



void goto_state(Sim1D* flame, double alpha, int loglevel)
{
    // Do not refine grid during solve
    bool refine_grid=false;

    // Pointers to various objects
    OneDim* OD = flame;
    //IdealGasMix* mixture=NULL;
    StFlow* flow=NULL;
    Inlet1D* inlet=NULL;
    Outlet1D* outlet=NULL;
    int flowdomain=-1;

    

    // Number of domains in simulation
    int nd = OD->nDomains();
    // Get back inlet, flow and outlet
    for (int n = 0; n < nd; n++)
    {

      if (OD->domain(n).domainType() == cFreeFlow)
      {
        flow = static_cast<StFlow*>(&OD->domain(n));
        flowdomain = n;
      }
      if (OD->domain(n).domainType() == cInletType)
      {
        inlet = static_cast<Inlet1D*>(&OD->domain(n));
      }
      if (OD->domain(n).domainType() == cOutletType)
      {
        outlet = static_cast<Outlet1D*>(&OD->domain(n));
      }
    }

    // Get mixture
    IdealGasMix* mixture = static_cast<IdealGasMix*>(&flow->phase());


    // Number of species in the mixture
    int nsp = mixture->nSpecies();
    // At start, last values are those of the reference flame
    double Tin_last = flame->value(flowdomain,flow->componentIndex("T"),0);
    double P_last = flow->pressure();
    double *Yin_last=new double[nsp];
    for (int i=0;i<nsp;i++)
    {
      string specname = mixture->speciesName(i);
      Yin_last[i] = flame->value(flowdomain,flow->componentIndex(specname),0);
    }

    // Get back final (wanted)  values from flame->mixture
    double Tin_final = mixture->temperature();
    double P_final = mixture->pressure();
    double *Yin_final = new double[nsp];
    mixture->getMassFractions(Yin_final);



    // Define current flame values (between last and final state)
    double Tin_curr = Tin_last + alpha*(Tin_final-Tin_last);
    double P_curr = P_last + alpha*(P_final-P_last);
    double *Yin_curr=new double[nsp];
    double *Xin_curr=new double[nsp];
    for (int i=0;i<nsp;i++)
    {
      string specname = mixture->speciesName(i);
      Yin_curr[i] = Yin_last[i]+ alpha*(Yin_final[i] - Yin_last[i]);
    }




    double alpha_last = alpha;
    double alpha_curr = alpha;


    // Main loop
    bool done=false;
    bool success=true;
    while (!done)
    {
      // Print distance info
      double dist0 = flame_distance(Tin_final,Tin_last,P_final,P_last,Yin_final,Yin_last,nsp);
      double dist1 = flame_distance(Tin_final,Tin_curr,P_final,P_curr,Yin_final,Yin_curr,nsp);
      double step_size  = flame_distance(Tin_last,Tin_curr,P_last,P_curr,Yin_last,Yin_curr,nsp);

      if (loglevel>1)
         cout<<"dist0 = "<<dist0<<" step = "<<step_size<<" dist1 = "<<dist1<<flush<<" alpha = "<<alpha_curr<<" |"<<flush;

      cout << "." << flush;
/*
    //Compute Phi for BFER
    double phi; 
    double phi_final; 
    double Y_CH4;
    double Y_O2;
    double Y_CH4_final;
    double Y_O2_final;

    for (int i=0; i<nsp; i++)
    {
      string specname = mixture->speciesName(i);
      if (specname == "CH4")
      {
        Y_CH4 = Yin_curr[i];
        Y_CH4_final = Yin_final[i];
      }
      if (specname == "O2")
      {
        Y_O2 = Yin_curr[i];
        Y_O2_final = Yin_final[i];
      }
    }

    phi = 3.98918889*Y_CH4/Y_O2;
    phi_final = 3.98918889*Y_CH4_final/Y_O2_final;
    cout << " phi " << phi << endl;
    cout << " phi_final " << phi_final << endl;
    m_phi = phi;
*/




      // Set flame properties
      flow->setPressure(P_curr);
      mixture->setState_TPY(Tin_curr, P_curr, Yin_curr);
      mixture->getMoleFractions(Xin_curr);
      inlet->setMoleFractions(Xin_curr);
      inlet->setMdot(0.0);
      inlet->setTemperature(Tin_curr);
      // Try to solve the problem
      int njac=0;
      int neval=0;
      try
      {
        if (success)
        {
          // Remesh
          if (loglevel>1)
          cout<<" remeshing"<<flush;
          double mesh_distance = flame->remesh(loglevel);
          if (loglevel>1)
          cout<<" (dist= "<<mesh_distance<<") |"<<flush;
        }
        // Set fixed temperature
        int np=flow->nPoints();
        double T_f = flame->value(flowdomain,flow->componentIndex("T"),np/2);
        flame->setFixedTemperature(T_f);
        // Reset the number of evaluation of the jacobian
        OD->jacobian().resetnEvals();
        // Solve
        if (loglevel>1)
           cout<<" solving mesh=["<<setw(4)<<np<<"]..."<<flush;
        flame->solve(loglevel,refine_grid);
        // Get back number of evaluation of the jacobian
        njac  = OD->jacobian().nEvals();
        //neval = OD->newton().getOptions();
      }
      // Check success or failure
      catch (CanteraError&)
      {
        cout << "cantera error" << endl;
        success=false;
      }
      if (flame->newton_failure())
      {
        cout << "newton_failed" << endl;
        success=false;
      }
      else
      {
        success=true;
      }


      // If success
      if (success)
      {

        // Save flame in case of crash
        //if (system("rm -f st_flame_tmp.cantera")==0)
        //{
        //  flame->save("st_flame_tmp.cantera", "st_flame", "Stoichiometric flame");
        //}

        // We are done!
        if (dist1 == 0.0)
        {
          done=true;
        }
        // Another step is needed
        else
        {
          done=false;

          // Adapt the next step
          double coeff_next;
          int niter= njac*neval;
          if (niter<30)
          {
            coeff_next=1.50;
          }
          else if (niter<60)
          {
            coeff_next=1.25;
          }
          else if (niter<90)
          {
            coeff_next=1.10;
          }
          else if (niter<110)
          {
            coeff_next=1.00;
          }
          else if (niter<160)
          {
            coeff_next=0.90;
          }
          else if (niter<210)
          {
            coeff_next=0.75;
          }
          else if (niter<310)
          {
            coeff_next=0.50;
          }
          else
          {
            coeff_next=0.25;
          }
          alpha_last = alpha_curr;
          alpha_curr = fminn(1.0,coeff_next*alpha_last);
          if ((alpha_last>0.5)&&(coeff_next>1.01))
          {
            alpha_curr=1.0;
          }

          // Define next flame values
          double delTin = alpha_curr*(Tin_final - Tin_curr);
          Tin_last = Tin_curr;
          Tin_curr = Tin_last + delTin;

          double delP   = alpha_curr*(P_final - P_curr);
          P_last = P_curr;
          P_curr = P_last + delP;

          for (int i=0;i<nsp;i++)
          {
            string specname = mixture->speciesName(i);
            double delY   = alpha_curr*(Yin_final[i] - Yin_curr[i]);
            Yin_last[i] = Yin_curr[i];
            Yin_curr[i] = Yin_last[i] + delY;
          }
        }
      }
      // if failure, change step size and restart
      else
      {  
        if (loglevel>1)
           cout<<" failure!!!!"<<endl;
        done=false;

        // Adapt the next step
        double coeff_next;
        if (alpha_curr>=alpha_last)
        {
          coeff_next = 0.5;
        }
        else
        {
          coeff_next = 0.2;
        }
        alpha_last = alpha_curr;
        alpha_curr = coeff_next*alpha_last;

        // The next step is too small: the flame is near extinction
        if (alpha_curr<1.0e-6)
        {
          throw -1;
        }

        // Define next flame values
        double delTin = alpha_curr*(Tin_final - Tin_last);
        Tin_curr = Tin_last + delTin;

        double delP   = alpha_curr*(P_final - P_last);
        P_curr = P_last + delP;

        for (int i=0;i<nsp;i++)
        {
          string specname = mixture->speciesName(i);
          double delY   = alpha_curr*(Yin_final[i] - Yin_last[i]);
          Yin_curr[i] = Yin_last[i] + delY;
        }
      }
    }

    if (loglevel>0)
    cout << endl <<"               Success !" << endl;
  }


  void createFlameFromScratch(Sim1D* flame, int nint0, doublereal lz, int loglevel)
  {

    // ==========================================================================
    // Various items
    bool do_not_refine_grid=false;
    double atol;
    double rtol;
    double ratio;
    double slope;
    double curve;
    double mesh_crit;
    double Uin=1.0; // in [m/sec]

    // ==========================================================================
    // Pointers to various objects
    OneDim* OD = flame;
    StFlow* flow=NULL;
    Inlet1D* inlet=NULL;
    Outlet1D* outlet=NULL;
    ChemEquil c;
    int flowdomain=-1;

    // ==========================================================================
    // Number of domains in simulation
    int nd = OD->nDomains();

    // ==========================================================================
    // Get back inlet, flow and outlet
    for (int n = 0; n < nd; n++)
    {
      if (OD->domain(n).domainType() == cFreeFlow)
      {
        flow = static_cast<StFlow*>(&OD->domain(n));
        flowdomain = n;
      }
      if (OD->domain(n).domainType() == cInletType)
      {
        inlet = static_cast<Inlet1D*>(&OD->domain(n));
      }
      if (OD->domain(n).domainType() == cOutletType)
      {
        outlet = static_cast<Outlet1D*>(&OD->domain(n));
      }
    }

    // ==========================================================================
    // Get mixture
    IdealGasMix* mixture = static_cast<IdealGasMix*>(&flow->phase());
    int nsp = mixture->nSpecies();
    double Tm = mixture->temperature();
    double Pm = mixture->pressure();
    double *Ym=new double[nsp];
    mixture->getMassFractions(Ym);

    // ==========================================================================
    // Mixture state (inlet)
    IdealGasMix* mixture_in  = new IdealGasMix(*mixture);
    mixture_in->setState_TPY(Tm,Pm,Ym);
    double *Yin=new double[nsp];
    double *Xin=new double[nsp];
    mixture_in->getMassFractions(Yin);
    mixture_in->getMoleFractions(Xin);
    double rho_in=mixture_in->density();
    double Tin = mixture_in->temperature();

    // ==========================================================================
    // Mixture state (outlet)
    IdealGasMix* mixture_out = new IdealGasMix(*mixture);
    mixture_out->setState_TPY(Tm,Pm,Ym);
    c.equilibrate(*mixture_out,"HP");
    double *Yout=new double[nsp];
    mixture_out->getMassFractions(Yout);
    double rho_out = mixture_out->density();
    double Tout = mixture_out->temperature();

    // ==========================================================================
    // Set the inlet
    double mdot=Uin*rho_in;
    inlet->setMdot(mdot);
    inlet->setMoleFractions(Xin);
    inlet->setTemperature(Tin);

    // ==========================================================================
    // Build the domain with an initial grid
    int np0 = nint0+3;
    double buff=10.0*lz;
    vector_fp zg(np0);
    doublereal dz=lz/((doublereal)(nint0));
    zg[0]=0.0;
    for(int iz=1;iz<np0-1;iz++)
    {
        zg[iz]=buff+((doublereal)(iz-1))*dz;
    }
    zg[np0-1]=buff+lz+buff;
    flow->setupGrid(np0,&zg[0]); //MAC
    flow->setPressure(Pm);
    // ==========================================================================
    // Supply an initial guess: ramp values from inlet to adiabatic flame conditions 
    vector_fp locs;
    vector_fp value;
    locs.resize(4);
    value.resize(4);

    double z1=0;
    double z2=buff/(buff+lz+buff);
    double z3=(buff+lz)/(buff+lz+buff);
    double z4=1.0;
    locs[0]=z1; locs[1]=z2; locs[2]=z3; locs[3]=z4;

    double Uout=inlet->mdot()/rho_out;
    value[0]=Uin; value[1]=Uin; value[2]=Uout; value[3]=Uout;
    flame->setInitialGuess("u",locs,value);

    value[0]=Tin; value[1]=Tin; value[2]=Tout; value[3]=Tout;
    flame->setInitialGuess("T",locs,value);

    for(int i=0;i<nsp;i++)
    {
        value[0]=Yin[i]; value[1]=Yin[i]; value[2]=Yout[i]; value[3]=Yout[i];
        flame->setInitialGuess(mixture->speciesName(i),locs,value);
    }

    cout<<"================================================================"<<endl;
    cout<<"                Fixed Temperature Problem                       "<<endl;
    cout<<"================================================================"<<endl;
//
    flow->fixTemperature();
    flame->setFixedTemperature(flame->value(flowdomain,flow->componentIndex("T"),flow->nPoints()/2));
//
    rtol=1.0e-3; atol=1.0e-6;
    cout<<"RTOL = "<<rtol<<" ATOL="<<atol<<endl;
   // flow->setTolerances(rtol,atol,0);
    flame->solve(loglevel,do_not_refine_grid);
//
    rtol=1.0e-4; atol=1.0e-8;
    cout<<"RTOL = "<<rtol<<" ATOL="<<atol<<endl;
   // flow->setTolerances(rtol,atol,0);
    flame->solve(loglevel,do_not_refine_grid);
//
    rtol=1.0e-5; atol=1.0e-10;
    cout<<"RTOL = "<<rtol<<" ATOL="<<atol<<endl;
   // flow->setTolerances(rtol,atol,0);
    flame->solve(loglevel,do_not_refine_grid);
//
    rtol=1.0e-6; atol=1.0e-12;
    cout<<"RTOL = "<<rtol<<" ATOL="<<atol<<endl;
   // flow->setTolerances(rtol,atol,0);
    flame->solve(loglevel,do_not_refine_grid);
//
    rtol=1.0e-9; atol=1.0e-15;
    cout<<"RTOL = "<<rtol<<" ATOL="<<atol<<endl;
   // flow->setTolerances(rtol,atol,0);
    flame->solve(loglevel,do_not_refine_grid);
//

    cout<<"================================================================"<<endl;
    cout<<"                Energy Equation Problem                         "<<endl;
    cout<<"================================================================"<<endl;
//
    flow->solveEnergyEqn();
//
    rtol=1.0e-0; atol=1.0e-0; ratio=50.0; slope=0.9; curve=0.9;
    cout<<"RTOL = "<<rtol<<" ATOL="<<atol<<" slope="<<slope<<" curve="<<curve<<endl;
    flame->setRefineCriteria(flowdomain,ratio,slope,curve,-1);
  //  flow->setTolerances(rtol,atol,0);
    flame->solve(loglevel,do_not_refine_grid);
//
    rtol=1.0e-1; atol=1.0e-2; ratio=50.0; slope=0.8; curve=0.8;
    cout<<"RTOL = "<<rtol<<" ATOL="<<atol<<" slope="<<slope<<" curve="<<curve<<endl;
    flame->setRefineCriteria(flowdomain,ratio,slope,curve,-1);
  //  flow->setTolerances(rtol,atol,0);
    flame->solve(loglevel,do_not_refine_grid);
//
    rtol=1.0e-2; atol=1.0e-4; ratio=50.0; slope=0.7; curve=0.7;
    cout<<"RTOL = "<<rtol<<" ATOL="<<atol<<" slope="<<slope<<" curve="<<curve<<endl;
    flame->setRefineCriteria(flowdomain,ratio,slope,curve,-1);
  //  flow->setTolerances(rtol,atol,0);
    flame->solve(loglevel,do_not_refine_grid);
//
    rtol=1.0e-3; atol=1.0e-6; ratio=10.0; slope=0.5; curve=0.5;
    cout<<"RTOL = "<<rtol<<" ATOL="<<atol<<" slope="<<slope<<" curve="<<curve<<endl;
    flame->setRefineCriteria(flowdomain,ratio,slope,curve,-1);
  //  flow->setTolerances(rtol,atol,0);
    flame->solve(loglevel,do_not_refine_grid);
//
    rtol=1.0e-4; atol=1.0e-8; ratio=10.0; slope=0.5; curve=0.5; mesh_crit=1.0;
    cout<<"RTOL = "<<rtol<<" ATOL="<<atol<<" slope="<<slope<<" curve="<<curve<<endl;
    flame->setRefineCriteria(flowdomain,ratio,slope,curve,-1);
   // flow->setTolerances(rtol,atol,0);
    remesh_and_converge(flame,mesh_crit,loglevel);
//
    rtol=1.0e-5; atol=1.0e-10; ratio=10.0; slope=0.25; curve=0.25; mesh_crit=1.0e-1;
    cout<<"RTOL = "<<rtol<<" ATOL="<<atol<<" slope="<<slope<<" curve="<<curve<<endl;
    flame->setRefineCriteria(flowdomain,ratio,slope,curve,-1);
   // flow->setTolerances(rtol,atol,0);
    remesh_and_converge(flame,mesh_crit,loglevel);
//
    rtol=1.0e-6; atol=1.0e-12; ratio=5.0; slope=0.10; curve=0.10; mesh_crit=1.0e-2;
    cout<<"RTOL = "<<rtol<<" ATOL="<<atol<<" slope="<<slope<<" curve="<<curve<<endl;
    flame->setRefineCriteria(flowdomain,ratio,slope,curve,-1);
   // flow->setTolerances(rtol,atol,0);
    remesh_and_converge(flame,mesh_crit,loglevel);
//
    rtol=1.0e-9; atol=1.0e-15; ratio=2.0; slope=0.05; curve=0.05; mesh_crit=1.0e-3;
    cout<<"RTOL = "<<rtol<<" ATOL="<<atol<<" slope="<<slope<<" curve="<<curve<<endl;
    flame->setRefineCriteria(flowdomain,ratio,slope,curve,-1);
   // flow->setTolerances(rtol,atol,0);
    remesh_and_converge(flame,mesh_crit,loglevel);
//
    cout<<"DONE!!!"<<endl;
  }

}
