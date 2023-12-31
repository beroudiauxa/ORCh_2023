//! @file ConstPressureReactor.cpp A constant pressure zero-dimensional reactor

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#include "cantera/zeroD/ConstPressureReactor.h"
#include "cantera/zeroD/FlowDevice.h"
#include "cantera/zeroD/Wall.h"
#include "cantera/thermo/SurfPhase.h"

using namespace std;

namespace Cantera
{

void ConstPressureReactor::getState(double* y)
{
    if (m_thermo == 0) {
        throw CanteraError("getState",
                           "Error: reactor is empty.");
    }
    m_thermo->restoreState(m_state);

    // set the first component to the total mass
    y[0] = m_thermo->density() * m_vol;

    // set the second component to the total enthalpy
    y[1] = m_thermo->enthalpy_mass() * m_thermo->density() * m_vol;

    // set components y+2 ... y+K+1 to the mass fractions Y_k of each species
    m_thermo->getMassFractions(y+2);

    // set the remaining components to the surface species
    // coverages on the walls
    getSurfaceInitialConditions(y + m_nsp + 2);
}

void ConstPressureReactor::initialize(doublereal t0)
{
    Reactor::initialize(t0);
    m_nv -= 1; // Constant pressure reactor has one fewer state variable
}

void ConstPressureReactor::updateState(doublereal* y)
{
    // The components of y are [0] the total mass, [1] the total enthalpy,
    // [2...K+2) are the mass fractions of each species, and [K+2...] are the
    // coverages of surface species on each wall.
    m_mass = y[0];
    m_thermo->setMassFractions_NoNorm(y+2);
    if (m_energy) {
        m_thermo->setState_HP(y[1]/m_mass, m_pressure);
    } else {
        m_thermo->setPressure(m_pressure);
    }
    m_vol = m_mass / m_thermo->density();
    updateSurfaceState(y + m_nsp + 2);

    // save parameters needed by other connected reactors
    m_enthalpy = m_thermo->enthalpy_mass();
    m_intEnergy = m_thermo->intEnergy_mass();
    m_thermo->saveState(m_state);
}

void ConstPressureReactor::evalEqs(doublereal time, doublereal* y,
                                   doublereal* ydot, doublereal* params)
{
    double dmdt = 0.0; // dm/dt (gas phase)
    double* dYdt = ydot + 2;
    m_thermo->restoreState(m_state);
    applySensitivity(params);
    evalWalls(time);
    double mdot_surf = evalSurfaces(time, ydot + m_nsp + 2);
    dmdt += mdot_surf;


    //const vector_fp& mw = m_thermo->molecularWeights();
    const doublereal* Y = m_thermo->massFractions();
    


    //------------------------------------------------------------------------------------------------------------
    //Implemented by Jaouen to call an analytic reduced mechanism
   
    //jaouen Jaouen 
    //Here is where you should call your Reduced function
    //just after ydot[n+2] *= mw[n]; the values obtained represents the omega dot in kg/s

    //Adapted for Cantera 2.4 by Andrea Seltz - 18/12/18
    const doublereal* mw = DATA_PTR(m_thermo->molecularWeights());
    //Mass of each species
    //m_thermo->getMassFractions(y);
    //scale(y+2, y + m_nsp + 2 , mss, m_mass);
    
    doublereal rrho = 1.0/m_thermo->density();


    switch(Mechanism)
    {
       case 0 :
          fill(dYdt, dYdt + m_nsp, 0.0);
          Reduced(dYdt, m_thermo->temperature(), y+2, m_thermo->density(), mw );
          for (size_t k = 0; k < m_nsp; k++)
          {
           dYdt[k] *= mw[k]*rrho;     //           moles/s/m^3 -> moles/s
          }

       break;

       default :
          /*
           * Species equation 
           *  Equation is:
           *  \dot M_k = \hat W_k \dot\omega_k + \dot m_{in} Y_{k,in}
           *            - \dot m_{out} Y_{k} + A \dot s_k.
           *
          */


          if (m_chem) {
              m_kin->getNetProductionRates(&m_wdot[0]); // "omega dot"
          }
          else
          {
            fill(dYdt, dYdt + m_nsp, 0.0);
          }


          for (size_t k = 0; k < m_nsp; k++) {
              // production in gas phase and from surfaces
              dYdt[k] = (m_wdot[k] * m_vol + m_sdot[k]) * mw[k] / m_mass;
              // dilution by net surface mass flux
              dYdt[k] -= Y[k] * mdot_surf / m_mass;
          }

    }
//--------------------------------------------------------------------------------------------------------------

   /*
    *   Energy equation.
    *  \f[
    *  \dot U = -P\dot V + A \dot q + \dot m_{in} h_{in}
    *  - \dot m_{out} h.
    *  \f]
    */  
      
    // external heat transfer
    double dHdt = - m_Q;

    // add terms for outlets
    for (size_t i = 0; i < m_outlet.size(); i++) {
        double mdot_out = m_outlet[i]->massFlowRate(time); // mass flow out of system
        dmdt -= mdot_out;
        dHdt -= mdot_out * m_enthalpy;
    }

    // add terms for inlets
    for (size_t i = 0; i < m_inlet.size(); i++) {
        double mdot_in = m_inlet[i]->massFlowRate(time);
        dmdt += mdot_in; // mass flow into system
        for (size_t n = 0; n < m_nsp; n++) {
            double mdot_spec = m_inlet[i]->outletSpeciesMassFlowRate(n);
            // flow of species into system and dilution by other species
            dYdt[n] += (mdot_spec - mdot_in * Y[n]) / m_mass;
        }
        dHdt += mdot_in * m_inlet[i]->enthalpy_mass();
    }

    ydot[0] = dmdt;
    if (m_energy) {
        ydot[1] = dHdt;
    } else {
        ydot[1] = 0.0;
    }

    // reset sensitivity parameters
    resetSensitivity(params);
}

size_t ConstPressureReactor::componentIndex(const string& nm) const
{
    size_t k = speciesIndex(nm);
    if (k != npos) {
        return k + 2;
    } else if (nm == "mass") {
        return 0;
    } else if (nm == "enthalpy") {
        return 1;
    } else {
        return npos;
    }
}

std::string ConstPressureReactor::componentName(size_t k) {
    if (k == 0) {
        return "mass";
    } else if (k == 1) {
        return "enthalpy";
    } else if (k >= 2 && k < neq()) {
        k -= 2;
        if (k < m_thermo->nSpecies()) {
            return m_thermo->speciesName(k);
        } else {
            k -= m_thermo->nSpecies();
        }
        for (auto& S : m_surfaces) {
            ThermoPhase* th = S->thermo();
            if (k < th->nSpecies()) {
                return th->speciesName(k);
            } else {
                k -= th->nSpecies();
            }
        }
    }
    throw CanteraError("ConstPressureReactor::componentName",
                       "Index is out of bounds.");
}

}
