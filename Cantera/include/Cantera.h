/**
 *  * @file Cantera.h
 *   *  Basic include file to be used in all Cantera application
 *    *  environments.
 *     */
/*
 * $Date: 2018-10-22 01:46:26 +0100 (Mon, 04 Jan 2010) $
 *  */

#ifndef CANTERA_H_INCL
#define CANTERA_H_INCL

// definitions
 #ifndef CANTERA_APP
 #define CANTERA_APP
 #endif
 



namespace Cantera_CXX{ }

using namespace Cantera_CXX;

#include "cantera/base/ct_defs.h"

 // some useful functions
 #include "cantera/base/global.h"

 // the CanteraError exception class
 #include "cantera/base/ctexceptions.h"

 //test
 //#include "cantera/equil/ChemEquil.h"

 //
 //#include "kernel/importCTML.h"

 // The Cantera logger class
 #include "cantera/base/logger.h"

 // Include the timer
 #include "cantera/base/clockWC.h"

 // Include routines for reading and writing XML files
 #include "cantera/base/xml.h"

 // Include string utility routines
 #include "cantera/base/stringUtils.h"

 #endif
 
