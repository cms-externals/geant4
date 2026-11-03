//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//! @file G4Exception.hh
// ********************************************************************
// Authors: G.Cosmo, M.Asai - May 1999 - First implementation
// --------------------------------------------------------------------
#ifndef G4EXCEPTION_HH
#define G4EXCEPTION_HH

#include "G4ExceptionSeverity.hh"

#include <sstream>

/** @addtogroup global_management
 * @{
 */

//! @typedef Output stream type for forming G4ExceptionDescriptions.
using G4ExceptionDescription = std::ostringstream;

/**
 * Report a Geant4 exception.
 *
 * Forwards parameters to G4VExceptionHandler held by G4StateManager if one
 * has been installed. Otherwise, will print the information to G4cerr (or G4cout
 * for warnings). May abort the program depending on severity.
 *
 * @param originOfException Name of the function or location where the exception occurred.
 * @param exceptionCode Unique code identifying the exception.
 * @param severity Severity of the exception (see G4ExceptionSeverity).
 * @param description Description of the exception.
 */
void G4Exception(const char* originOfException, const char* exceptionCode,
                 G4ExceptionSeverity severity, const char* description);

//! @copydoc G4Exception
void G4Exception(const char* originOfException, const char* exceptionCode,
                 G4ExceptionSeverity severity, G4ExceptionDescription& description);

/**
 * @copydoc G4Exception
 * @param comments Additional comments to append to the description.
 */
void G4Exception(const char* originOfException, const char* exceptionCode,
                 G4ExceptionSeverity severity, G4ExceptionDescription& description,
                 const char* comments);

/** @} */
#endif
