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
//! @file G4ExceptionSeverity.hh
// ********************************************************************
// Author: M.Asai, 19 August 2002
// --------------------------------------------------------------------
#ifndef G4EXCEPTIONSEVERITY_HH
#define G4EXCEPTIONSEVERITY_HH

/**
 * @ingroup global_management
 * @enum G4ExceptionSeverity
 * @brief Specifies the severity of a G4Exception
 *
 * @var FatalException
 *   Error is severe or it happens at the initialization time.
 *   Program should be aborted and core dump will be generated.
 *
 * @var FatalErrorInArgument
 *   Fatal error caused by most likely the mis-use of interfaces
 *   by the user's code. Program should be aborted and core dump
 *   will be generated.
 *
 * @var RunMustBeAborted
 *   Error happens at initialization of a run. This might be at the
 *   moment of closing geometry, or some unpleasant situation
 *   occurs during the event loop. Current run will be aborted
 *   and the application returns to "Idle" state.
 *
 * @var EventMustBeAborted
 *   Error happens during tracking a particle. The event currently
 *   being processed should be aborted, run will not be aborted.
 *
 * @var JustWarning
 *   Just display messages.
 *
 * @var IgnoreTheIssue
 *   No message generated.
 */
enum G4ExceptionSeverity
{
  FatalException,
  FatalErrorInArgument,
  RunMustBeAborted,
  EventMustBeAborted,
  JustWarning,
  IgnoreTheIssue
};
#endif
