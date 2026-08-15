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
//! @file G4ApplicationState.hh
// ********************************************************************
#ifndef G4APPLICATIONSTATE_HH
#define G4APPLICATIONSTATE_HH

/**
 * @ingroup global_management
 * @enum G4ApplicationState
 * Specifies the state of the Geant4 application.
 *
 * The global/per-thread state is managed by @ref G4StateManager
 * with possible state transitions:
 *
 * ~~~
 *  PreInit
 *    |
 *    v
 *  Init
 *    |
 *    v
 *  Idle ------> Quit
 *   | ^
 *   v |
 *  GeomClosed (at each run)
 *   | ^
 *   v |
 *  EventProc (at each event)
 * ~~~
 *
 * @check Move descriptions of state transitions/behaviour to @ref G4StateManager?
 *
 * @var G4State_PreInit
 * At the very beginning of the Application G4StateManager starts
 * with this state. G4Initializer changes this state to Init when
 * G4Initializer::Initialize() method starts. At the moment of
 * the state change of PreInit->Init, no material, geometrical,
 * particle or physics process has been initialized.
 *
 * @var G4State_Init
 * State during the G4Initializer::Initialize() method.
 * G4Initializer changes this state to @ref G4State_Idle when all initialization procedures
 * are successfully Done.
 *
 * @var G4State_Idle
 * State when ready to start "BeamOn".
 * G4RunManager changes this state to GeomClosed when G4RunManager::BeamOn()
 * method starts and G4GeometryManager::CloseGeometry() is Done.
 * At the end of BeamOn() method, G4RunManager will reset the application state
 * to Idle after G4GeometryManager::OpenGeometry() is Done.
 *
 * @var G4State_GeomClosed
 * State between G4GeometryManager::CloseGeometry() and G4GeometryManager::OpenGeometry(),
 * but no event is in progress.
 * At the begining of each event (construction of a G4Event object and primary particle generation),
 * G4RunManager changes this state to EventProc and resets to GeomClosed state
 * when G4EventManager::ProcessOneEvent() is over.
 *
 * @var G4State_EventProc
 * State when processing an event.
 *
 * @var G4State_Quit
 * State when the destructor of G4RunManager is invoked.
 *
 * @var G4State_Abort
 * State when G4Exception is invoked.
 *
 */
enum G4ApplicationState
{
  G4State_PreInit,
  G4State_Init,
  G4State_Idle,
  G4State_GeomClosed,
  G4State_EventProc,
  G4State_Quit,
  G4State_Abort
};
#endif
