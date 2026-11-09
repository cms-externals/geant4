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
//! @file G4StateManager.hh
// ********************************************************************
// Authors: G.Cosmo, M.Asai - November 1996
// --------------------------------------------------------------------
#ifndef G4STATEMANAGER_HH
#define G4STATEMANAGER_HH

#include "G4ApplicationState.hh"
#include "G4String.hh"
#include "G4Types.hh"
#include "G4VExceptionHandler.hh"
#include "G4VStateDependent.hh"

#include <vector>

class G4Run;
class G4Event;

/**
 * @ingroup global_management
 * Handles and updates the running Geant4 application state.
 *
 * The class is a singleton, it can be accessed via the public
 * method G4StateManager::GetStateManager().
 *
 * @check Document use of G4VStateDependent and G4VExceptionHandler.
 *
 * @sa G4ApplicationState for possible states.
 */
class G4StateManager
{
  public:

    /** Returns the instance on the callee's thread.
     *
     * @post Returned pointer is not null.
     */
    static G4StateManager* GetStateManager();

    ~G4StateManager();

    G4StateManager(const G4StateManager&) = delete;
    G4StateManager& operator=(const G4StateManager&) = delete;
    G4bool operator==(const G4StateManager&) const = delete;
    G4bool operator!=(const G4StateManager&) const = delete;

    //! Returns the current state
    const G4ApplicationState& GetCurrentState() const;

    //! Returns the previous state
    const G4ApplicationState& GetPreviousState() const;

    /**
     * Set Geant4 to a new state.
     *
     * In case the request is illegal, false will be returned
     * and the state of Geant4 will not be changed.
     */
    G4bool SetNewState(const G4ApplicationState& requestedState);

    /**
     * Set Geant4 to a new state.
     *
     * In case the request is illegal, false will be returned
     * and the state of Geant4 will not be changed.
     * "msg" is the information associated to the state change
     */
    G4bool SetNewState(const G4ApplicationState& requestedState, const char* msg);

    /**
     * Register a concrete class of G4VStateDependent.
     * Registered concrete classes will be notified via
     * G4VStateDependent::Notify() method when the state of Geant4 changes.
     *
     * @return `true` if registration successful.
     */
    G4bool RegisterDependent(G4VStateDependent* aDependent, G4bool bottom = false);

    /**
     * Remove the registration.
     *
     * @return `false` if @p aDependent is not already registered.
     */
    G4bool DeregisterDependent(G4VStateDependent* aDependent);

    /** Remove the registration.
     *
     * Removed pointer is returned
     */
    G4VStateDependent* RemoveDependent(const G4VStateDependent* aDependent);

    //! Return string of the state name
    G4String GetStateString(const G4ApplicationState& aState) const;

    /**
     * Notify observers of event deletion.
     */
    void NotifyDeletion(const G4Event*);

    /**
     * Notify observers of run deletion.
     *
     * Notifying when an event or a run is deleted for the sake of avoiding
     * a state-dependent class from accessing to an obsolete event/run object.
     */
    void NotifyDeletion(const G4Run*);

    inline void SetSuppressAbortion(G4int i);
    inline G4int GetSuppressAbortion() const;
    inline const char* GetMessage() const;

    //! Install handler for G4Exception.
    inline void SetExceptionHandler(G4VExceptionHandler* eh);

    //! Return pointer to installed handler for G4Exception.
    inline G4VExceptionHandler* GetExceptionHandler() const;

    static void SetVerboseLevel(G4int val);

  private:

    G4StateManager();

  private:

    static G4ThreadLocal G4StateManager* theStateManager;
    G4ApplicationState theCurrentState = G4State_PreInit;
    G4ApplicationState thePreviousState = G4State_PreInit;
    std::vector<G4VStateDependent*> theDependentsList;
    G4VStateDependent* theBottomDependent = nullptr;
    G4int suppressAbortion = 0;
    const char* msgptr = nullptr;
    G4VExceptionHandler* exceptionHandler = nullptr;
    static G4int verboseLevel;
};

#include "G4StateManager.icc"

#endif
