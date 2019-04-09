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
//
//
//

#ifndef Test2Run_h
#define Test2Run_h 1

#include "globals.hh"
#include "G4Run.hh"

#include "Test2RunSD.hh"
#include "Test2RunPS.hh"

class G4Event;

class Test2Run : public G4Run
{
public:
  Test2Run();
  virtual ~Test2Run();

public:
  virtual void RecordEvent(const G4Event*);

public:
  inline Test2RunSD* GetMassRunSD() const{ return MassRunSD; }
  inline Test2RunSD* GetParaRunSD() const{ return ParaRunSD; }
  inline Test2RunPS* GetMassRunPS() const{ return MassRunPS; }
  inline Test2RunPS* GetParaRunPS() const{ return ParaRunPS; }
  void DumpQuantitiesToFile();

private:
  Test2RunSD* MassRunSD;
  Test2RunSD* ParaRunSD;
  Test2RunPS* MassRunPS;
  Test2RunPS* ParaRunPS;

};
#endif
