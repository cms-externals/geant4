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
// 
// ----------------------------------------------------------------------
#include "G4Timer.hh"
#include "globals.hh"
#include "G4ios.hh"

int main()
{
  G4Timer timer;
  G4double j=1.0;
  G4cout.precision(20);

  timer.Start();
  G4cout << "Result: ";
  for (size_t i=0; i<100000000; i++) j=(i+1)*(j+1)/(5+j*i);
  G4cout << j << G4endl;
  timer.Stop();  

  G4cout << "System time: " << timer.GetSystemElapsed() << G4endl
         << "User time:   " << timer.GetUserElapsed() << G4endl
         << "Real time:   " << timer.GetRealElapsed() << G4endl;
  G4cout << "Clock time:  " << timer.GetClockTime() << G4endl;

  return 0;
}
