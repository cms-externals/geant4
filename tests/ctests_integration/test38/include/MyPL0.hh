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
//---------------------------------------------------------------------------
//
// ClassName:   MyPL0 - tst physics list
//
// Author: 2015-05-11 R. Hatcher
//
//----------------------------------------------------------------------------
//
#ifndef TMyPL0_h
#  define TMyPL0_h 1

#  include "G4VModularPhysicsList.hh"
#  include "globals.hh"

#  include "CompileTimeConstraints.hh"
#  include <CLHEP/Units/SystemOfUnits.h>

template<class T>
class TMyPL0 : public T
{
  public:

    TMyPL0(G4int ver = 1);
    virtual ~TMyPL0();

  public:

    // SetCuts()
    virtual void SetCuts();

  private:

    enum
    {
      ok = CompileTimeConstraints::IsA<T, G4VModularPhysicsList>::ok
    };
};

#  include "MyPL0.icc"
typedef TMyPL0<G4VModularPhysicsList> MyPL0;

#endif
