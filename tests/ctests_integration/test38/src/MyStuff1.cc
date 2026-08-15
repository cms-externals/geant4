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
// MyStuff1() function
// registers physics lists:  MyPL1
//
// Author: 2015-05-11 R. Hatcher
//
// Modified:
//
//----------------------------------------------------------------------------
//
#include "G4PhysListStamper.hh"

#include "MyPL1.hh"

#include <iostream>
G4_DECLARE_PHYSLIST_FACTORY(MyPL1);

int MyStuff1()
{
  // need to make a library with some code in it
  std::cout << "MyStuff1() was called " << std::endl;
  return 1;
}

// static int globalMyStuff1 = MyStuff1();
