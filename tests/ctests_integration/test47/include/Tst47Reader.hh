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
#ifndef Tst47Reader_H
#define Tst47Reader_H

#include "TstDiscreteProcessReader.hh"

class Tst47Reader : public  TstDiscreteProcessReader
{

   public:
      
      // ctor & dtor
      Tst47Reader() : TstDiscreteProcessReader() /*, fAngleX(0.), fdAngleX(0.), fAngleY(0.), fdAngleY(0.), fdpByp(0.) */ {}
      virtual ~Tst47Reader() {}
      
//      G4double GetAngleX() const { return fAngleX; }
//      G4double GetAngleY() const { return fAngleY; }
//      G4double GetdpByp()  const { return fdpByp;  }
   
   protected:
   
      virtual void ProcessLine( G4String line );

   private:
   
//      G4double fAngleX;
//      G4double fdAngleX;
//      G4double fAngleY;
//      G4double fdAngleY;
//      G4double fdpByp; // don't know what it's for, will figure out later

};

#endif
