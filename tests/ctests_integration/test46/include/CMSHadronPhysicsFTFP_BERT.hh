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
//---------------------------------------------------------------------------
// Author: Vladimir Ivanchenko
// Date:   March 2018
//
// Hadron physics for the new CMS physics list FTFP_BERT_EMM_TRK.
// The hadron physics of FTFP_BERT has the transition between Bertini
// (BERT) intra-nuclear cascade model and Fritiof (FTF) string model in the
// energy region [4, 5] GeV (instead of the default for Geant4 10.4).
//---------------------------------------------------------------------------
//
#ifndef SimG4Core_PhysicsLists_CMSHadronPhysicsFTFP_BERT_h
#define SimG4Core_PhysicsLists_CMSHadronPhysicsFTFP_BERT_h

#include "G4HadronPhysicsFTFP_BERT.hh"
#include "globals.hh"

class CMSHadronPhysicsFTFP_BERT : public G4HadronPhysicsFTFP_BERT
{
  public:

    explicit CMSHadronPhysicsFTFP_BERT(G4int verb);
    explicit CMSHadronPhysicsFTFP_BERT(G4double e1, G4double e2, G4double e3);
    ~CMSHadronPhysicsFTFP_BERT() override = default;

    // copy constructor and hide assignment operator
    CMSHadronPhysicsFTFP_BERT(CMSHadronPhysicsFTFP_BERT&) = delete;
    CMSHadronPhysicsFTFP_BERT& operator=(const CMSHadronPhysicsFTFP_BERT& right) = delete;

  protected:

    void DumpBanner() override;
};

#endif
