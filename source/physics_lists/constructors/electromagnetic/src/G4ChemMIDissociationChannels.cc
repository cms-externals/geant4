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
// G4ChemMIDissociationChannels.cc
//
// Created at 2026/1/19 (Mon.)
// Author: Shogo OKADA @KEK-CRC (shogo.okada@kek.jp)
//
// Reference: J.Meesungnoen et. al, DOI: 10.1021/jp058037z
//

#include "G4ChemMIDissociationChannels.hh"

#include "G4DNAWaterDissociationDisplacer.hh"
#include "G4H2O.hh"
#include "G4H3O.hh"
#include "G4Hydrogen.hh"
#include "G4MolecularConfiguration.hh"
#include "G4MoleculeTable.hh"
#include "G4OH.hh"
#include "G4Oxygen.hh"

void G4ChemMIDissociationChannels::ConstructDissociationChannels()
{
  // Get the molecular configuration
  auto molTable = G4MoleculeTable::Instance();
  G4MolecularConfiguration* OH = molTable->GetConfiguration("°OH");
  G4MolecularConfiguration* H3O = molTable->GetConfiguration("H3Op");
  G4MolecularConfiguration* H = molTable->GetConfiguration("H");
  G4MolecularConfiguration* O3P = molTable->GetConfiguration("Oxy");

  // Define the decay channels
  G4MoleculeDefinition* water = G4H2O::Definition();
  G4MolecularDissociationChannel* decCh1{nullptr};
  G4MolecularDissociationChannel* decCh2{nullptr};
  G4MolecularDissociationChannel* decCh3{nullptr};

  G4ElectronOccupancy* occ = new G4ElectronOccupancy(*(water->GetGroundStateElectronOccupancy()));

  //////////////////////////////////////////////////////////////////////////////
  // Double ionisation
  //////////////////////////////////////////////////////////////////////////////
  decCh1 = new G4MolecularDissociationChannel("DoubleIonisation_Channel1");
  decCh2 = new G4MolecularDissociationChannel("DoubleIonisation_Channel2");
  decCh3 = new G4MolecularDissociationChannel("DoubleIonisation_Channel3");

  // Decay Channel #1: H2O^2+ -> 2H+ + O(3P) -> 2H3O+ + O(3P)
  decCh1->AddProduct(H3O);
  decCh1->AddProduct(H3O);
  decCh1->AddProduct(O3P);
  decCh1->SetProbability(0.29);
  decCh1->SetDisplacementType(G4DNAWaterDissociationDisplacer::DoubleIonisation_DissociationDecay1);

  // Decay Channel #2: H2O^2+ -> H+ + H* + O+ -> 2H3O+ + H* + *OH + O(3P)
  decCh2->AddProduct(H3O);
  decCh2->AddProduct(H3O);
  decCh2->AddProduct(H);
  decCh2->AddProduct(OH);
  decCh2->AddProduct(O3P);
  decCh2->SetProbability(0.16);
  decCh2->SetDisplacementType(G4DNAWaterDissociationDisplacer::DoubleIonisation_DissociationDecay2);

  // Decay Channel #3: H2O^2+ -> H+ + OH+ -> 2H3O+ + O(3P)
  decCh3->AddProduct(H3O);
  decCh3->AddProduct(H3O);
  decCh3->AddProduct(O3P);
  decCh3->SetProbability(0.55);
  decCh3->SetDisplacementType(G4DNAWaterDissociationDisplacer::DoubleIonisation_DissociationDecay3);

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(4, 2);
  water->NewConfigurationWithElectronOccupancy("DoubleIonisation1", *occ);
  water->AddDecayChannel("DoubleIonisation1", decCh1);
  water->AddDecayChannel("DoubleIonisation1", decCh2);
  water->AddDecayChannel("DoubleIonisation1", decCh3);

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(3, 1);
  occ->RemoveElectron(4, 1);
  water->NewConfigurationWithElectronOccupancy("DoubleIonisation2", *occ);
  water->AddDecayChannel("DoubleIonisation2", new G4MolecularDissociationChannel(*decCh1));
  water->AddDecayChannel("DoubleIonisation2", new G4MolecularDissociationChannel(*decCh2));
  water->AddDecayChannel("DoubleIonisation2", new G4MolecularDissociationChannel(*decCh3));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(3, 2);
  water->NewConfigurationWithElectronOccupancy("DoubleIonisation3", *occ);
  water->AddDecayChannel("DoubleIonisation3", new G4MolecularDissociationChannel(*decCh1));
  water->AddDecayChannel("DoubleIonisation3", new G4MolecularDissociationChannel(*decCh2));
  water->AddDecayChannel("DoubleIonisation3", new G4MolecularDissociationChannel(*decCh3));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(2, 1);
  occ->RemoveElectron(4, 1);
  water->NewConfigurationWithElectronOccupancy("DoubleIonisation4", *occ);
  water->AddDecayChannel("DoubleIonisation4", new G4MolecularDissociationChannel(*decCh1));
  water->AddDecayChannel("DoubleIonisation4", new G4MolecularDissociationChannel(*decCh2));
  water->AddDecayChannel("DoubleIonisation4", new G4MolecularDissociationChannel(*decCh3));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(2, 1);
  occ->RemoveElectron(3, 1);
  water->NewConfigurationWithElectronOccupancy("DoubleIonisation5", *occ);
  water->AddDecayChannel("DoubleIonisation5", new G4MolecularDissociationChannel(*decCh1));
  water->AddDecayChannel("DoubleIonisation5", new G4MolecularDissociationChannel(*decCh2));
  water->AddDecayChannel("DoubleIonisation5", new G4MolecularDissociationChannel(*decCh3));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(2, 2);
  water->NewConfigurationWithElectronOccupancy("DoubleIonisation6", *occ);
  water->AddDecayChannel("DoubleIonisation6", new G4MolecularDissociationChannel(*decCh1));
  water->AddDecayChannel("DoubleIonisation6", new G4MolecularDissociationChannel(*decCh2));
  water->AddDecayChannel("DoubleIonisation6", new G4MolecularDissociationChannel(*decCh3));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(1, 1);
  occ->RemoveElectron(4, 1);
  water->NewConfigurationWithElectronOccupancy("DoubleIonisation7", *occ);
  water->AddDecayChannel("DoubleIonisation7", new G4MolecularDissociationChannel(*decCh1));
  water->AddDecayChannel("DoubleIonisation7", new G4MolecularDissociationChannel(*decCh2));
  water->AddDecayChannel("DoubleIonisation7", new G4MolecularDissociationChannel(*decCh3));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(1, 1);
  occ->RemoveElectron(3, 1);
  water->NewConfigurationWithElectronOccupancy("DoubleIonisation8", *occ);
  water->AddDecayChannel("DoubleIonisation8", new G4MolecularDissociationChannel(*decCh1));
  water->AddDecayChannel("DoubleIonisation8", new G4MolecularDissociationChannel(*decCh2));
  water->AddDecayChannel("DoubleIonisation8", new G4MolecularDissociationChannel(*decCh3));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(1, 1);
  occ->RemoveElectron(2, 1);
  water->NewConfigurationWithElectronOccupancy("DoubleIonisation9", *occ);
  water->AddDecayChannel("DoubleIonisation9", new G4MolecularDissociationChannel(*decCh1));
  water->AddDecayChannel("DoubleIonisation9", new G4MolecularDissociationChannel(*decCh2));
  water->AddDecayChannel("DoubleIonisation9", new G4MolecularDissociationChannel(*decCh3));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(1, 2);
  water->NewConfigurationWithElectronOccupancy("DoubleIonisation10", *occ);
  water->AddDecayChannel("DoubleIonisation10", new G4MolecularDissociationChannel(*decCh1));
  water->AddDecayChannel("DoubleIonisation10", new G4MolecularDissociationChannel(*decCh2));
  water->AddDecayChannel("DoubleIonisation10", new G4MolecularDissociationChannel(*decCh3));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(0, 1);
  occ->RemoveElectron(4, 1);
  water->NewConfigurationWithElectronOccupancy("DoubleIonisation11", *occ);
  water->AddDecayChannel("DoubleIonisation11", new G4MolecularDissociationChannel(*decCh1));
  water->AddDecayChannel("DoubleIonisation11", new G4MolecularDissociationChannel(*decCh2));
  water->AddDecayChannel("DoubleIonisation11", new G4MolecularDissociationChannel(*decCh3));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(0, 1);
  occ->RemoveElectron(3, 1);
  water->NewConfigurationWithElectronOccupancy("DoubleIonisation12", *occ);
  water->AddDecayChannel("DoubleIonisation12", new G4MolecularDissociationChannel(*decCh1));
  water->AddDecayChannel("DoubleIonisation12", new G4MolecularDissociationChannel(*decCh2));
  water->AddDecayChannel("DoubleIonisation12", new G4MolecularDissociationChannel(*decCh3));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(0, 1);
  occ->RemoveElectron(2, 1);
  water->NewConfigurationWithElectronOccupancy("DoubleIonisation13", *occ);
  water->AddDecayChannel("DoubleIonisation13", new G4MolecularDissociationChannel(*decCh1));
  water->AddDecayChannel("DoubleIonisation13", new G4MolecularDissociationChannel(*decCh2));
  water->AddDecayChannel("DoubleIonisation13", new G4MolecularDissociationChannel(*decCh3));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(0, 1);
  occ->RemoveElectron(1, 1);
  water->NewConfigurationWithElectronOccupancy("DoubleIonisation14", *occ);
  water->AddDecayChannel("DoubleIonisation14", new G4MolecularDissociationChannel(*decCh1));
  water->AddDecayChannel("DoubleIonisation14", new G4MolecularDissociationChannel(*decCh2));
  water->AddDecayChannel("DoubleIonisation14", new G4MolecularDissociationChannel(*decCh3));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(0, 2);
  water->NewConfigurationWithElectronOccupancy("DoubleIonisation15", *occ);
  water->AddDecayChannel("DoubleIonisation15", new G4MolecularDissociationChannel(*decCh1));
  water->AddDecayChannel("DoubleIonisation15", new G4MolecularDissociationChannel(*decCh2));
  water->AddDecayChannel("DoubleIonisation15", new G4MolecularDissociationChannel(*decCh3));

  //////////////////////////////////////////////////////////////////////////////
  // Triple ionisation
  //////////////////////////////////////////////////////////////////////////////
  decCh1 = new G4MolecularDissociationChannel("TripleIonisation_Channel");

  // Decay Channel 1 : H2O^3+ -> 3H_3Op + OH + O(3P)
  decCh1->AddProduct(H3O);
  decCh1->AddProduct(H3O);
  decCh1->AddProduct(H3O);
  decCh1->AddProduct(OH);
  decCh1->AddProduct(O3P);
  decCh1->SetProbability(1.0);
  decCh1->SetDisplacementType(G4DNAWaterDissociationDisplacer::TripleIonisation_DissociationDecay);

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(4, 1);
  occ->RemoveElectron(3, 1);
  occ->RemoveElectron(2, 1);
  water->NewConfigurationWithElectronOccupancy("TripleIonisation1", *occ);
  water->AddDecayChannel("TripleIonisation1", decCh1);

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(4, 1);
  occ->RemoveElectron(3, 1);
  occ->RemoveElectron(1, 1);
  water->NewConfigurationWithElectronOccupancy("TripleIonisation2", *occ);
  water->AddDecayChannel("TripleIonisation2", new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(4, 1);
  occ->RemoveElectron(2, 1);
  occ->RemoveElectron(1, 1);
  water->NewConfigurationWithElectronOccupancy("TripleIonisation3", *occ);
  water->AddDecayChannel("TripleIonisation3", new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(3, 1);
  occ->RemoveElectron(2, 1);
  occ->RemoveElectron(1, 1);
  water->NewConfigurationWithElectronOccupancy("TripleIonisation4", *occ);
  water->AddDecayChannel("TripleIonisation4", new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(4, 1);
  occ->RemoveElectron(3, 1);
  occ->RemoveElectron(0, 1);
  water->NewConfigurationWithElectronOccupancy("TripleIonisation5", *occ);
  water->AddDecayChannel("TripleIonisation5", new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(4, 1);
  occ->RemoveElectron(2, 1);
  occ->RemoveElectron(0, 1);
  water->NewConfigurationWithElectronOccupancy("TripleIonisation6", *occ);
  water->AddDecayChannel("TripleIonisation6", new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(3, 1);
  occ->RemoveElectron(2, 1);
  occ->RemoveElectron(0, 1);
  water->NewConfigurationWithElectronOccupancy("TripleIonisation7", *occ);
  water->AddDecayChannel("TripleIonisation7", new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(4, 1);
  occ->RemoveElectron(1, 1);
  occ->RemoveElectron(0, 1);
  water->NewConfigurationWithElectronOccupancy("TripleIonisation8", *occ);
  water->AddDecayChannel("TripleIonisation8", new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(3, 1);
  occ->RemoveElectron(1, 1);
  occ->RemoveElectron(0, 1);
  water->NewConfigurationWithElectronOccupancy("TripleIonisation9", *occ);
  water->AddDecayChannel("TripleIonisation9", new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(2, 1);
  occ->RemoveElectron(1, 1);
  occ->RemoveElectron(0, 1);
  water->NewConfigurationWithElectronOccupancy("TripleIonisation10", *occ);
  water->AddDecayChannel("TripleIonisation10", new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(4, 2);
  occ->RemoveElectron(3, 1);
  water->NewConfigurationWithElectronOccupancy("TripleIonisation11", *occ);
  water->AddDecayChannel("TripleIonisation11", new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(4, 2);
  occ->RemoveElectron(2, 1);
  water->NewConfigurationWithElectronOccupancy("TripleIonisation12", *occ);
  water->AddDecayChannel("TripleIonisation12", new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(4, 2);
  occ->RemoveElectron(1, 1);
  water->NewConfigurationWithElectronOccupancy("TripleIonisation13", *occ);
  water->AddDecayChannel("TripleIonisation13", new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(4, 2);
  occ->RemoveElectron(0, 1);
  water->NewConfigurationWithElectronOccupancy("TripleIonisation14", *occ);
  water->AddDecayChannel("TripleIonisation14", new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(3, 2);
  occ->RemoveElectron(4, 1);
  water->NewConfigurationWithElectronOccupancy("TripleIonisation15", *occ);
  water->AddDecayChannel("TripleIonisation15", new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(3, 2);
  occ->RemoveElectron(2, 1);
  water->NewConfigurationWithElectronOccupancy("TripleIonisation16", *occ);
  water->AddDecayChannel("TripleIonisation16", new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(3, 2);
  occ->RemoveElectron(1, 1);
  water->NewConfigurationWithElectronOccupancy("TripleIonisation17", *occ);
  water->AddDecayChannel("TripleIonisation17", new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(3, 2);
  occ->RemoveElectron(0, 1);
  water->NewConfigurationWithElectronOccupancy("TripleIonisation18", *occ);
  water->AddDecayChannel("TripleIonisation18", new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(2, 2);
  occ->RemoveElectron(4, 1);
  water->NewConfigurationWithElectronOccupancy("TripleIonisation19", *occ);
  water->AddDecayChannel("TripleIonisation19", new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(2, 2);
  occ->RemoveElectron(3, 1);
  water->NewConfigurationWithElectronOccupancy("TripleIonisation20", *occ);
  water->AddDecayChannel("TripleIonisation20", new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(2, 2);
  occ->RemoveElectron(1, 1);
  water->NewConfigurationWithElectronOccupancy("TripleIonisation21", *occ);
  water->AddDecayChannel("TripleIonisation21", new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(2, 2);
  occ->RemoveElectron(0, 1);
  water->NewConfigurationWithElectronOccupancy("TripleIonisation22", *occ);
  water->AddDecayChannel("TripleIonisation22", new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(1, 2);
  occ->RemoveElectron(4, 1);
  water->NewConfigurationWithElectronOccupancy("TripleIonisation23", *occ);
  water->AddDecayChannel("TripleIonisation23", new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(1, 2);
  occ->RemoveElectron(3, 1);
  water->NewConfigurationWithElectronOccupancy("TripleIonisation24", *occ);
  water->AddDecayChannel("TripleIonisation24", new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(1, 2);
  occ->RemoveElectron(2, 1);
  water->NewConfigurationWithElectronOccupancy("TripleIonisation25", *occ);
  water->AddDecayChannel("TripleIonisation25", new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(1, 2);
  occ->RemoveElectron(0, 1);
  water->NewConfigurationWithElectronOccupancy("TripleIonisation26", *occ);
  water->AddDecayChannel("TripleIonisation26", new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(0, 2);
  occ->RemoveElectron(4, 1);
  water->NewConfigurationWithElectronOccupancy("TripleIonisation27", *occ);
  water->AddDecayChannel("TripleIonisation27", new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(0, 2);
  occ->RemoveElectron(3, 1);
  water->NewConfigurationWithElectronOccupancy("TripleIonisation28", *occ);
  water->AddDecayChannel("TripleIonisation28", new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(0, 2);
  occ->RemoveElectron(2, 1);
  water->NewConfigurationWithElectronOccupancy("TripleIonisation29", *occ);
  water->AddDecayChannel("TripleIonisation29", new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(0, 2);
  occ->RemoveElectron(1, 1);
  water->NewConfigurationWithElectronOccupancy("TripleIonisation30", *occ);
  water->AddDecayChannel("TripleIonisation30", new G4MolecularDissociationChannel(*decCh1));

  //////////////////////////////////////////////////////////////////////////////
  // Quadruple ionisation
  //////////////////////////////////////////////////////////////////////////////
  decCh1 = new G4MolecularDissociationChannel("QuadrupleIonisation_Channel");

  // Decay Channel 1 : H2O^4+ -> 4H_3Op + 2OH + O(3P)
  decCh1->AddProduct(H3O);
  decCh1->AddProduct(H3O);
  decCh1->AddProduct(H3O);
  decCh1->AddProduct(H3O);
  decCh1->AddProduct(OH);
  decCh1->AddProduct(OH);
  decCh1->AddProduct(O3P);
  decCh1->SetProbability(1.0);
  decCh1->SetDisplacementType(
    G4DNAWaterDissociationDisplacer::QuadrupleIonisation_DissociationDecay);

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(0, 1);
  occ->RemoveElectron(1, 1);
  occ->RemoveElectron(2, 1);
  occ->RemoveElectron(3, 1);
  water->NewConfigurationWithElectronOccupancy("QuadrupleIonisation1", *occ);
  water->AddDecayChannel("QuadrupleIonisation1", decCh1);

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(0, 1);
  occ->RemoveElectron(1, 1);
  occ->RemoveElectron(2, 1);
  occ->RemoveElectron(4, 1);
  water->NewConfigurationWithElectronOccupancy("QuadrupleIonisation2", *occ);
  water->AddDecayChannel("QuadrupleIonisation2", new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(0, 1);
  occ->RemoveElectron(1, 1);
  occ->RemoveElectron(3, 1);
  occ->RemoveElectron(4, 1);
  water->NewConfigurationWithElectronOccupancy("QuadrupleIonisation3", *occ);
  water->AddDecayChannel("QuadrupleIonisation3", new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(0, 1);
  occ->RemoveElectron(2, 1);
  occ->RemoveElectron(3, 1);
  occ->RemoveElectron(4, 1);
  water->NewConfigurationWithElectronOccupancy("QuadrupleIonisation4", *occ);
  water->AddDecayChannel("QuadrupleIonisation4", new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(1, 1);
  occ->RemoveElectron(2, 1);
  occ->RemoveElectron(3, 1);
  occ->RemoveElectron(4, 1);
  water->NewConfigurationWithElectronOccupancy("QuadrupleIonisation5", *occ);
  water->AddDecayChannel("QuadrupleIonisation5", new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(0, 1);
  occ->RemoveElectron(1, 1);
  occ->RemoveElectron(2, 2);
  water->NewConfigurationWithElectronOccupancy("QuadrupleIonisation6", *occ);
  water->AddDecayChannel("QuadrupleIonisation6", new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(0, 1);
  occ->RemoveElectron(1, 1);
  occ->RemoveElectron(3, 2);
  water->NewConfigurationWithElectronOccupancy("QuadrupleIonisation7", *occ);
  water->AddDecayChannel("QuadrupleIonisation7", new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(0, 1);
  occ->RemoveElectron(1, 1);
  occ->RemoveElectron(4, 2);
  water->NewConfigurationWithElectronOccupancy("QuadrupleIonisation8", *occ);
  water->AddDecayChannel("QuadrupleIonisation8", new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(0, 1);
  occ->RemoveElectron(2, 1);
  occ->RemoveElectron(1, 2);
  water->NewConfigurationWithElectronOccupancy("QuadrupleIonisation9", *occ);
  water->AddDecayChannel("QuadrupleIonisation9", new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(0, 1);
  occ->RemoveElectron(2, 1);
  occ->RemoveElectron(3, 2);
  water->NewConfigurationWithElectronOccupancy("QuadrupleIonisation10", *occ);
  water->AddDecayChannel("QuadrupleIonisation10", new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(0, 1);
  occ->RemoveElectron(2, 1);
  occ->RemoveElectron(4, 2);
  water->NewConfigurationWithElectronOccupancy("QuadrupleIonisation11", *occ);
  water->AddDecayChannel("QuadrupleIonisation11", new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(0, 1);
  occ->RemoveElectron(3, 1);
  occ->RemoveElectron(1, 2);
  water->NewConfigurationWithElectronOccupancy("QuadrupleIonisation12", *occ);
  water->AddDecayChannel("QuadrupleIonisation12", new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(0, 1);
  occ->RemoveElectron(3, 1);
  occ->RemoveElectron(2, 2);
  water->NewConfigurationWithElectronOccupancy("QuadrupleIonisation13", *occ);
  water->AddDecayChannel("QuadrupleIonisation13", new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(0, 1);
  occ->RemoveElectron(3, 1);
  occ->RemoveElectron(4, 2);
  water->NewConfigurationWithElectronOccupancy("QuadrupleIonisation14", *occ);
  water->AddDecayChannel("QuadrupleIonisation14", new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(0, 1);
  occ->RemoveElectron(4, 1);
  occ->RemoveElectron(1, 2);
  water->NewConfigurationWithElectronOccupancy("QuadrupleIonisation15", *occ);
  water->AddDecayChannel("QuadrupleIonisation15", new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(0, 1);
  occ->RemoveElectron(4, 1);
  occ->RemoveElectron(2, 2);
  water->NewConfigurationWithElectronOccupancy("QuadrupleIonisation16", *occ);
  water->AddDecayChannel("QuadrupleIonisation16", new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(0, 1);
  occ->RemoveElectron(4, 1);
  occ->RemoveElectron(3, 2);
  water->NewConfigurationWithElectronOccupancy("QuadrupleIonisation17", *occ);
  water->AddDecayChannel("QuadrupleIonisation17", new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(1, 1);
  occ->RemoveElectron(2, 1);
  occ->RemoveElectron(0, 2);
  water->NewConfigurationWithElectronOccupancy("QuadrupleIonisation18", *occ);
  water->AddDecayChannel("QuadrupleIonisation18", new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(1, 1);
  occ->RemoveElectron(2, 1);
  occ->RemoveElectron(3, 2);
  water->NewConfigurationWithElectronOccupancy("QuadrupleIonisation19", *occ);
  water->AddDecayChannel("QuadrupleIonisation19", new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(1, 1);
  occ->RemoveElectron(2, 1);
  occ->RemoveElectron(4, 2);
  water->NewConfigurationWithElectronOccupancy("QuadrupleIonisation20", *occ);
  water->AddDecayChannel("QuadrupleIonisation20", new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(1, 1);
  occ->RemoveElectron(3, 1);
  occ->RemoveElectron(0, 2);
  water->NewConfigurationWithElectronOccupancy("QuadrupleIonisation21", *occ);
  water->AddDecayChannel("QuadrupleIonisation21", new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(1, 1);
  occ->RemoveElectron(3, 1);
  occ->RemoveElectron(2, 2);
  water->NewConfigurationWithElectronOccupancy("QuadrupleIonisation22", *occ);
  water->AddDecayChannel("QuadrupleIonisation22", new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(1, 1);
  occ->RemoveElectron(3, 1);
  occ->RemoveElectron(4, 2);
  water->NewConfigurationWithElectronOccupancy("QuadrupleIonisation23", *occ);
  water->AddDecayChannel("QuadrupleIonisation23", new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(1, 1);
  occ->RemoveElectron(4, 1);
  occ->RemoveElectron(0, 2);
  water->NewConfigurationWithElectronOccupancy("QuadrupleIonisation24", *occ);
  water->AddDecayChannel("QuadrupleIonisation24", new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(1, 1);
  occ->RemoveElectron(4, 1);
  occ->RemoveElectron(2, 2);
  water->NewConfigurationWithElectronOccupancy("QuadrupleIonisation25", *occ);
  water->AddDecayChannel("QuadrupleIonisation25", new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(1, 1);
  occ->RemoveElectron(4, 1);
  occ->RemoveElectron(3, 2);
  water->NewConfigurationWithElectronOccupancy("QuadrupleIonisation26", *occ);
  water->AddDecayChannel("QuadrupleIonisation26", new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(2, 1);
  occ->RemoveElectron(3, 1);
  occ->RemoveElectron(0, 2);
  water->NewConfigurationWithElectronOccupancy("QuadrupleIonisation27", *occ);
  water->AddDecayChannel("QuadrupleIonisation27", new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(2, 1);
  occ->RemoveElectron(3, 1);
  occ->RemoveElectron(1, 2);
  water->NewConfigurationWithElectronOccupancy("QuadrupleIonisation28", *occ);
  water->AddDecayChannel("QuadrupleIonisation28", new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(2, 1);
  occ->RemoveElectron(3, 1);
  occ->RemoveElectron(4, 2);
  water->NewConfigurationWithElectronOccupancy("QuadrupleIonisation29", *occ);
  water->AddDecayChannel("QuadrupleIonisation29", new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(2, 1);
  occ->RemoveElectron(4, 1);
  occ->RemoveElectron(0, 2);
  water->NewConfigurationWithElectronOccupancy("QuadrupleIonisation30", *occ);
  water->AddDecayChannel("QuadrupleIonisation30", new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(2, 1);
  occ->RemoveElectron(4, 1);
  occ->RemoveElectron(1, 2);
  water->NewConfigurationWithElectronOccupancy("QuadrupleIonisation31", *occ);
  water->AddDecayChannel("QuadrupleIonisation31", new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(2, 1);
  occ->RemoveElectron(4, 1);
  occ->RemoveElectron(3, 2);
  water->NewConfigurationWithElectronOccupancy("QuadrupleIonisation32", *occ);
  water->AddDecayChannel("QuadrupleIonisation32", new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(3, 1);
  occ->RemoveElectron(4, 1);
  occ->RemoveElectron(0, 2);
  water->NewConfigurationWithElectronOccupancy("QuadrupleIonisation33", *occ);
  water->AddDecayChannel("QuadrupleIonisation33", new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(3, 1);
  occ->RemoveElectron(4, 1);
  occ->RemoveElectron(1, 2);
  water->NewConfigurationWithElectronOccupancy("QuadrupleIonisation34", *occ);
  water->AddDecayChannel("QuadrupleIonisation34", new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(3, 1);
  occ->RemoveElectron(4, 1);
  occ->RemoveElectron(2, 2);
  water->NewConfigurationWithElectronOccupancy("QuadrupleIonisation35", *occ);
  water->AddDecayChannel("QuadrupleIonisation35", new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(0, 2);
  occ->RemoveElectron(1, 2);
  water->NewConfigurationWithElectronOccupancy("QuadrupleIonisation36", *occ);
  water->AddDecayChannel("QuadrupleIonisation36", new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(0, 2);
  occ->RemoveElectron(2, 2);
  water->NewConfigurationWithElectronOccupancy("QuadrupleIonisation37", *occ);
  water->AddDecayChannel("QuadrupleIonisation37", new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(0, 2);
  occ->RemoveElectron(3, 2);
  water->NewConfigurationWithElectronOccupancy("QuadrupleIonisation38", *occ);
  water->AddDecayChannel("QuadrupleIonisation38", new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(0, 2);
  occ->RemoveElectron(4, 2);
  water->NewConfigurationWithElectronOccupancy("QuadrupleIonisation39", *occ);
  water->AddDecayChannel("QuadrupleIonisation39", new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(1, 2);
  occ->RemoveElectron(2, 2);
  water->NewConfigurationWithElectronOccupancy("QuadrupleIonisation40", *occ);
  water->AddDecayChannel("QuadrupleIonisation40", new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(1, 2);
  occ->RemoveElectron(3, 2);
  water->NewConfigurationWithElectronOccupancy("QuadrupleIonisation41", *occ);
  water->AddDecayChannel("QuadrupleIonisation41", new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(1, 2);
  occ->RemoveElectron(4, 2);
  water->NewConfigurationWithElectronOccupancy("QuadrupleIonisation42", *occ);
  water->AddDecayChannel("QuadrupleIonisation42", new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(2, 2);
  occ->RemoveElectron(3, 2);
  water->NewConfigurationWithElectronOccupancy("QuadrupleIonisation43", *occ);
  water->AddDecayChannel("QuadrupleIonisation43", new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(2, 2);
  occ->RemoveElectron(4, 2);
  water->NewConfigurationWithElectronOccupancy("QuadrupleIonisation44", *occ);
  water->AddDecayChannel("QuadrupleIonisation44", new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(3, 2);
  occ->RemoveElectron(4, 2);
  water->NewConfigurationWithElectronOccupancy("QuadrupleIonisation45", *occ);
  water->AddDecayChannel("QuadrupleIonisation45", new G4MolecularDissociationChannel(*decCh1));

  delete occ;
}
