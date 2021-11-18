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
// -------------------------------------------------------------
//
// 
// -------------------------------------------------------------
//      GEANT4 hTest
//
//      History: based on object model of
//      2nd December 1995, G.Cosmo
//      ---------- hTestDetectorMessenger -------
//              
//  Modified: 05.04.01 Vladimir Ivanchenko new design of hTest 
// 
// -------------------------------------------------------------

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "hTestDetectorMessenger.hh"
#include "hTestDetectorConstruction.hh"
#include "hTestHisto.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

hTestDetectorMessenger::hTestDetectorMessenger(hTestDetectorConstruction* h):
  hDet(h)
{ 
  hTestdetDir = new G4UIdirectory("/hTest/");
  hTestdetDir->SetGuidance("General hTest commands");
  hTestdetDir1= new G4UIdirectory("/hTest/physics/");
  hTestdetDir1->SetGuidance("hTest commands to define physics");
  hTestdetDir2= new G4UIdirectory("/hTest/gun/");
  hTestdetDir2->SetGuidance("hTest commands to define gun");
  if(hDet->GetVerbose() > 0) {
    G4cout << "hTestDetectorMessenger: Is constructed" << G4endl;
  }
      
  AbsMaterCmd = new G4UIcmdWithAString("/hTest/AbsorberMaterial",this);
  AbsMaterCmd->SetGuidance("Select Material of the Absorber.");
  AbsMaterCmd->SetParameterName("AbsoberMaterial",false);
  AbsMaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  WorldMaterCmd = new G4UIcmdWithAString("/hTest/WorldMaterial",this);
  WorldMaterCmd->SetGuidance("Select Material of the World.");
  WorldMaterCmd->SetParameterName("WorldMaterial",false);
  WorldMaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  NumOfAbsCmd = new G4UIcmdWithAnInteger("/hTest/NumberOfAbsorbers",this);
  NumOfAbsCmd->SetGuidance("Set number of absorbers");
  NumOfAbsCmd->SetParameterName("Nabs",false);
  NumOfAbsCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  AbsThickCmd = new G4UIcmdWithADoubleAndUnit("/hTest/AbsorberThick",this);
  AbsThickCmd->SetGuidance("Set Thickness of the Absorber");
  AbsThickCmd->SetParameterName("SizeZ",false);  
  AbsThickCmd->SetRange("SizeZ>0.");
  AbsThickCmd->SetUnitCategory("Length");  
  AbsThickCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  AbsGapCmd = new G4UIcmdWithADoubleAndUnit("/hTest/AbsorberGap",this);
  AbsGapCmd->SetGuidance("Set gap between absorbers");
  AbsGapCmd->SetParameterName("SizeZ",false);  
  AbsGapCmd->SetRange("SizeZ>0.");
  AbsGapCmd->SetUnitCategory("Length");  
  AbsGapCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  AbsSizYZCmd = new G4UIcmdWithADoubleAndUnit("/hTest/AbsorberXY",this);
  AbsSizYZCmd->SetGuidance("Set sizeXY of the Absorber");
  AbsSizYZCmd->SetParameterName("SizeYZ",false);
  AbsSizYZCmd->SetRange("SizeYZ>0.");
  AbsSizYZCmd->SetUnitCategory("Length");
  AbsSizYZCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
    
  WorldXCmd = new G4UIcmdWithADoubleAndUnit("/hTest/WorldZ",this);
  WorldXCmd->SetGuidance("Set Z size of the World");
  WorldXCmd->SetParameterName("WSizeX",false);
  WorldXCmd->SetRange("WSizeX>0.");
  WorldXCmd->SetUnitCategory("Length");
  WorldXCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
    
  UpdateCmd = new G4UIcmdWithoutParameter("/hTest/update",this);
  UpdateCmd->SetGuidance("Update calorimeter geometry.");
  UpdateCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
  UpdateCmd->SetGuidance("if you changed geometrical value(s).");
  UpdateCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
      
  XMagFieldCmd = new G4UIcmdWithADoubleAndUnit("/hTest/FieldX",this);  
  XMagFieldCmd->SetGuidance("Define magnetic field along X");
  XMagFieldCmd->SetGuidance("Magnetic field will be in X direction.");
  XMagFieldCmd->SetParameterName("Bx",false);
  XMagFieldCmd->SetUnitCategory("Magnetic flux density");
  XMagFieldCmd->AvailableForStates(G4State_PreInit,G4State_Idle);  

  YMagFieldCmd = new G4UIcmdWithADoubleAndUnit("/hTest/FieldY",this);  
  YMagFieldCmd->SetGuidance("Define magnetic field along Y");
  YMagFieldCmd->SetGuidance("Magnetic field will be in Y direction.");
  YMagFieldCmd->SetParameterName("By",false);
  YMagFieldCmd->SetUnitCategory("Magnetic flux density");
  YMagFieldCmd->AvailableForStates(G4State_PreInit,G4State_Idle);  

  ZMagFieldCmd = new G4UIcmdWithADoubleAndUnit("/hTest/FieldZ",this);  
  ZMagFieldCmd->SetGuidance("Define magnetic field along Z");
  ZMagFieldCmd->SetGuidance("Magnetic field will be in Z direction.");
  ZMagFieldCmd->SetParameterName("Bz",false);
  ZMagFieldCmd->SetUnitCategory("Magnetic flux density");
  ZMagFieldCmd->AvailableForStates(G4State_PreInit,G4State_Idle);  

  HistoCmd = new G4UIcmdWithAString("/hTest/HistoName",this);
  HistoCmd->SetGuidance("Set the name of the histo file");
  HistoCmd->SetParameterName("histo",false);
  HistoCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  NumOfEvt = new G4UIcmdWithAnInteger("/hTest/NumberOfEvents",this);
  NumOfEvt->SetGuidance("Set number of event to be simulated");
  NumOfEvt->SetParameterName("Nevt",false);
  NumOfEvt->AvailableForStates(G4State_PreInit,G4State_Idle);

  verbCmd = new G4UIcmdWithAnInteger("/hTest/verbose",this);
  verbCmd->SetGuidance("Set verbose for hTest");
  verbCmd->SetParameterName("verb",false);
  verbCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  intCmd = new G4UIcmdWithAnInteger("/hTest/numberAbsToSave",this);
  intCmd->SetGuidance("Set number of absorbers for which "); 
  intCmd->SetGuidance("the energy is saved to tuple");
  intCmd->SetParameterName("numberAbsToSave",false);
  intCmd->AvailableForStates(G4State_PreInit);

  nhistCmd = new G4UIcmdWithAnInteger("/hTest/HistoNumber",this);
  nhistCmd->SetGuidance("Set number of histograms to fill"); 
  nhistCmd->SetParameterName("HistoNumber",false);
  nhistCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  ntupCmd = new G4UIcmdWithABool("/hTest/ntuple",this);
  ntupCmd->SetGuidance("Set number ntuple to fill"); 
  ntupCmd->SetParameterName("ntuple",false);
  ntupCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  nDebugSCmd = new G4UIcmdWithAnInteger("/hTest/nFirstEventToDebug",this);
  nDebugSCmd->SetGuidance("Set number of the first event to debug"); 
  nDebugSCmd->SetParameterName("nFirstEventToDebug",false);
  nDebugSCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  nDebugECmd = new G4UIcmdWithAnInteger("/hTest/nLastEventToDebug",this);
  nDebugECmd->SetGuidance("Set number of the last event to debug"); 
  nDebugECmd->SetParameterName("nLastEventToDebug",false);
  nDebugECmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  DeltaECmd = new G4UIcmdWithADoubleAndUnit("/hTest/maxDeltaEnergy",this);  
  DeltaECmd->SetGuidance("Define scale of delta-Energy histogram");
  DeltaECmd->SetParameterName("DeltaE",false);
  DeltaECmd->SetUnitCategory("Energy");
  DeltaECmd->AvailableForStates(G4State_PreInit,G4State_Idle);  
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

hTestDetectorMessenger::~hTestDetectorMessenger()
{
  delete NumOfAbsCmd; 
  delete AbsMaterCmd; 
  delete AbsThickCmd; 
  delete AbsGapCmd; 
  delete AbsSizYZCmd;  
  delete WorldMaterCmd;
  delete WorldXCmd;
  delete UpdateCmd;
  delete XMagFieldCmd;
  delete YMagFieldCmd;
  delete ZMagFieldCmd;
  delete HistoCmd;
  delete NumOfEvt;
  delete verbCmd;
  delete intCmd;
  delete nhistCmd;
  delete nDebugSCmd;
  delete nDebugECmd;
  delete hTestdetDir;
  delete hTestdetDir1;
  delete hTestdetDir2;
  delete DeltaECmd;
  delete ntupCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestDetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  if(hDet->GetVerbose() > 1) {
    G4cout << "hTestDetectorMessenger: new value = " << newValue << G4endl;
  }

  if( command == NumOfAbsCmd )
   { hDet->SetNumberOfAbsorbers(NumOfAbsCmd->GetNewIntValue(newValue));
     (hTestHisto::GetPointer())->SetNumberOfAbsorbers(NumOfAbsCmd->GetNewIntValue(newValue));
   }

  if( command == AbsMaterCmd )
   { hDet->SetAbsorberMaterial(newValue);}
   
  if( command == AbsThickCmd )
   { hDet->SetAbsorberThickness(AbsThickCmd->GetNewDoubleValue(newValue));
    (hTestHisto::GetPointer())->SetAbsorberThickness(AbsThickCmd->GetNewDoubleValue(newValue));
   }

  if( command == AbsGapCmd )
   { hDet->SetGap(AbsGapCmd->GetNewDoubleValue(newValue));
    (hTestHisto::GetPointer())->SetGap(AbsGapCmd->GetNewDoubleValue(newValue));
   }

  if( command == WorldMaterCmd )
   { hDet->SetWorldMaterial(newValue);}
   
  if( command == AbsSizYZCmd )
   { hDet->SetAbsorberSizeXY(AbsSizYZCmd->GetNewDoubleValue(newValue));}
      
  if( command == WorldXCmd )
   { hDet->SetWorldSizeZ(WorldXCmd->GetNewDoubleValue(newValue));}
      
  if( command == UpdateCmd )
   { hDet->UpdateGeometry(); }

  if( command == XMagFieldCmd )
   { hDet->SetMagField(XMagFieldCmd->GetNewDoubleValue(newValue),1);}

  if( command == YMagFieldCmd )
   { hDet->SetMagField(YMagFieldCmd->GetNewDoubleValue(newValue),2);}

  if( command == ZMagFieldCmd )
   { hDet->SetMagField(ZMagFieldCmd->GetNewDoubleValue(newValue),3);}

  if( command == HistoCmd ) 
   { (hTestHisto::GetPointer())->SetHistoName(newValue);}

  if( command == ntupCmd ) 
   { (hTestHisto::GetPointer())->SetNtuple(ntupCmd->GetNewBoolValue(newValue));}

  if( command == NumOfEvt )
   { hDet->SetNumberOfEvents(NumOfAbsCmd->GetNewIntValue(newValue));}

  if( command == verbCmd ){ 
     G4int ver = verbCmd->GetNewIntValue(newValue);
     hDet->SetVerbose(ver);
     (hTestHisto::GetPointer())->SetVerbose(ver);
   }

  if( command == intCmd )
   { hDet->SetNumAbsorbersSaved(intCmd->GetNewIntValue(newValue));}

  if( command == nhistCmd )
   { (hTestHisto::GetPointer())->SetHistoNumber(nhistCmd->GetNewIntValue(newValue));}

  if( command == nDebugSCmd )
   { hDet->SetFirstEventToDebug(nDebugSCmd->GetNewIntValue(newValue));}

  if( command == nDebugECmd )
   { hDet->SetLastEventToDebug(nDebugECmd->GetNewIntValue(newValue));}

  if( command == DeltaECmd )
   { (hTestHisto::GetPointer())
      ->SetMaxEnergy(DeltaECmd->GetNewDoubleValue(newValue));}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
