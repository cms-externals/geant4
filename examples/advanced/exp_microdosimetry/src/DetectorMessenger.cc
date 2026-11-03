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
// Authors: Susanna Guatelli, Francesco Romano and Giuliana Milluzzo
// susanna@uow.edu.au, francesco.romano@ct.infn.it, giuliana.milluzzo@ct.infn.it
//
// Created by Jacopo Magini: j.magini@surrey.ac.uk

#include "DetectorMessenger.hh"

#include "DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAnInteger.hh"

DetectorMessenger::DetectorMessenger(AnalysisManager* analysis_manager) : geometryHasChanged(false)
{
  changeTheGeometryDir = new G4UIdirectory("/geometrySetup/");
  changeTheGeometryDir->SetGuidance("Geometry setup");

  changeTheDetectorCmd = new G4UIcmdWithAString("/geometrySetup/selectDetector", this);
  changeTheDetectorCmd->SetGuidance("Select the type of detector you wish to use");
  changeTheDetectorCmd->SetParameterName("Material", false);
  changeTheDetectorCmd->AvailableForStates(G4State_PreInit);

  changeDetectorDimensionDir = new G4UIdirectory("/geometrySetup/detectorDimension/");
  changeDetectorDimensionDir->SetGuidance("Modify the dimensions of the detector");

  changeDetectorSizeWidthCmd =
    new G4UIcmdWithADoubleAndUnit("/geometrySetup/detectorDimension/setWidth", this);
  changeDetectorSizeWidthCmd->SetGuidance("Set the width of the detector");
  changeDetectorSizeWidthCmd->SetParameterName("Width", false);
  changeDetectorSizeWidthCmd->SetRange("Width >= 0.1 && Width <= 150.");
  changeDetectorSizeWidthCmd->SetUnitCategory("Length");
  changeDetectorSizeWidthCmd->SetDefaultUnit("um");
  changeDetectorSizeWidthCmd->AvailableForStates(G4State_PreInit);

  changeDetectorSizeThicknessCmd =
    new G4UIcmdWithADoubleAndUnit("/geometrySetup/detectorDimension/setThickness", this);
  changeDetectorSizeThicknessCmd->SetGuidance("Set the thickness of the detector");
  changeDetectorSizeThicknessCmd->SetParameterName("Thickness", false);
  changeDetectorSizeThicknessCmd->SetRange("Thickness >= 0.1 && Thickness <= 50.");
  changeDetectorSizeThicknessCmd->SetUnitCategory("Length");
  changeDetectorSizeThicknessCmd->SetDefaultUnit("um");
  changeDetectorSizeThicknessCmd->AvailableForStates(G4State_PreInit);

  changeSecondStageDimensionDir =
    new G4UIdirectory("/geometrySetup/detectorDimension/secondStage/");
  changeSecondStageDimensionDir->SetGuidance(
    "Modify the dimensions of the E stage for the telescope detector");

  changeSecondStageSizeWidthCmd =
    new G4UIcmdWithADoubleAndUnit("/geometrySetup/detectorDimension/secondStage/setWidth", this);
  changeSecondStageSizeWidthCmd->SetGuidance("Set the width of the E-stage for telescope detector");
  changeSecondStageSizeWidthCmd->SetParameterName("Width", false);
  changeSecondStageSizeWidthCmd->SetRange("Width >= 10 && Width <= 1000.");
  changeSecondStageSizeWidthCmd->SetUnitCategory("Length");
  changeSecondStageSizeWidthCmd->SetDefaultUnit("um");
  changeSecondStageSizeWidthCmd->AvailableForStates(G4State_PreInit);

  changeSecondStageSizeThicknessCmd = new G4UIcmdWithADoubleAndUnit(
    "/geometrySetup/detectorDimension/secondStage/setThickness", this);
  changeSecondStageSizeThicknessCmd->SetGuidance(
    "Set the thickness of the E-stage for telescope detector");
  changeSecondStageSizeThicknessCmd->SetParameterName("Thickness", false);
  changeSecondStageSizeThicknessCmd->SetRange("Thickness >= 10 && Thickness <= 1000.");
  changeSecondStageSizeThicknessCmd->SetUnitCategory("Length");
  changeSecondStageSizeThicknessCmd->SetDefaultUnit("um");
  changeSecondStageSizeThicknessCmd->AvailableForStates(G4State_PreInit);

  siliconArrayDir =
    new G4UIdirectory("/geometrySetup/siliconArray/");
  siliconArrayDir->SetGuidance("Configure the silicon detector array");

  setNumSV_XCmd = new G4UIcmdWithAnInteger("/geometrySetup/siliconArray/setNumSV_X", this);
  setNumSV_XCmd->SetGuidance("Set the number of silicon detectors along X");
  setNumSV_XCmd->SetParameterName("numSV_X", false);
  setNumSV_XCmd->SetRange("numSV_X >= 1");
  setNumSV_XCmd->AvailableForStates(G4State_PreInit);

  setNumSV_YCmd = new G4UIcmdWithAnInteger("/geometrySetup/siliconArray/setNumSV_Y", this);
  setNumSV_YCmd->SetGuidance("Set the number of silicon detectors along Y");
  setNumSV_YCmd->SetParameterName("numSV_Y", false);
  setNumSV_YCmd->SetRange("numSV_Y >= 1");
  setNumSV_YCmd->AvailableForStates(G4State_PreInit);
  
  setNumPMMABlocksXCmd = new G4UIcmdWithAnInteger("/geometrySetup/siliconArray/setNumPMMABlocksX", this);
  setNumPMMABlocksXCmd->SetGuidance("Set the number of PMMA blocks along X");
  setNumPMMABlocksXCmd->SetParameterName("numPMMABlocksX", false);
  setNumPMMABlocksXCmd->SetRange("numPMMABlocksX >= 1");
  setNumPMMABlocksXCmd->AvailableForStates(G4State_PreInit);
  
  setNumPMMABlocksYCmd = new G4UIcmdWithAnInteger("/geometrySetup/siliconArray/setNumPMMABlocksY", this);
  setNumPMMABlocksYCmd->SetGuidance("Set the number of PMMA blocks along Y");
  setNumPMMABlocksYCmd->SetParameterName("numPMMABlocksY", false);
  setNumPMMABlocksYCmd->SetRange("numPMMABlocksY >= 1");
  setNumPMMABlocksYCmd->AvailableForStates(G4State_PreInit);

  setSV_XspacingCmd = new G4UIcmdWithADoubleAndUnit("/geometrySetup/siliconArray/setSV_Xspacing", this);
  setSV_XspacingCmd->SetGuidance("Set the centre-to-centre detector pitch along X");
  setSV_XspacingCmd->SetParameterName("SV_Xspacing", false);
  setSV_XspacingCmd->SetRange("SV_Xspacing >= 0.");
  setSV_XspacingCmd->SetUnitCategory("Length");
  setSV_XspacingCmd->SetDefaultUnit("um");
  setSV_XspacingCmd->AvailableForStates(G4State_PreInit);

  setSV_YspacingCmd = new G4UIcmdWithADoubleAndUnit("/geometrySetup/siliconArray/setSV_Yspacing", this);
  setSV_YspacingCmd->SetGuidance("Set the centre-to-centre detector pitch along Y");
  setSV_YspacingCmd->SetParameterName("SV_Yspacing", false);
  setSV_YspacingCmd->SetRange("SV_Yspacing >= 0.");
  setSV_YspacingCmd->SetUnitCategory("Length");
  setSV_YspacingCmd->SetDefaultUnit("um");
  setSV_YspacingCmd->AvailableForStates(G4State_PreInit);
  
  setBlockGapXCmd = new G4UIcmdWithADoubleAndUnit("/geometrySetup/siliconArray/setBlockGapX", this);
  setBlockGapXCmd->SetGuidance("Set the edge-to-edge gap between PMMA blocks along X");
  setBlockGapXCmd->SetParameterName("blockGapX", false);
  setBlockGapXCmd->SetRange("blockGapX >= 0.");
  setBlockGapXCmd->SetUnitCategory("Length");
  setBlockGapXCmd->SetDefaultUnit("um");
  setBlockGapXCmd->AvailableForStates(G4State_PreInit);

  setBlockGapYCmd =new G4UIcmdWithADoubleAndUnit("/geometrySetup/siliconArray/setBlockGapY", this);
  setBlockGapYCmd->SetGuidance("Set the edge-to-edge gap between PMMA blocks along Y");
  setBlockGapYCmd->SetParameterName("blockGapY", false);
  setBlockGapYCmd->SetRange("blockGapY >= 0.");
  setBlockGapYCmd->SetUnitCategory("Length");
  setBlockGapYCmd->SetDefaultUnit("um");
  setBlockGapYCmd->AvailableForStates(G4State_PreInit);
  
  enableWaterPhantomCmd = new G4UIcmdWithABool("/geometrySetup/enableWaterPhantom", this);
  enableWaterPhantomCmd->SetGuidance("If true, the detector is placed inside a water phantom");
  enableWaterPhantomCmd->SetParameterName("Phantom", false);
  enableWaterPhantomCmd->AvailableForStates(G4State_PreInit);

  changeDetectorPositionDir = new G4UIdirectory("/geometrySetup/detectorPosition/");
  changeDetectorPositionDir->SetGuidance("Modify the placement of the detector");

  changeDetectorPositionDepthCmd =
    new G4UIcmdWithADoubleAndUnit("/geometrySetup/detectorPosition/setDepth", this);
  changeDetectorPositionDepthCmd->SetGuidance("Set the detector depth inside the water phantom");
  changeDetectorPositionDepthCmd->SetParameterName("Depth", false);
  changeDetectorPositionDepthCmd->SetRange("Depth >= 1. && Depth <= 250.");
  changeDetectorPositionDepthCmd->SetUnitCategory("Length");
  changeDetectorPositionDepthCmd->SetDefaultUnit("mm");
  changeDetectorPositionDepthCmd->AvailableForStates(G4State_PreInit);

  applyChangesToGeometryCmd = new G4UIcmdWithoutParameter("/geometrySetup/applyChanges", this);
  applyChangesToGeometryCmd->SetGuidance("Apply selected changes to the geometry");
  applyChangesToGeometryCmd->AvailableForStates(G4State_PreInit);

  // default values
  detectorType = "Diamond";
  detectorWidth = 30. * um;
  detectorThickness = 10. * um;
  phantomEnabled = false;
  detectorDepth = 0 * mm;
  secondStageWidth = 800. * um;
  secondStageThickness = 500. * um;

  numSV_X = 3;
  numSV_Y = 3;
  numPMMABlocksX = 1;
  numPMMABlocksY = 1;
  SV_Xspacing = detectorWidth;
  SV_Yspacing = detectorWidth;
  blockGapX = 0.;
  blockGapY = 0.;

  analysis = analysis_manager;
}

DetectorMessenger::~DetectorMessenger()
{
  delete changeTheDetectorCmd;
  delete changeDetectorSizeWidthCmd;
  delete changeDetectorSizeThicknessCmd;
  delete changeSecondStageSizeWidthCmd;
  delete changeSecondStageSizeThicknessCmd;
  delete enableWaterPhantomCmd;
  delete changeDetectorPositionDepthCmd;

  delete applyChangesToGeometryCmd;

  delete changeSecondStageDimensionDir;
  delete changeDetectorDimensionDir;
  delete changeDetectorPositionDir;
  delete setNumSV_XCmd;
  delete setNumSV_YCmd;
  delete setNumPMMABlocksXCmd;
  delete setNumPMMABlocksYCmd;
  delete setSV_XspacingCmd;
  delete setSV_YspacingCmd;
  delete setBlockGapXCmd;
  delete setBlockGapYCmd;
  delete siliconArrayDir;

  delete changeTheGeometryDir;
}

void DetectorMessenger::SetNewValue(G4UIcommand* command, G4String commandContent)
{
  if (command == changeTheDetectorCmd)
  {
    if (commandContent == "Diamond" || commandContent == "MicroDiamond"
        || commandContent == "Silicon"|| commandContent == "SiliconArray"
	|| commandContent == "SiliconBridge" || commandContent == "SiCDetector"
	|| commandContent == "DiamondTelescope" || commandContent == "TEPC")
    {
      detectorType = commandContent;
      geometryHasChanged = true;

      G4cout << "Detector type changed to " << commandContent << G4endl;
      G4cout << "Run /geometrySetup/applyChanges to apply" << G4endl;
    }

    else
    {
      G4cout << "Unknown detector type: " << commandContent << ". Geometry not changed." << G4endl;
    }
  }

  else if (command == changeDetectorSizeWidthCmd)
  {
    detectorWidth = G4UIcmdWithADoubleAndUnit::GetNewDoubleValue(commandContent);
    geometryHasChanged = true;

    G4cout << "Detector width changed to " << commandContent << G4endl;
    G4cout << "Run /geometrySetup/applyChanges to apply" << G4endl;
  }

  else if (command == changeDetectorSizeThicknessCmd)
  {
    detectorThickness = G4UIcmdWithADoubleAndUnit::GetNewDoubleValue(commandContent);
    geometryHasChanged = true;

    G4cout << "Detector thickness changed to " << commandContent << G4endl;
    G4cout << "Run /geometrySetup/applyChanges to apply" << G4endl;
  }

  else if (command == changeSecondStageSizeWidthCmd)
  {
    secondStageWidth = G4UIcmdWithADoubleAndUnit::GetNewDoubleValue(commandContent);
    geometryHasChanged = true;

    G4cout << "Telescope E-stage width changed to " << commandContent << G4endl;
    if (detectorType != "DiamondTelescope")
    {
      G4cout << "However this setting only takes effect when using the telescope detector.";
      G4cout << "Select it with '/geometrySetup/selectDetector DiamondTelescope' or this command "
                "will be ignored"
             << G4endl;
    }
    else
      G4cout << "Run /geometrySetup/applyChanges to apply" << G4endl;
  }

  else if (command == changeSecondStageSizeThicknessCmd)
  {
    secondStageThickness = G4UIcmdWithADoubleAndUnit::GetNewDoubleValue(commandContent);
    geometryHasChanged = true;

    G4cout << "Telescope E-stage thickness changed to " << commandContent << G4endl;
    if (detectorType != "DiamondTelescope")
    {
      G4cout << "However this setting only takes effect when using the telescope detector.";
      G4cout << "Select it with '/geometrySetup/selectDetector DiamondTelescope' or this command "
                "will be ignored"
             << G4endl;
    }
    else
      G4cout << "Run /geometrySetup/applyChanges to apply" << G4endl;
  }

  else if (command == setNumSV_XCmd)
  {
    numSV_X = G4UIcmdWithAnInteger::GetNewIntValue(commandContent);
    geometryHasChanged = true;
		
    G4cout << "Number of SV along X changed to " << commandContent << G4endl;
    G4cout << "Run /geometrySetup/applyChanges to apply" << G4endl;
  }
	
 else if (command == setNumSV_YCmd)
  {
    numSV_Y = G4UIcmdWithAnInteger::GetNewIntValue(commandContent);
    geometryHasChanged = true;
		
    G4cout << "Number of SV along Y changed to " << commandContent << G4endl;
    G4cout << "Run /geometrySetup/applyChanges to apply" << G4endl;
  }
	
 else if (command == setNumPMMABlocksXCmd)
  {
    numPMMABlocksX = G4UIcmdWithAnInteger::GetNewIntValue(commandContent);
    geometryHasChanged = true;
		
    G4cout << "Number of PMMA blocks along X changed to " << commandContent << G4endl;
    G4cout << "Run /geometrySetup/applyChanges to apply" << G4endl;
  }
 else if (command == setNumPMMABlocksYCmd)
  {
    numPMMABlocksY = G4UIcmdWithAnInteger::GetNewIntValue(commandContent);
    geometryHasChanged = true;
		
    G4cout << "Number of PMMA blocks along Y changed to " << commandContent << G4endl;
    G4cout << "Run /geometrySetup/applyChanges to apply" << G4endl;
   }
	
 else if (command == setSV_XspacingCmd)
  {
    SV_Xspacing = setSV_XspacingCmd->GetNewDoubleValue(commandContent);
    geometryHasChanged = true;
		
    G4cout << "SV distance along X changed to " << commandContent << G4endl;
    G4cout << "Run /geometrySetup/applyChanges to apply" << G4endl;
  }
	
 else if (command == setSV_YspacingCmd)
  {
    SV_Yspacing =setSV_YspacingCmd->GetNewDoubleValue(commandContent);
    geometryHasChanged = true;
		
    G4cout << "SV distance along Y changed to " << commandContent << G4endl;
    G4cout << "Run /geometrySetup/applyChanges to apply" << G4endl;
  }
	
 else if (command == setBlockGapXCmd)
  {
    blockGapX = setBlockGapXCmd->GetNewDoubleValue(commandContent);
    geometryHasChanged = true;
		
    G4cout << "Block gap along X changed to " << commandContent << G4endl;
    G4cout << "Run /geometrySetup/applyChanges to apply" << G4endl;
  }
	
 else if (command == setBlockGapYCmd)
  {
    blockGapY =setBlockGapYCmd->GetNewDoubleValue(commandContent);
    geometryHasChanged = true;
		
    G4cout << "Block gap along Y changed to " << commandContent << G4endl;
    G4cout << "Run /geometrySetup/applyChanges to apply" << G4endl;
  }
  
  else if (command == enableWaterPhantomCmd)
  {
    phantomEnabled = G4UIcmdWithABool::GetNewBoolValue(commandContent);
    geometryHasChanged = true;

    if (phantomEnabled == true)
      G4cout << "Water phantom enabled" << G4endl;
    else
      G4cout << "Water phantom disabled" << G4endl;
    G4cout << "Run /geometrySetup/applyChanges to apply" << G4endl;
  }

  else if (command == changeDetectorPositionDepthCmd)
  {
    detectorDepth = G4UIcmdWithADoubleAndUnit::GetNewDoubleValue(commandContent);
    geometryHasChanged = true;

    G4cout << "Detector depth in water changed to " << commandContent << G4endl;
    if (phantomEnabled == false)
    {
      G4cout << "However the water phantom is not enabled.";
      G4cout
        << "Enable it with '/geometrySetup/enableWaterPhantom true' or this command will be ignored"
        << G4endl;
    }
    else
      G4cout << "Run /geometrySetup/applyChanges to apply" << G4endl;
  }

  else if (command == applyChangesToGeometryCmd)
  {
    if (geometryHasChanged == true)
    {
      G4RunManager* runManager = G4RunManager::GetRunManager();

      if (detectorType == "SiliconArray" && ((numSV_X > 1 && SV_Xspacing < detectorWidth) ||
     								(numSV_Y > 1 && SV_Yspacing < detectorWidth)))
      {
    	 G4Exception("DetectorMessenger::SetNewValue",
        	      "InvalidSiliconArrayConfiguration",
        	      FatalException,
        	      "SiliconArray configuration causes overlapping silicon detectors: "
        	      "the centre-to-centre pitch must be at least the detector diameter.");
    		      return;
      }
      
      DetectorConstruction* newDetector = new DetectorConstruction(analysis, this);

      runManager->SetUserInitialization(newDetector);
      runManager->GeometryHasBeenModified();

      geometryHasChanged = false;

      G4cout << "All pending changes to geometry have been applied" << G4endl;
    }

    else
    {
      G4cout << "No pending changes to geometry. This command will be ignored" << G4endl;
    }
  }
}
