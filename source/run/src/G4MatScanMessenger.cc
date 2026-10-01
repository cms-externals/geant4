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
// G4MatScanMessenger implementation
//
// Author: M.Asai, May 2006
// --------------------------------------------------------------------

#include "G4MatScanMessenger.hh"

#include "G4MaterialScanner.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4Tokenizer.hh"
#include "G4UIcmdWith3Vector.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcommand.hh"
#include "G4UIdirectory.hh"
#include "G4UIparameter.hh"

// --------------------------------------------------------------------
G4MatScanMessenger::G4MatScanMessenger(G4MaterialScanner* p1)
{
  theScanner = p1;
  G4UIparameter* par;
  msDirectory = new G4UIdirectory("/control/matScan/");
  msDirectory->SetGuidance("Material scanner commands.");

  scanCmd = new G4UIcmdWithoutParameter("/control/matScan/scan", this);
  scanCmd->SetGuidance("Start material scanning.");
  scanCmd->SetGuidance("Scanning range should be defined with");
  scanCmd->SetGuidance("/control/matScan/theta and /control/matSca/phi commands.");
  scanCmd->AvailableForStates(G4State_Idle);
  scanCmd->SetToBeBroadcasted(false);

  thetaCmd = new G4UIcommand("/control/matScan/theta", this);
  thetaCmd->SetGuidance("Define theta range.");
  thetaCmd->SetGuidance("Usage : /control/matScan/theta [nbin] [thetaMin] [thetaSpan] [unit]");
  thetaCmd->SetGuidance("Notation of angles :");
  thetaCmd->SetGuidance(" theta --- +Z axis : +90 deg. / X-Y plane : 0 deg. / -Z axis : -90 deg.");
  par = new G4UIparameter("nbin", 'i', false);
  par->SetParameterRange("nbin>0");
  thetaCmd->SetParameter(par);
  par = new G4UIparameter("thetaMin", 'd', false);
  thetaCmd->SetParameter(par);
  par = new G4UIparameter("thetaSpan", 'd', true);
  par->SetParameterRange("thetaSpan>=0.");
  par->SetDefaultValue(0.);
  thetaCmd->SetParameter(par);
  par = new G4UIparameter("unit", 'c', true);
  par->SetDefaultValue("deg");
  par->SetParameterCandidates(thetaCmd->UnitsList(thetaCmd->CategoryOf("deg")));
  thetaCmd->SetParameter(par);
  thetaCmd->SetToBeBroadcasted(false);

  phiCmd = new G4UIcommand("/control/matScan/phi", this);
  phiCmd->SetGuidance("Define phi range.");
  phiCmd->SetGuidance("Usage : /control/matScan/phi [nbin] [phiMin] [phiSpan] [unit]");
  phiCmd->SetGuidance("Notation of angles :");
  phiCmd->SetGuidance(
    " phi   --- +X axis : 0 deg. / +Y axis : 90 deg. / -X axis : 180 "
    "deg. / -Y axis : 270 deg.");
  par = new G4UIparameter("nbin", 'i', false);
  par->SetParameterRange("nbin>0");
  phiCmd->SetParameter(par);
  par = new G4UIparameter("phiMin", 'd', false);
  phiCmd->SetParameter(par);
  par = new G4UIparameter("phiSpan", 'd', true);
  par->SetParameterRange("phiSpan>=0.");
  par->SetDefaultValue(0.);
  phiCmd->SetParameter(par);
  par = new G4UIparameter("unit", 'c', true);
  par->SetDefaultValue("deg");
  par->SetParameterCandidates(phiCmd->UnitsList(phiCmd->CategoryOf("deg")));
  phiCmd->SetParameter(par);
  phiCmd->SetToBeBroadcasted(false);

  singleCmd = new G4UIcommand("/control/matScan/singleMeasure", this);
  singleCmd->SetGuidance("Measure thickness for one particular direction.");
  singleCmd->SetGuidance("Notation of angles :");
  singleCmd->SetGuidance(" theta --- +Z axis : +90 deg. / X-Y plane : 0 deg. / -Z axis : -90 deg.");
  singleCmd->SetGuidance(
    " phi   --- +X axis : 0 deg. / +Y axis : 90 deg. / -X axis : "
    "180 deg. / -Y axis : 270 deg.");
  singleCmd->AvailableForStates(G4State_Idle);
  par = new G4UIparameter("theta", 'd', false);
  singleCmd->SetParameter(par);
  par = new G4UIparameter("phi", 'd', false);
  singleCmd->SetParameter(par);
  par = new G4UIparameter("unit", 'c', true);
  par->SetDefaultValue("deg");
  par->SetParameterCandidates(singleCmd->UnitsList(singleCmd->CategoryOf("deg")));
  singleCmd->SetParameter(par);
  singleCmd->SetToBeBroadcasted(false);

  single2Cmd = new G4UIcmdWith3Vector("/control/matScan/singleTo", this);
  single2Cmd->SetGuidance("Measure thickness for one direction defined by a unit vector.");
  single2Cmd->SetParameterName("X", "Y", "Z", false);
  single2Cmd->SetToBeBroadcasted(false);

  eyePosCmd = new G4UIcmdWith3VectorAndUnit("/control/matScan/eyePosition", this);
  eyePosCmd->SetGuidance("Define the eye position.");
  eyePosCmd->SetParameterName("X", "Y", "Z", true);
  eyePosCmd->SetDefaultValue(G4ThreeVector(0., 0., 0.));
  eyePosCmd->SetDefaultUnit("m");
  eyePosCmd->SetToBeBroadcasted(false);

  regSenseCmd = new G4UIcmdWithABool("/control/matScan/regionSensitive", this);
  regSenseCmd->SetGuidance("Set region sensitivity.");
  regSenseCmd->SetGuidance("This command is automatically set to TRUE");
  regSenseCmd->SetGuidance(" if /control/matScan/region command is issued.");
  regSenseCmd->SetParameterName("senseFlag", true);
  regSenseCmd->SetDefaultValue(false);
  regSenseCmd->SetToBeBroadcasted(false);

  regionCmd = new G4UIcmdWithAString("/control/matScan/region", this);
  regionCmd->SetGuidance("Define region name to be scanned.");
  regionCmd->SetGuidance("/control/matScan/regionSensitive command is automatically");
  regionCmd->SetGuidance("set to TRUE with this command.");
  regionCmd->SetParameterName("region", true);
  regionCmd->SetDefaultValue("DefaultRegionForTheWorld");
  regionCmd->SetToBeBroadcasted(false);

  thinCmd = new G4UIcommand("/control/matScan/ignoreThinMaterials", this);
  thinCmd->SetGuidance("Ignore thin density materials from accumulation of mass length,");
  thinCmd->SetGuidance("radiation length and interaction length.");
  thinCmd->SetGuidance("If material is specified, material whose density is equal to or");
  thinCmd->SetGuidance("thinner than the specified material is not accumulated.");
  thinCmd->SetGuidance("If material is not specified, G4_AIR is used.");
  thinCmd->SetGuidance("Specifying material must be defined in G4NISTMaterialManager");
  thinCmd->SetGuidance("or already defined in the detector construction.");
  par = new G4UIparameter("flag", 'b', true);
  par->SetDefaultValue(true);
  thinCmd->SetParameter(par);
  par = new G4UIparameter("material", 's', true);
  par->SetDefaultValue("G4_AIR");
  thinCmd->SetParameter(par);
  thinCmd->SetToBeBroadcasted(false);

  verboseCmd = new G4UIcmdWithAnInteger("/control/matScan/verbose", this);
  verboseCmd->SetGuidance("Set verbose level of material scan");
  verboseCmd->SetGuidance("0: default, properties integrated over the scan");
  verboseCmd->SetGuidance("1: integrated properties per material");
  verboseCmd->SetGuidance("2: detailed properties per material crossed");
  verboseCmd->SetParameterName("verbose_level", true);
  verboseCmd->SetDefaultValue(0);
  verboseCmd->SetToBeBroadcasted(false);

  defPartCmd = new G4UIcommand("/control/matScan/particleForIntLength", this);
  defPartCmd->SetGuidance("Alternate particle type and its kinetic energy");
  defPartCmd->SetGuidance("Hadronic inelastic cross-section, and hense the hadronic");
  defPartCmd->SetGuidance("interaction length depend on particle type and its kinetic");
  defPartCmd->SetGuidance("energy. By default the material scanner uses 1 GeV/c proton.");
  defPartCmd->SetGuidance("This command alternates the particle type and its kinetic energy.");
  par = new G4UIparameter("pName", 's', false);
  defPartCmd->SetParameter(par);
  par = new G4UIparameter("eKin", 'd', false);
  defPartCmd->SetParameter(par);
  par = new G4UIparameter("unit", 's', true);
  par->SetDefaultUnit("GeV");
  defPartCmd->SetParameter(par);
  defPartCmd->SetToBeBroadcasted(false);

  oFileCmd = new G4UIcmdWithAString("/control/matScan/outputFile", this);
  oFileCmd->SetGuidance("Define the file name of the material scanner output");
  oFileCmd->SetGuidance("If this command is not used, the output is shown on the display");
  oFileCmd->SetGuidance("To reset the setting, use \"**Screen**\" as the file name.");
  oFileCmd->SetGuidance("If more than one scan is made without changing the file name,");
  oFileCmd->SetGuidance("results are appended to the file.");
  oFileCmd->SetParameterName("fName", true);
  oFileCmd->SetDefaultValue("**Screen**");
  oFileCmd->SetToBeBroadcasted(false);
}

// --------------------------------------------------------------------
G4MatScanMessenger::~G4MatScanMessenger()
{
  delete scanCmd;
  delete thetaCmd;
  delete phiCmd;
  delete singleCmd;
  delete single2Cmd;
  delete eyePosCmd;
  delete regSenseCmd;
  delete regionCmd;
  delete thinCmd;
  delete verboseCmd;
  delete defPartCmd;
  delete oFileCmd;
  delete msDirectory;
}

// --------------------------------------------------------------------
G4String G4MatScanMessenger::GetCurrentValue(G4UIcommand* command)
{
  G4String currentValue;
  if (command == thetaCmd)
  {
    currentValue = thetaCmd->ConvertToString(theScanner->GetNTheta());
    currentValue += " ";
    currentValue += thetaCmd->ConvertToString((theScanner->GetThetaMin()) / deg);
    currentValue += " ";
    currentValue += thetaCmd->ConvertToString((theScanner->GetThetaSpan()) / deg);
  }
  else if (command == phiCmd)
  {
    currentValue = phiCmd->ConvertToString(theScanner->GetNPhi());
    currentValue += " ";
    currentValue += phiCmd->ConvertToString((theScanner->GetPhiMin()) / deg);
    currentValue += " ";
    currentValue += phiCmd->ConvertToString((theScanner->GetPhiSpan()) / deg);
  }
  else if (command == eyePosCmd)
  {
    currentValue = eyePosCmd->ConvertToString(theScanner->GetEyePosition(), "m");
  }
  else if (command == regSenseCmd)
  {
    currentValue = regSenseCmd->ConvertToString(theScanner->GetRegionSensitive());
  }
  else if (command == regionCmd)
  {
    currentValue = theScanner->GetRegionName();
  }
  return currentValue;
}

// --------------------------------------------------------------------
void G4MatScanMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if (command == scanCmd)
  {
    theScanner->Scan();
  }
  else if (command == thetaCmd)
  {
    G4Tokenizer next(newValue);
    G4int nbin = StoI(next());
    G4double thetaMin = StoD(next());
    G4double thetaSpan = StoD(next());
    G4String unit = next();
    thetaMin *= thetaCmd->ValueOf(unit);
    thetaSpan *= thetaCmd->ValueOf(unit);
    theScanner->SetNTheta(nbin);
    theScanner->SetThetaMin(thetaMin);
    theScanner->SetThetaSpan(thetaSpan);
  }
  else if (command == phiCmd)
  {
    G4Tokenizer next(newValue);
    G4int nbin = StoI(next());
    G4double phiMin = StoD(next());
    G4double phiSpan = StoD(next());
    G4String unit = next();
    phiMin *= phiCmd->ValueOf(unit);
    phiSpan *= phiCmd->ValueOf(unit);
    theScanner->SetNPhi(nbin);
    theScanner->SetPhiMin(phiMin);
    theScanner->SetPhiSpan(phiSpan);
  }
  else if (command == eyePosCmd)
  {
    theScanner->SetEyePosition(eyePosCmd->GetNew3VectorValue(newValue));
  }
  else if (command == regSenseCmd)
  {
    theScanner->SetRegionSensitive(regSenseCmd->GetNewBoolValue(newValue));
  }
  else if (command == regionCmd)
  {
    if (theScanner->SetRegionName(newValue))
    {
      theScanner->SetRegionSensitive(true);
    }
    else
    {
      G4ExceptionDescription ed;
      ed << "Region <" << newValue << "> is not defined. Command ignored.";
      command->CommandFailed(ed);
    }
  }
  else if (command == verboseCmd)
  {
    theScanner->SetVerbosity(StoI(newValue));
  }
  else if (command == thinCmd)
  {
    G4Tokenizer next(newValue);
    G4bool flg = StoB(next());
    G4String mat = next();
    auto rv = theScanner->SetThinMaterial(flg, mat);
    if (!rv)
    {
      G4ExceptionDescription ed;
      ed << "Material <" << mat << "> is not defined. Command ignored.";
      command->CommandFailed(ed);
    }
  }
  else if (command == singleCmd || command == single2Cmd)
  {
    G4int ntheta = theScanner->GetNTheta();
    G4double thetaMin = theScanner->GetThetaMin();
    G4double thetaSpan = theScanner->GetThetaSpan();
    G4int nphi = theScanner->GetNPhi();
    G4double phiMin = theScanner->GetPhiMin();
    G4double phiSpan = theScanner->GetPhiSpan();

    G4double theta = 0.;
    G4double phi = 0.;
    if (command == singleCmd)
    {
      G4Tokenizer next(newValue);
      theta = StoD(next());
      phi = StoD(next());
      G4String unit = next();
      theta *= singleCmd->ValueOf(unit);
      phi *= singleCmd->ValueOf(unit);
    }
    else if (command == single2Cmd)
    {
      G4ThreeVector v = single2Cmd->GetNew3VectorValue(newValue);
      theta = 90. * deg - v.theta();
      phi = v.phi();
    }
    theScanner->SetNTheta(1);
    theScanner->SetThetaMin(theta);
    theScanner->SetThetaSpan(0.);
    theScanner->SetNPhi(1);
    theScanner->SetPhiMin(phi);
    theScanner->SetPhiSpan(0.);
    theScanner->Scan();

    theScanner->SetNTheta(ntheta);
    theScanner->SetThetaMin(thetaMin);
    theScanner->SetThetaSpan(thetaSpan);
    theScanner->SetNPhi(nphi);
    theScanner->SetPhiMin(phiMin);
    theScanner->SetPhiSpan(phiSpan);
  }
  else if (command == defPartCmd)
  {
    G4Tokenizer next(newValue);
    G4String pName = next();
    auto* ptcl = G4ParticleTable::GetParticleTable()->FindParticle(pName);
    if (ptcl == nullptr)
    {
      G4ExceptionDescription ed;
      ed << "Particle <" << pName << "> is not found. Command ignored.";
      command->CommandFailed(ed);
    }
    else
    {
      G4double eKin = StoD(next());
      G4String unit = next();
      eKin *= defPartCmd->ValueOf(unit);
      theScanner->SetDetDefaultParticle(ptcl, eKin);
    }
  }
  else if (command == oFileCmd)
  {
    theScanner->SetOutFile(newValue);
  }
}
