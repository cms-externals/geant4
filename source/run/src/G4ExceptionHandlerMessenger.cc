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
// G4ExceptionHandlerMessenger implementation
//
// Original author: M.Asai, 1997
// --------------------------------------------------------------------

#include "G4ExceptionHandlerMessenger.hh"

#include "G4ExceptionHandler.hh"
#include "G4Tokenizer.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcommand.hh"
#include "G4UIdirectory.hh"
#include "G4UImanager.hh"
#include "G4UIparameter.hh"

G4ExceptionHandlerMessenger::G4ExceptionHandlerMessenger(G4ExceptionHandler* expH)
  : expHandler(expH)
{
  EHDirectory = new G4UIdirectory("/control/exception/");
  EHDirectory->SetGuidance("Managing number of G4Exception warning messages");

  nwCmd = new G4UIcommand("/control/exception/maxExceptionWarning", this);
  nwCmd->SetGuidance("Set maximum number of G4Exception warning messages to be displayed");
  nwCmd->SetGuidance("for the specified error code.");
  nwCmd->SetGuidance("Number of warnings is counted for each thread individually.");
  nwCmd->SetGuidance("Once the number reaches to the maximum, the warning message of specified");
  nwCmd->SetGuidance("error code won't be printed out any more.");
  nwCmd->SetGuidance("Number is counted through the program execution if more than one run");
  nwCmd->SetGuidance("are executed. Count can be reset with the same command.");
  nwCmd->SetGuidance("If number is set to zero, warning message of the specified error code");
  nwCmd->SetGuidance("won't be displayed at all.");
  nwCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  auto p1 = new G4UIparameter("errorCode", 's', false);
  nwCmd->SetParameter(p1);
  p1 = new G4UIparameter("maxNumber", 'i', false);
  p1->SetParameterRange("maxNumber >= 0");
  nwCmd->SetParameter(p1);

  tnwCmd = new G4UIcmdWithAnInteger("/control/exception/maxTotalExceptionWarning", this);
  tnwCmd->SetGuidance("Set total maximum number of G4Exception warning messages to be displayed.");
  tnwCmd->SetGuidance("Number of warnings is counted for each thread individually.");
  tnwCmd->SetGuidance("Once the number reaches to the maximum, warning messages won't be");
  tnwCmd->SetGuidance("printed out any more.");
  tnwCmd->SetGuidance("Number is counted through the program execution if more than one run");
  tnwCmd->SetGuidance("are executed. Count can be reset with the same command.");
  tnwCmd->SetGuidance("If number is set to zero, warning message won't be displayed at all.");
  tnwCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  tnwCmd->SetParameterName("maxNumber", false);
  tnwCmd->SetRange("maxNumber >= 0");

  wfCmd = new G4UIcommand("/control/exception/writeToFile", this);
  wfCmd->SetGuidance("Write G4Exception warning messages of specified error code to a file.");
  wfCmd->SetGuidance("In the multithreaded mode, each thread writes to a file with thread");
  wfCmd->SetGuidance("id appended to the file name.");
  wfCmd->SetGuidance("To reset, use \"**Screen**\" as the file name.");
  wfCmd->SetGuidance("Please note that if number of warning message exceeds the number");
  wfCmd->SetGuidance("defined by /control/exception/maxExceptionWarning, no additional");
  wfCmd->SetGuidance("message will be written.");
  wfCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  p1 = new G4UIparameter("errorCode", 's', false);
  wfCmd->SetParameter(p1);
  p1 = new G4UIparameter("fileName", 's', true);
  p1->SetDefaultValue("**Screen**");
  wfCmd->SetParameter(p1);

  awfCmd = new G4UIcmdWithAString("/control/exception/writeAllToFile", this);
  awfCmd->SetGuidance("Write all G4Exception warning messages to a file instead of G4cout.");
  awfCmd->SetGuidance("In the multithreaded mode, each thread writes to a file with thread");
  awfCmd->SetGuidance("id appended to the file name.");
  awfCmd->SetGuidance("To reset, use \"**Screen**\" as the file name.");
  awfCmd->SetGuidance("Please note that if number of warning message exceeds the number");
  awfCmd->SetGuidance("defined by /control/exception/maxTotalExceptionWarning, no additional");
  awfCmd->SetGuidance("message will be written.");
  awfCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  awfCmd->SetParameterName("fileName", true);
  awfCmd->SetDefaultValue("**Screen**");
}

// --------------------------------------------------------------------
G4ExceptionHandlerMessenger::~G4ExceptionHandlerMessenger()
{
  delete nwCmd;
  delete tnwCmd;
  delete wfCmd;
  delete awfCmd;
  delete EHDirectory;
}

// --------------------------------------------------------------------
void G4ExceptionHandlerMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if (command == nwCmd)
  {
    G4Tokenizer next(newValue);
    G4String ec = next();
    G4int mc = StoI(next());
    expHandler->SetMaxWarning(ec, mc);
  }
  else if (command == tnwCmd)
  {
    expHandler->SetMaxTotalWarning(tnwCmd->GetNewIntValue(newValue));
  }
  else if (command == wfCmd)
  {
    G4Tokenizer next(newValue);
    G4String ec = next();
    G4String fn = next();
    expHandler->SetWarningFile(ec, fn);
  }
  else if (command == awfCmd)
  {
    expHandler->SetAllWarningFile(newValue);
  }
}

// --------------------------------------------------------------------
G4String G4ExceptionHandlerMessenger::GetCurrentValue(G4UIcommand* /*command*/)
{
  G4String cv;
  return cv;
}
