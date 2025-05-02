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
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "HistoManager.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::HistoManager()
{
  fHistoId.resize(MaxHisto);
  fExist.resize(MaxHisto);
  fLabel.resize(MaxHisto);
  fTitle.resize(MaxHisto);
  fNbins.resize(MaxHisto);
  fVmin.resize(MaxHisto);
  fVmax.resize(MaxHisto);
  fUnit.resize(MaxHisto);
  fWidth.resize(MaxHisto);
  fascii.resize(MaxHisto);

  // histograms
  for (G4int k=0; k<MaxHisto; k++) {
    fHistoId[k] = 0;
    fExist[k] = false;
    fUnit[k]  = 1.0;
    fWidth[k] = 1.0;
    fascii[k] = false;       
  }

  Book();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::~HistoManager()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::Book()
{
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetDefaultFileType("root");

  // Create or get analysis manager
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetActivation(false);    // enable inactivation of histograms

  // Define the histogramms 
  if(analysisManager->GetNofH1s()==0)
    {
      // create absorber histogramms
      G4int nbins=100;
      G4double vmin=0., vmax=10.;
      for (G4int k=0; k<MaxHisto; k++) 
	{
	  const G4String vunit=(k!=6)?"MeV":"mm";
	  SetHisto( k,  nbins, vmin, vmax, vunit);
	}

     // Creating an 1-dimensional histograms in the root directory of the tree
      for(G4int i=0; i<MaxHisto; i++) 
	{
	  fHistoId[i] = analysisManager->CreateH1(fLabel[i], fTitle[i], fNbins[i], fVmin[i], fVmax[i]);
	  analysisManager->SetH1Activation(fHistoId[i],fExist[i]);
	}

    }

  // Added to catch the SetActivation parameters set througj UI interface
  for (G4int k=0; k<MaxHisto; k++) fExist[k] = analysisManager->GetH1Activation(fHistoId[k]);

   // Check if a filename is set
  if(analysisManager->GetFileName()=="") return;
  
  // Activate the analysisManage ronly if a filename is set
  analysisManager->SetActivation(true);     // enable inactivation of histograms

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::SetHisto(G4int ih,
                 G4int nbins, G4double valmin, G4double valmax, const G4String& unit)
{
  if (ih > MaxHisto) {
    G4cout << "---> warning from HistoManager::SetHisto() : histo " << ih
           << "does not exist" << G4endl;
    return;
  }

  const G4String id[] = { "0", "1", "2", "3", "4", "5", "6" };
  const G4String title[] =
                { "dummy",					//0
                  "continuous energy loss along primary track",	//1
                  "energy from secondaries",			//2
                  "total energy lost by primary track",		//3
		  "energy spectrum of e-+",			//4
		  "energy spectrum of gamma",			//5
		  "step size"					//6
                 };

  G4String titl = fTitle[ih];
  G4double vmin = valmin, vmax = valmax;
  fUnit[ih] = 1.;

  if (unit != "none") {
    titl = title[ih];
    fUnit[ih] = G4UnitDefinition::GetValueOf(unit);
    vmin = valmin/fUnit[ih]; vmax = valmax/fUnit[ih];
  }

  fExist[ih] = true;
  fLabel[ih] = id[ih];
  fTitle[ih] = titl;
  fNbins[ih] = nbins;
  fVmin[ih]  = vmin;
  fVmax[ih]  = vmax;
  fWidth[ih] = (valmax-valmin)/nbins;

  if(titl.substr(0,5)=="dummy") fExist[ih]=false;

  G4cout << "----> SetHisto " << ih << ": " << titl << ";  "
         << nbins << " bins from "
         << vmin << " " << unit << " to " << vmax << " " << unit << " - "<<fExist[ih]<<G4endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::FillHisto(G4int ih, G4double e, G4double weight)
{
  if (ih > MaxHisto) {
    G4cout << "---> warning from HistoManager::FillHisto() : histo " << ih
           << "does not exist; e= " << e << " w= " << weight << G4endl;
    return;
  }

  if(!fExist[ih]) return;
  G4AnalysisManager::Instance()->FillH1( fHistoId[ih],e/fUnit[ih], weight);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


