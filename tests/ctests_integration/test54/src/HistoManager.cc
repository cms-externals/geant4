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
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

#include <sstream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::HistoManager()
{
  fHistoId.assign(MaxHisto+1,-1);
  fLabel.resize(MaxHisto+1);
  fTitle.resize(MaxHisto+1);
  fNbins.resize(MaxHisto+1);
  fVmin.resize(MaxHisto+1);
  fVmax.resize(MaxHisto+1);        
  fUnit.resize(MaxHisto+1,1.);
  fExist.resize(MaxHisto+1,false);
  fMaster = false;

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

  // Create or get analysis manager
  analysisManager->SetDefaultFileType("root");
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetActivation(false);    // enable inactivation of histograms

  // Define the histogramms 
  if(analysisManager->GetNofH1s()==0)
    {
      // create absorber histogramms
      G4int nbins=100;
      G4double vmin=0., vmax=10.;
      const G4String vunit=G4String("none");
      for (G4int k=0; k<MaxHisto; k++) SetHisto( k,  nbins, vmin, vmax, vunit);
      
      // create selected histograms
      for (G4int k=0; k<MaxHisto; k++) {	
	bool bActive=fExist[k];
	fHistoId[k] = analysisManager->CreateH1(fLabel[k], fTitle[k],
						fNbins[k], fVmin[k], fVmax[k]);
	analysisManager->SetH1Activation(fHistoId[k], bActive);
      }
    }

  // Added to catch the SetActivation parameters set througj UI interface
  for (G4int k=0; k<MaxHisto; k++) {
    fExist[k] = analysisManager->GetH1Activation(fHistoId[k]);
  }
  // Check if a filename is set
  if(analysisManager->GetFileName()=="") return;
  
  // Activate the analysisManage ronly if a filename is set
  analysisManager->SetActivation(true);     // enable inactivation of histograms
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::SetHisto(G4int ih,
            G4int nbins, G4double valmin, G4double valmax, const G4String& unit)
{
  if (ih > MaxHisto) {
    G4cout << "---> warning from HistoManager::SetHisto() : histo " << ih
           << "does not exist" << G4endl;
    return;
  }

  const G4String id[] = { "0", "1", "2", "3", "4", "5", "6", "7", "8", "9",
                         "10","11","12","13","14","15","16","17","18","19",
			 "20","21","22","23","24","25","26","27","28","29",
			 "30","31","32","33","34","35","36","37","38","39",
			 "40","41","42","43","44","45","46","47","48","49" 
			};
			
  const G4String title[] =
                { "dummy",						//0
                  "energy deposit in absorber",				//1
                  "energy of charged secondaries at creation",		//2
                  "energy of gammas at creation (std::log10(ekin/MeV))",//3
                  "x_vertex of charged secondaries (all)",		//4
                  "x_vertex of charged secondaries (not absorbed)",	//5
		  "dummy","dummy","dummy","dummy",			//6-9
		  "(transmit, charged) : kinetic energy at exit",	//10
		  "(transmit, charged) : ener fluence: dE(MeV)/dOmega",	//11
		  "(transmit, charged) : space angle: dN/dOmega",	//12
		  "(transmit, charged) : projected angle at exit",	//13
		  "(transmit, charged) : projected position at exit",	//14
		  "(transmit, charged) : radius at exit",		//15
		  "dummy","dummy","dummy","dummy",			//16-19
		  "(transmit, neutral) : kinetic energy at exit",	//20
		  "(transmit, neutral) : ener fluence: dE(MeV)/dOmega",	//21
		  "(transmit, neutral) : space angle: dN/dOmega",	//22
		  "(transmit, neutral) : projected angle at exit",	//23
		  "dummy","dummy","dummy","dummy","dummy","dummy",	//24-29
		  "(reflect , charged) : kinetic energy at exit",	//30
		  "(reflect , charged) : ener fluence: dE(MeV)/dOmega",	//31
		  "(reflect , charged) : space angle: dN/dOmega",	//32
		  "(reflect , charged) : projected angle at exit",	//33
		  "dummy","dummy","dummy","dummy","dummy","dummy",	//34-39
		  "(reflect , neutral) : kinetic energy at exit",	//40
		  "(reflect , neutral) : ener fluence: dE(MeV)/dOmega",	//41
		  "(reflect , neutral) : space angle: dN/dOmega",	//42
		  "(reflect , neutral) : projected angle at exit",	//43
		  "dummy","dummy","dummy","dummy","dummy","dummy"	//44-49
                 };

  G4String titl = title[ih];
  G4double vmin = valmin, vmax = valmax;
  fUnit[ih] = 1.;

  if (ih == 3) { 
    valmin=std::max(valmin,0.01);  // to avoid log(0.)
    vmin = std::log10(valmin/MeV); vmax = std::log10(valmax/MeV);
  }
  else if (unit != "none") {
    titl = title[ih] + " (" + unit + ")";
    fUnit[ih] = G4UnitDefinition::GetValueOf(unit);
    vmin = valmin/fUnit[ih]; vmax = valmax/fUnit[ih];
  }

  fExist[ih] = true;
  fLabel[ih] = id[ih];
  fTitle[ih] = titl;
  fNbins[ih] = nbins;
  fVmin[ih]  = vmin;
  fVmax[ih]  = vmax;

  if(titl.substr(0,5)=="dummy")fExist[ih]=false;

  if(fMaster) {
    G4cout << "----> SetHisto " << ih << ": " << titl << ";  "
	   << nbins << " bins from "
	   << vmin << " " << unit << " to " << vmax << " unit : " << unit << 
      " -> active "<<fExist[ih]<<G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double HistoManager::GetHistoBinWidth(G4int id) const   
{ 
  if(fNbins[id]<1) return -1.;
  return (fVmax[id]-fVmin[id])/(G4double)fNbins[id]; 
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

void HistoManager::Scale(G4int ih, G4double fac)
{
 if (ih > MaxHisto) {
    G4cout << "---> warning from HistoManager::Scale() : histo " << ih
           << "does not exist.  (fac = " << fac << ")" << G4endl;
    return;
  }

  if(!fExist[ih]) return;
  G4AnalysisManager::Instance()->GetH1(fHistoId[ih])->scale(fac);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
/*
#include <fstream>

void HistoManager::saveAscii()
{

  if (!ascii[0]) return;
 
  G4String name = fileName[0] + ".ascii";
  std::ofstream File(name, std::ios::out);
  File.setf( std::ios::scientific, std::ios::floatfield );
 
  //write selected histograms
  for (G4int ih=0; ih<MaxHisto; ih++) {
    if (exist[ih] && ascii[ih]) {
      File << "\n  1D histogram " << ih << ": " << Title[ih] 
           << "\n \n \t     X \t\t     Y" << G4endl;
     
      for (G4int iBin=0; iBin<Nbins[ih]; iBin++) {
         File << "  " << iBin << "\t" 
              << 0.5*(histo[ih]->axis().binLowerEdge(iBin) +
	              histo[ih]->axis().binUpperEdge(iBin)) << "\t"	      
	      << histo[ih]->binHeight(iBin) 
	      << G4endl;
      } 
    }
  }
}
*/
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


