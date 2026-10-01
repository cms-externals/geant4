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
//Created: Ji Won Choi, Adrien Paillet

#include "PSFReader.hh"
#include <filesystem>      // for internal helper GetFileList
#include "G4IonTable.hh"     // for GetIonTable
#include "G4MTRunManager.hh" // for AbortRun
#include "AnalysisManager.hh"
#include <random>          // for G4UniformRand()

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
namespace fs = std::filesystem;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

static std::vector<std::string> GetFileList(const std::string &pattern) {
  std::vector<std::string> files;

  for (const auto &entry: fs::directory_iterator(fs::path(pattern).parent_path())) {
    if (entry.path().string().find(pattern) != std::string::npos) {
      files.push_back(entry.path().string());
    }
  }
  return files;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PSFReader::InitializePhaseSpace() {
  if (g_isLoaded) {return;}

  if (fpattern.empty()) {
    G4cout << "[INFO] /psf/pattern was not set in the macro." << G4endl;
    return;
  }

  if (fpattern.find(".root", fpattern.size() - 5) != std::string::npos) //5 is size of .root
  {
    G4cout << "[INFO] Loading ROOT PSF (MASTER ONLY)" << G4endl;
    ReadRootFiles(fpattern);
    g_isLoaded = true;

    G4cout << "[INFO] Total loaded particles: "
        << g_particles.size() << G4endl;
  } else if (fpattern.find(".txt", fpattern.size() - 4) != std::string::npos) //4 is size of .txt
  {
    G4cout << "[INFO] Loading TEXT PSF (MASTER ONLY)" << G4endl;
    ReadTextFiles(fpattern);
    g_isLoaded = true;

    G4cout << "[INFO] Total loaded particles: "
        << g_particles.size() << G4endl;
  } 
  else if (fpattern.find(".csv", fpattern.size() - 4) != std::string::npos) //4 is size of .csv
  {
    G4cout << "[INFO] Loading CSV PSF (MASTER ONLY)" << G4endl;
    Read_GRAS_CSV_Files(fpattern);
    g_isLoaded = true;

    G4cout << "[INFO] Total loaded particles: "
        << g_particles.size() << G4endl;
  } else {
    G4cerr << "[ERROR] PSF format not supported" << G4endl;
    G4MTRunManager::GetMasterRunManager()->AbortRun();
  }

  // You may call a position-modifying function here.

  // You may use one of the index-modifying functions here.
  // A messenger is used to control fPSFsamplingMode value.

  if (fPSFsamplingMode == "Sampling") {
    Shuffle();
  } else if (fPSFsamplingMode == "Sampling_wRplc") {
    OversampleWithReplacement(
      G4RunManager::GetRunManager()->GetNumberOfEventsToBeProcessed()
    );
  }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PSFReader::ReadRootFiles(const G4String &pattern) {
  const auto files = GetFileList(pattern);

  G4cout << "[INFO] Found " << files.size() << " ROOT files" << G4endl;

  auto analysisReader = G4RootAnalysisReader::Instance();
  analysisReader->SetVerboseLevel(0);

  PhaseSpaceParticle PSP;

  for (auto &fname: files) {
    G4cout << "Reading file : " << fname << G4endl;

    const G4int ntupleId = analysisReader->GetNtuple(fPSFNTupleName.c_str(), fname.c_str());

    if (ntupleId >= 0) {
      G4int partID; // PDG code
      G4double x, y, z, px, py, pz, e, w;

      analysisReader->SetNtupleIColumn("partID", partID);
      analysisReader->SetNtupleDColumn("x", x);
      analysisReader->SetNtupleDColumn("y", y);
      analysisReader->SetNtupleDColumn("z", z);
      analysisReader->SetNtupleDColumn("px", px);
      analysisReader->SetNtupleDColumn("py", py);
      analysisReader->SetNtupleDColumn("pz", pz);
      analysisReader->SetNtupleDColumn("e", e);
      analysisReader->SetNtupleDColumn("w", w);

      while (analysisReader->GetNtupleRow(ntupleId)) {
        PSP.pdg = partID;
        PSP.pos.set(x, y, z);
        PSP.dir.set(px, py, pz);
        PSP.eng = e;
        PSP.weight = w;

        g_particles.push_back(PSP);
      }
    }
  }

  G4cout << "End of root file(s) reading" << G4endl;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PSFReader::ReadTextFiles(const G4String &pattern) {
  // the text file header must looks like :
  // partID x y z px py pz e w
  // the data are arranged line by line in this way :
  // 22 -1000 -430.714 -345.084 0.0287535 -0.0266418 0.0607417 0.0722918 1

  // Parser adapted from Mistral AI
  // Open the text file
  auto files = GetFileList(pattern);

  G4cout << "[INFO] Found " << files.size() << " TEXT files" << G4endl;
  PhaseSpaceParticle PSP;

  for (auto &fname: files) {
    std::ifstream inputFile(fname);
    if (!inputFile.is_open()) {
      G4cerr << "Error: Could not open the text file." << G4endl;
    }

    // Read and discard the header line
    std::string header;
    std::getline(inputFile, header);

    // Read the file line by line
    std::string line;

    G4int partID;
    G4double x, y, z, px, py, pz, e, w;

    while (std::getline(inputFile, line)) {
      std::istringstream iss(line);

      // Read the values from the line
      if (iss >> partID >> x >> y >> z >> px >> py >> pz >> e >> w) {
        // Store the values in the respective vectors

        PSP.pdg = partID;
        PSP.pos.set(x, y, z);
        PSP.dir.set(px, py, pz);
        PSP.eng = e;
        PSP.weight = w;

        g_particles.push_back(PSP);
      } else {
        G4cerr << "Error: Could not read all values from line: " << line << G4endl;
      }
    }

    // Close the file
    inputFile.close();
  }
  // Print the size of each vector to verify
  G4cout << "Number of entries read: " << g_particles.size() << G4endl;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PSFReader::Read_GRAS_CSV_Files(const G4String &pattern)
{
  // Adapted from  : PrimaryGeneratorSourceGRASCSV

  // Open the file(s)
  auto files = GetFileList(pattern);

  G4cout << "[INFO] Found " << files.size() << " CSV files" << G4endl;
  PhaseSpaceParticle PSP;

  for (auto &fname: files) {
    std::ifstream inputFile(fname);
    if (!inputFile.is_open()) {
      G4cerr << "Error: Could not open the CSV file." << G4endl;
    }

    // Prepare variables
    std::string line;
    G4int partID;
    G4double x, y, z, px, py, pz, e, w;  // e is kinetic energy here (not used in generate primary)

    // Not used yet
    // G4int eID,tID,A,Z,Q,instanceID;
    // G4double Ex,T,xl,yl,zl,pxl,pyl,pzl,D;

    std::list<G4String> tempList;
    G4String tempString = "";

    G4bool headerRead = false;
    G4bool EndOfFile = false;

    // Parsing loop
    while (std::getline(inputFile, line) && !EndOfFile) {
      std::istringstream iss(line);

      // Header reading
      if (line.find
          ("TWO-STAGE PARTICLE PHASE SPACE OUTPUT FILE: PARTICLE PHASE-SPACE INFORMATION")
          != std::string::npos)
      {
        for (G4int i = 0; i < 18; i++) {
          std::getline(inputFile, line);
          // Print the header information
          G4cout << line << G4endl;
        }
        headerRead = true;      
        std::getline(inputFile, line); // Go to next line after header : first particle
      }

      // PSF reading 
      if (headerRead == true )
      {
        // Break down the line
        std::stringstream ss(line);
        while(ss>>tempString)
           tempList.push_front(tempString);                       
           
        if (line.find("End") != std::string::npos) 
        { 
          G4cout << "End of file reached" << G4endl;
          EndOfFile = true ;
        }         
        else
        {  
          // Convert and store all the values in a variable (see header for description)
          tempList.pop_back(); // event ID
          tempList.pop_back(); // track ID
          partID = std::stod(tempList.back()); tempList.pop_back(); // particle PDG code
          tempList.pop_back(); // atomic number
          tempList.pop_back(); // nucleon number
          tempList.pop_back(); // charge
          tempList.pop_back(); // Nuclear excitation energy
          w = std::stod(tempList.back()); tempList.pop_back(); // weight
          e = std::stod(tempList.back()); tempList.pop_back(); // kinetic energy
          tempList.pop_back(); // time
          x = std::stod(tempList.back()); tempList.pop_back(); // global position
          y = std::stod(tempList.back()); tempList.pop_back();
          z = std::stod(tempList.back()); tempList.pop_back();
          tempList.pop_back(); // local position
          tempList.pop_back();
          tempList.pop_back();
          px = std::stod(tempList.back()); tempList.pop_back(); // global momentum ("Direction cosine (global)")
          py = std::stod(tempList.back()); tempList.pop_back();
          pz = std::stod(tempList.back()); tempList.pop_back();
          tempList.pop_back(); // local momentum
          tempList.pop_back();
          tempList.pop_back();
          tempList.pop_back();  // "Distance from last volume/scatter"
          tempList.pop_back(); // "Instance ID of PV"

        // Store the needed values in the PSP vector
        PSP.pdg = partID;
        PSP.pos.set(x, y, z);
        PSP.dir.set(px, py, pz);
        PSP.eng = e;
        PSP.weight = w;

        g_particles.push_back(PSP);
        } 
      }    
    }
    // Close the file
    inputFile.close();
  }
  // Print the size of each vector to verify
  G4cout << "Number of entries read: " << g_particles.size() << G4endl;  
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// this function should be called by the primary generator to get the particle data from the PSF
G4PrimaryVertex *PSFReader::GetNextPrimary() {
  // PSF sampling taken from ICRP145 example.
  size_t id = g_nextIndex.fetch_add(1, std::memory_order_relaxed);
  // the function will loop over the PSF if
  // the number of events to be processed is higher than what is in the PSF
  const size_t nbOfParticles = g_particles.size();
  // nbOfParticles > 0 is checked in PrimaryGeneratorAction
  id = id % nbOfParticles;

  if (id == nbOfParticles - 1) {
    // inform the user that end of PSF is reached
    if (const auto eventNumber = static_cast<std::size_t>(G4RunManager::GetRunManager()->
      GetNumberOfEventsToBeProcessed());
      eventNumber > nbOfParticles) {
      G4cerr << "[WARNING] PSF end reached at entry " << id << ". Run continues with entry 0." <<
         G4endl;
    }
  }

  const PhaseSpaceParticle PSP = g_particles[id];
  const auto particleTable = G4ParticleTable::GetParticleTable();
  auto particle = particleTable->FindParticle(PSP.pdg);
  if (!particle) {
    const G4int Z = (PSP.pdg / 10000) % 1000;
    const G4int A = (PSP.pdg / 10) % 1000;
    particle = particleTable->GetIonTable()->GetIon(Z, A, 0.);
  }

  // Vertex at position t=0
  auto *vertex = new G4PrimaryVertex(PSP.pos, 0.0);

  // Primary particle with momentum vector (px, py, pz)
  auto *primary = new G4PrimaryParticle(particle,
                                        PSP.dir.x(),
                                        PSP.dir.y(),
                                        PSP.dir.z());
  primary->SetWeight(PSP.weight);

  vertex->SetPrimary(primary);
  return vertex;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PSFReader::mapOnASphere(const G4double offset_x, const G4double offset_y, const G4double offset_z,
                             const G4double radius_x, const G4double radius_y, const G4double radius_z) {
  // units passed to function must be mm
  // maps particles positions to a sphere
  // parameters could be passed via a PSF messenger
  // (to be implemented if needed)

  for (PhaseSpaceParticle &PSP: g_particles) {
    // offset_x/y/z is the centroid of the points
    // this will center the resulting sphere on 0
    const G4double x = PSP.pos.x() + offset_x;
    const G4double y = PSP.pos.y() + offset_y;
    const G4double z = PSP.pos.z() + offset_z;

    const G4double vecLength = sqrt(x * x + y * y + z * z);

    G4double scale_x = radius_x / vecLength;
    G4double scale_y = radius_y / vecLength;
    G4double scale_z = radius_z / vecLength;

    PSP.pos.set(x * scale_x, y * scale_y, z * scale_z);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// The functions below are examples of how to reorganize the data
// before the run begins
// random sampling of a vector with MT can be challenging,
// so a simple workaround is to modify the indexation before sampling

// Auxiliary function for Shuffle
template<class data>
static void swapVectorWithIndex(std::vector<data> v, std::vector<G4int> newIdx) {
  data temp;

  for (size_t i = 0; i < newIdx.size(); i++) {
    temp = v[i];
    v[i] = v[newIdx[i]];
    v[newIdx[i]] = temp;
  }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// auxiliary function for SampleWithReplacement
template<class data>
static void sampleVectorWithIndex(std::vector<data> v, std::vector<G4int> newIdx) {
  std::vector<data> copy_v(v);

  for (size_t i = 0; i < newIdx.size(); i++) {
    v[i] = copy_v[newIdx[i]];
  }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// auxiliary function for OversampleWithReplacement
template<class data>
static void rewriteVectorWithIndex(std::vector<data> v, std::vector<G4int> newIdx, size_t n) {
  std::vector<data> copy_v(v);
  v.clear();

  for (size_t i = 0; i < n; i++) {
    v.push_back(copy_v[newIdx[i]]);
  }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// index-modifying functions
void PSFReader::SampleWithReplacement() const {
  // Generate a new indexation for all the vectors
  // The same particle can be sampled several times
  size_t i = 0;
  const size_t n = g_particles.size();
  std::vector<G4int> newIdx;

  for (i = 0; i < n; i++) {
    newIdx.push_back(static_cast<size_t>(G4UniformRand() * (n - 1)));
  }

  sampleVectorWithIndex(g_particles, newIdx);

  G4cerr << "[INFO] data have been sampled with replacement" << G4endl;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PSFReader::Shuffle() const {
  // Generate a new indexation for all the vectors
  // all particles are sampled one time in a
  // random order
  size_t i = 0;
  const size_t n = g_particles.size();
  std::vector<G4int> newIdx;

  for (i = 0; i < n; i++) {
    newIdx.push_back(static_cast<size_t>(G4UniformRand() * (n - 1)));
  }

  swapVectorWithIndex(g_particles, newIdx);

  G4cerr << "[INFO] data have been shuffled" << G4endl;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PSFReader::OversampleWithReplacement(const size_t &nOfEvents) const
// n is intended to be the number of events to be processed, accessible using :
// G4RunManager::GetRunManager()->GetNumberOfEventsToBeProcessed()
{
  // If we want more particles than what is in the PSF,
  // and also that the particles are randomly sampled
  // N.B. : the same particle can be sampled several times
  // Warning : increases the size of the vectors in the RAM
  // Note : this is not actually an oversampling
  const size_t nOfParticles = g_particles.size();

  if (nOfParticles >= nOfEvents) { SampleWithReplacement(); } else {
    size_t i = 0;
    std::vector<G4int> newIdx;

    // For each particle to be shot, we
    // associate a random index of the PSF vector
    for (i = 0; i < nOfEvents; i++) {
      newIdx.push_back(static_cast<size_t>(G4UniformRand() * (nOfParticles - 1)));
    }

    rewriteVectorWithIndex(g_particles, newIdx, nOfEvents);

    G4cerr << "[INFO] data have been oversampled with replacement" << G4endl;
  }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
