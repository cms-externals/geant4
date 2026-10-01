// A. Paillet
// email: paillet@lp2ib.in2p3.fr

#include <iostream>
#include <fstream>
#include <TFile.h>
#include <TTree.h>

void ConvertRootToText(const char* rootFileName, const char* textFileName) {
    // Open the ROOT file
    TFile* rootFile = new TFile(rootFileName, "READ");
    if (!rootFile || rootFile->IsZombie()) {
        std::cerr << "Error: Could not open ROOT file: " << rootFileName << std::endl;
        return;
    }

    // Access the TTree
    TTree* tree = (TTree*)rootFile->Get("PSF"); // Replace "tree" with your TTree name if different
    if (!tree) {
        std::cerr << "Error: Could not find TTree in the ROOT file." << std::endl;
        rootFile->Close();
        return;
    }

    // Define variables for the branches
    Int_t partID;
    Double_t x, y, z, px, py, pz, e, w;

    // Set branch addresses
    tree->SetBranchAddress("partID", &partID);
    tree->SetBranchAddress("x", &x);
    tree->SetBranchAddress("y", &y);
    tree->SetBranchAddress("z", &z);
    tree->SetBranchAddress("px", &px);
    tree->SetBranchAddress("py", &py);
    tree->SetBranchAddress("pz", &pz);
    tree->SetBranchAddress("e", &e);
    tree->SetBranchAddress("w", &w);

    // Open the text file for writing
    std::ofstream textFile(textFileName);
    if (!textFile.is_open()) {
        std::cerr << "Error: Could not open text file for writing: " << textFileName << std::endl;
        rootFile->Close();
        return;
    }

    // Write header to the text file
    textFile << "partID x[mm] y[mm] z[mm] px[MeV/c] py[MeV/c] pz[MeV/c] e[MeV] w\n"; 
    // partID is PDG code
    // e is total energy (kinetic + mass)

    // Loop over all entries in the TTree
    Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);

        // Write data to the text file
        textFile << partID << " "
                 << x << " " << y << " " << z << " "
                 << px << " " << py << " " << pz << " "
                 << e << " " << w << "\n";
    }

    // Close files
    textFile.close();
    rootFile->Close();

    std::cout << "Conversion completed. Data written to: " << textFileName << std::endl;
}

void rootToTxt() {
    // Example usage
    ConvertRootToText("example.root", "example.txt");
    
}