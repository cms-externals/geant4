// A. Paillet
// email: paillet@lp2ib.in2p3.fr

void txtToRoot()
{
        std::ifstream inputFile("example.txt");
        if (!inputFile.is_open()) {
        cout << "Error: Could not open the text file." << endl;
        }

        // Read and discard the header line
        std::string header;
        std::getline(inputFile, header);

        // Read the file line by line
        std::string line;

        int partID;
        double x, y, z, px, py, pz, e, w;

        // Units : x[mm] y[mm] z[mm] px[MeV/c] py[MeV/c] pz[MeV/c] e[MeV]

        // Create the root TTree
        TFile *hfile = TFile::Open("example.root","RECREATE");
        TTree *tree = new TTree("PSF","test");
        tree->Branch("partID", &partID,"partID/I");
        tree->Branch("x", &x,"x/D");
        tree->Branch("y", &y,"y/D");
        tree->Branch("z", &z,"z/D");
        tree->Branch("px", &px,"px/D");
        tree->Branch("py", &py,"py/D");
        tree->Branch("pz", &pz,"pz/D");
        tree->Branch("e", &e,"e/D");
        tree->Branch("w", &w,"w/D");

        while (std::getline(inputFile, line)) {
            std::istringstream iss(line);

            // Read the values from the line
            if (iss >> partID >> x >> y >> z >> px >> py >> pz >> e >> w) {
                // Store the values in the respective vectors
              tree->Fill(); 
            } 
            else {
                cout << "Error: Could not read all values from line: " << line << endl;
            }
        }
        tree->Print();
        tree->Write();

        // Close the file
        inputFile.close();
}