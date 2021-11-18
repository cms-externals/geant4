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
//
// 
// ------------------------------------------------------------

#include "G4MaterialPropertyVector.hh"

void LoopUntilPressEnter();

int main ()
{
        const G4int NUMENTRIES = 32;
	G4double aPhotonMomentum, anAbsorptionCoefficient;

        G4double PPCKOV[NUMENTRIES] =
	          { 4.010E-9, 4.144E-9, // see (**)
		    2.038E-9, 2.072E-9, 2.107E-9, 2.143E-9, 2.181E-9,
                    2.220E-9, 2.260E-9, 2.302E-9, 2.346E-9, 2.391E-9,
                    2.438E-9, 2.486E-9, 2.537E-9, 2.590E-9, 2.645E-9,
                    2.702E-9, 2.763E-9, 2.825E-9, 2.891E-9, 2.960E-9,
                    3.032E-9, 3.108E-9, 3.188E-9, 3.271E-9, 3.360E-9,
                    3.453E-9, 3.552E-9, 3.656E-9, 3.767E-9, 3.884E-9,
                    };

        G4double ABSCO[NUMENTRIES] =
	         {  1750.0, 1450.0, // see (**)
		    344.8,  408.2,  632.9,  917.4, 1234.6, 1388.9,
                    1515.2, 1724.1, 1886.8, 2000.0, 2631.6, 3571.4,
                    4545.5, 4761.9, 5263.2, 5263.2, 5555.6, 5263.2,
                    5263.2, 4761.9, 4545.5, 4166.7, 3703.7, 3333.3,
                    3000.0, 2850.0, 2700.0, 2450.0, 2200.0, 1950.0,
                    };

	// (**) Were at the end of the lists: I put them at the beginning to
	// test reordering.


	// Test Vector creation
	// --------------------
	G4cout << "Test Vector creation" << G4endl;
	G4cout << "--------------------" << G4endl << G4endl;

	G4MaterialPropertyVector absco(PPCKOV, ABSCO, NUMENTRIES);
	absco.DumpVector();
	LoopUntilPressEnter();
		
	// Test GetProperty
	// ----------------	
	G4cout << "Test GetProperty" << G4endl;
	G4cout << "----------------" << G4endl << G4endl; 
        G4cout << "Input a photon momentum value (MeV) for which you wish" << G4endl;
        G4cout << "to find the absorption coefficient:  " << G4endl;

        G4cin >> aPhotonMomentum;
	G4cout << G4endl;

	aPhotonMomentum *= aPhotonMomentum*MeV;

	anAbsorptionCoefficient = absco.GetProperty(aPhotonMomentum);

	G4cout << "The absorption coefficient is:  " 
	     << anAbsorptionCoefficient  
	     << G4endl;
        LoopUntilPressEnter();

	// Test AddElement
	// ---------------
	G4cout << "Test AddElement" << G4endl;
	G4cout << "---------------" << G4endl << G4endl;
	G4cout << "Add the entry just created" << G4endl;
	LoopUntilPressEnter();

	absco.AddElement(aPhotonMomentum, anAbsorptionCoefficient); 
	absco.DumpVector();

	// Test RemoveElement
	// ------------------
	G4cout << "Test RemoveElement" << G4endl;
        G4cout << "------------------" << G4endl << G4endl;
        G4cout << "Input a photon momentum value for which you wish" << G4endl;
        G4cout << "to remove the OPVEntry:  " << G4endl;

        G4cin >> aPhotonMomentum;
	G4cout << G4endl;

	aPhotonMomentum *= aPhotonMomentum*MeV;

	absco.RemoveElement(aPhotonMomentum);
	absco.DumpVector(); 

	// Test the iterator
	// -----------------
	G4cout << "Test the iterator" << G4endl;
	G4cout << "-----------------" << G4endl << G4endl;
        LoopUntilPressEnter();

	absco.ResetIterator();

	while (++absco)
	{
		G4cout << absco.GetPhotonMomentum()/MeV << "MeV \t"
		     << absco.GetProperty() << G4endl;
	}
	G4cout << "\n\n <END OF TEST>" << G4endl;

	return EXIT_SUCCESS;
}

// LoopUntilPressEnter
// -------------------
//
void LoopUntilPressEnter()
{
        char ch;
        G4cout << "Press <Enter> to continue ... ";
        while ( G4cin.get(ch) )
        {
                if (ch == '\n') break;
        }
        G4cout << G4endl;
}
