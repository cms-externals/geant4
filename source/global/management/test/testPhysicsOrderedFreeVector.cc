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
#include "globals.hh"
#include "G4PhysicsOrderedFreeVector.hh"

void LoopUntilPressEnter();

int main ()
{
        const G4int NUMENTRIES = 32;
	G4double anEnergy, aValue;

        G4double PPCKOV[NUMENTRIES] =
                  { 2.038E-9, 2.072E-9, 2.107E-9, 2.143E-9, 2.181E-9,
                    2.220E-9, 2.260E-9, 2.302E-9, 2.346E-9, 2.391E-9,
                    2.438E-9, 2.486E-9, 2.537E-9, 2.590E-9, 2.645E-9,
                    2.702E-9, 2.763E-9, 2.825E-9, 2.891E-9, 2.960E-9,
                    3.032E-9, 3.108E-9, 3.188E-9, 3.271E-9, 3.360E-9,
                    3.453E-9, 3.552E-9, 3.656E-9, 3.767E-9, 3.884E-9,
                    4.010E-9, 4.144E-9 };

        G4double RINDEX[NUMENTRIES] =
                 {  1.33, 1.33, 1.33, 1.33, 1.33, 1.33, 1.33,
                    1.33, 1.33, 1.34, 1.34, 1.34, 1.34, 1.34,
                    1.34, 1.34, 1.34, 1.34, 1.34, 1.34, 1.34,
                    1.34, 1.34, 1.35, 1.35, 1.35, 1.35, 1.35,
                    1.35, 1.35, 1.35, 1.35 };

	// Test Vector creation
	// --------------------
	G4cout << "Test Vector creation" << G4endl;
	G4cout << "--------------------" << G4endl << G4endl; 

	G4PhysicsOrderedFreeVector aVector(PPCKOV, RINDEX, NUMENTRIES);
	aVector.DumpValues();
        G4cout << "Press <Enter> to continue ... " << G4endl;
	LoopUntilPressEnter();

	// Test GetEnergy
	// --------------
	G4cout << "Test GetEnergy" << G4endl;
	G4cout << "--------------" << G4endl;
	G4cout << "Input a value within the vector range for which you" << G4endl;
	G4cout << "wish to find the corresponding energy:  " << G4endl;
	G4cin >> aValue;

	anEnergy = aVector.GetEnergy(aValue);
	G4cout << "The corresponding energy is " << anEnergy << G4endl;

	// Test GetMaxValue 
	// ----------------	
	G4cout << "Test GetMaxValue" << G4endl;
	G4cout << "----------------" << G4endl << G4endl
               << "Press <Enter> to continue ... " << G4endl;
        LoopUntilPressEnter();

	aValue = aVector.GetMaxValue();

	G4cout << "The Max Value is:  " << aValue << G4endl
               << "Press <Enter> to continue ... " << G4endl;
	LoopUntilPressEnter();

	// Test GetMinValue 
	// ----------------
	G4cout << "Test GetMinValue" << G4endl;
	G4cout << "----------------" << G4endl << G4endl;
 
	aValue = aVector.GetMinValue();

	G4cout << "The Max Value is:  " << aValue << G4endl
               << "Press <Enter> to continue ... " << G4endl;
	LoopUntilPressEnter();

	// Test GetMaxLowEdgeEnergy 
	// ------------------------
	G4cout << "Test GetMaxLowEdgeEnergy" << G4endl;
	G4cout << "------------------------" << G4endl << G4endl;
 
	anEnergy = aVector.GetMaxLowEdgeEnergy();

	G4cout << "The Max Value is:  " << anEnergy << G4endl
               << "Press <Enter> to continue ... " << G4endl;
	LoopUntilPressEnter();

	// Test GetMinLowEdgeEnergy 
	// ------------------------
	G4cout << "Test GetMinLowEdgeEnergy" << G4endl;
	G4cout << "------------------------" << G4endl << G4endl;
 
	anEnergy = aVector.GetMinLowEdgeEnergy();

	G4cout << "The Max Value is:  " << anEnergy << G4endl;

        return EXIT_SUCCESS;
}

// LoopUntilPressEnter
// -------------------
//
void LoopUntilPressEnter()
{
        char ch;
        while ( G4cin.get(ch) )
        {
                if (ch == '\n') break;
        }
        G4cout << G4endl;
}
