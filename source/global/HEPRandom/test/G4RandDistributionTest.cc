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
// ----------------------------------------------------------------------
#include "Randomize.hh"
#include "G4ios.hh"

using namespace CLHEP;

HepJamesRandom theJamesEngine;
MixMaxRng theMixMaxEngine;
RanluxEngine theRanluxEngine(19780503,4);
RanluxppEngine theRanluxppEngine(19780503);
RanecuEngine theRanecuEngine;

void init()
{
   char Pause;

   G4cout << G4endl << G4endl;
   G4cout << "---------------------------- Random shooting test -----------------------------" << G4endl;
   G4cout << "                             --------------------                              " << G4endl;
   G4cout << " >>> Random Engines available <<<" << G4endl << G4endl;
   G4cout << "   > HepJamesRandom" << G4endl;
   G4cout << "   > MixMax (default)" << G4endl;
   G4cout << "   > Ranlux" << G4endl;
   G4cout << "   > Ranlux++" << G4endl;
   G4cout << "   > Ranecu" << G4endl << G4endl;
   G4cout << "                   -----  Press <ENTER> to continue  -----";
   if ( (Pause = G4cin.get()) != '\n') exit(0);
   G4cout << G4endl;

}  // end init()


void layout()
{
   float m=3.0;
   const int size=5;
   double vect[size];

   G4cout << " Flat ]0,1[       : " << RandFlat::shoot() << G4endl;
   G4cout << " Flat ]0,5[       : " << RandFlat::shoot(5) << G4endl;
   G4cout << " Flat ]-5,3[      : " << RandFlat::shoot(-5,3) << G4endl;
   G4cout << " Exp (m=1)        : " << RandExponential::shoot() << G4endl;
   G4cout << " Exp (m=3)        : " << RandExponential::shoot(3) << G4endl;
   G4cout << " Gauss (m=1)      : " << RandGauss::shoot() << G4endl;
   G4cout << " Gauss (m=3,v=1)  : " << RandGauss::shoot(3,1) << G4endl;
   G4cout << " Wigner(1,0.2)    : " << RandBreitWigner::shoot(1,0.2) << G4endl;
   G4cout << " Wigner(1,0.2,1)  : " << RandBreitWigner::shoot(1,0.2,1) << G4endl;
   G4cout << " Wigner2(1,0.2)   : " << RandBreitWigner::shootM2(1,0.2) << G4endl;
   G4cout << " Wigner2(1,0.2,1) : " << RandBreitWigner::shootM2(1,0.2,1) << G4endl;
   G4cout << " IntFlat [0,99[   : " << RandFlat::shootInt(99) << G4endl;
   G4cout << " IntFlat [-99,37[ : " << RandFlat::shootInt(-99,37) << G4endl;
   G4cout << " Poisson (m=3.0)  : " << RandPoisson::shoot(m) << G4endl;
   G4cout << G4endl;
   G4cout << " Shooting an array of 5 flat numbers ..." << G4endl << G4endl;
   RandFlat::shootArray(size,vect);
   for ( int i=0; i<size; ++i )
     G4cout << " " << vect[i];
   G4cout << G4endl << G4endl;
}   // end layout() 

void dist_layout()
{
   float m=3.0;
   const int size=5;
   double vect[size];
   char Pause;

   HepJamesRandom aJamesEngine;
   MixMaxRng aMixMaxEngine;
   RanluxEngine aRanluxEngine(19780503,4);
   RanluxppEngine aRanluxppEngine(19780503);
   RanecuEngine aRanecuEngine;

   RandFlat aFlatObj(aJamesEngine);
   RandExponential anExponentialObj(aMixMaxEngine);
   RandGauss aGaussObj(aMixMaxEngine);
   RandBreitWigner aBreitObj(aRanluxEngine);
   RandPoisson aPoissonObj(aRanecuEngine);

   G4cout << "                   -----  Press <ENTER> to continue  -----";
   if ( (Pause = G4cin.get()) != '\n') exit(0);
   G4cout << G4endl;
   G4cout << "-------------------- Shooting test on distribution objects --------------------" << G4endl;
   G4cout << G4endl;
   G4cout << " Flat ]0,1[       : " << aFlatObj.fire() << G4endl;
   G4cout << " Flat ]0,5[       : " << aFlatObj.fire(5) << G4endl;
   G4cout << " Flat ]-5,3[      : " << aFlatObj.fire(-5,3) << G4endl;
   G4cout << " Exp (m=1)        : " << anExponentialObj.fire() << G4endl;
   G4cout << " Exp (m=3)        : " << anExponentialObj.fire(3) << G4endl;
   G4cout << " Gauss (m=1)      : " << aGaussObj.fire() << G4endl;
   G4cout << " Gauss (m=3,v=1)  : " << aGaussObj.fire(3,1) << G4endl;
   G4cout << " Wigner(1,0.2)    : " << aBreitObj.fire(1,0.2) << G4endl;
   G4cout << " Wigner(1,0.2,1)  : " << aBreitObj.fire(1,0.2,1) << G4endl;
   G4cout << " Wigner2(1,0.2)   : " << aBreitObj.fireM2(1,0.2) << G4endl;
   G4cout << " Wigner2(1,0.2,1) : " << aBreitObj.fireM2(1,0.2,1) << G4endl;
   G4cout << " IntFlat [0,99[   : " << aFlatObj.fireInt(99) << G4endl;
   G4cout << " IntFlat [-99,37[ : " << aFlatObj.fireInt(-99,37) << G4endl;
   G4cout << " Poisson (m=3.0)  : " << aPoissonObj.fire(m) << G4endl;
   G4cout << G4endl;
   G4cout << " Shooting an array of 5 flat numbers ..." << G4endl << G4endl;
   aFlatObj.fireArray(size,vect);
   for ( int i=0; i<size; ++i )
     G4cout << " " << vect[i];
   G4cout << G4endl << G4endl;
   G4cout << "                   -----  Press <ENTER> to continue  -----";
   if ( (Pause = G4cin.get()) != '\n') exit(0);
}   // end dist_layout() 

void user_layout()
{
   float m=3.0;
   const int size=5;
   double vect[size];
   char sel;
   HepRandomEngine* anEngine;

   G4cout << G4endl << G4endl;
   G4cout << "-------------------- Shooting test skeeping the generator ---------------------" << G4endl;
   G4cout << G4endl;
   G4cout << " >>> Select a Random Engine <<<" << G4endl << G4endl;
   G4cout << "   1. HepJamesRandom" << G4endl;
   G4cout << "   2. MixMax" << G4endl;
   G4cout << "   3. Ranlux" << G4endl;
   G4cout << "   4. Ranlux++" << G4endl;
   G4cout << "   5. Ranecu" << G4endl << G4endl;
   G4cout << " > ";
   G4cin >> sel;
   while ((sel!='1')&&(sel!='2')&&(sel!='3')&&(sel!='4')&&(sel!='5')) {
     G4cout << G4endl << " >>> Choice not legal !!  [1..5]<<<" << G4endl;
     G4cin >> sel;
   }

   switch (sel) {
     case '1':
       anEngine = &theJamesEngine;
       break;
     case '2':
       anEngine = &theMixMaxEngine;
       break;
     case '3':
       anEngine = &theRanluxEngine;
       break;
     case '4':
       anEngine = &theRanluxppEngine;
       break;
     case '5':
       anEngine = &theRanecuEngine;
       break;

     default:
       anEngine = &theJamesEngine;
       break;
   }
   G4cout << G4endl;
 
   G4cout << " Flat ]0,1[       : " << RandFlat::shoot(anEngine) << G4endl;
   G4cout << " Flat ]0,5[       : " << RandFlat::shoot(anEngine,5) << G4endl;
   G4cout << " Flat ]-5,3[      : " << RandFlat::shoot(anEngine,-5,3) << G4endl;
   G4cout << " Exp (m=1)        : " << RandExponential::shoot(anEngine) << G4endl;
   G4cout << " Exp (m=3)        : " << RandExponential::shoot(anEngine,3) << G4endl;
   G4cout << " Gauss (m=1)      : " << RandGauss::shoot(anEngine) << G4endl;
   G4cout << " Gauss (m=3,v=1)  : " << RandGauss::shoot(anEngine,3,1) << G4endl;
   G4cout << " Wigner(1,0.2)    : " << RandBreitWigner::shoot(anEngine,1,0.2) << G4endl;
   G4cout << " Wigner(1,0.2,1)  : " << RandBreitWigner::shoot(anEngine,1,0.2,1) << G4endl;
   G4cout << " Wigner2(1,0.2)   : " << RandBreitWigner::shootM2(anEngine,1,0.2) << G4endl;
   G4cout << " Wigner2(1,0.2,1) : " << RandBreitWigner::shootM2(anEngine,1,0.2,1) << G4endl;
   G4cout << " IntFlat [0,99[   : " << RandFlat::shootInt(anEngine,99) << G4endl;
   G4cout << " IntFlat [-99,37[ : " << RandFlat::shootInt(anEngine,-99,37) << G4endl;
   G4cout << " Poisson (m=3.0)  : " << RandPoisson::shoot(anEngine,m) << G4endl;
   G4cout << G4endl;
   G4cout << " Shooting an array of 5 flat numbers ..." << G4endl << G4endl;
   RandFlat::shootArray(anEngine,size,vect);
   for ( int i=0; i<size; ++i )
     G4cout << " " << vect[i];
   G4cout << G4endl << G4endl;
}   // end layout() 

void start_test()
{
   char Pause;

   G4cout << "-------------------------  Test on HepJamesRandom  ----------------------------" << G4endl;
   G4cout << G4endl;
   HepRandom::setTheEngine(&theJamesEngine);
   layout();
   G4cout << "                   -----  Press <ENTER> to continue  -----";
   if ( (Pause = G4cin.get()) != '\n') exit(0);
   G4cout << G4endl;
   G4cout << "---------------------------  Test on MixMaxEngine  ------------------------------" << G4endl;
   G4cout << G4endl;
   HepRandom::setTheEngine(&theMixMaxEngine);
   layout();
   G4cout << "                   -----  Press <ENTER> to continue  -----";
   if ( (Pause = G4cin.get()) != '\n') exit(0);
   G4cout << G4endl;
   G4cout << "---------------------  Test on RanluxEngine (luxury 4) ------------------------" << G4endl;
   G4cout << G4endl;
   HepRandom::setTheEngine(&theRanluxEngine);
   layout();
   G4cout << "                   -----  Press <ENTER> to continue  -----";
   if ( (Pause = G4cin.get()) != '\n') exit(0);
   G4cout << G4endl;
   G4cout << "------------------------  Test on Ranlux++Engine ---------------------------" << G4endl;
   G4cout << G4endl;
   HepRandom::setTheEngine(&theRanluxEngine);
   layout();
   G4cout << "                   -----  Press <ENTER> to continue  -----";
   if ( (Pause = G4cin.get()) != '\n') exit(0);
   G4cout << G4endl;
   G4cout << "--------------------------  Test on RanecuEngine ------------------------------" << G4endl;
   G4cout << G4endl;
   HepRandom::setTheEngine(&theRanecuEngine);
   layout();
   dist_layout();
   user_layout();
}  // end start_test()


int main() {

   init();
   start_test();
   
   return 0;
}

