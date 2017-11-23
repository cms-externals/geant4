#ifndef G4NRESP71M03_HH
#define G4NRESP71M03_HH

//#include "G4PhysicalConstants.hh"
//#include "G4SystemOfUnits.hh"
//#include "Randomize.hh"

#include "G4ReactionProduct.hh"
#include "globals.hh"

class G4NRESP71M03 
{
   public:
     
      G4NRESP71M03(){;}; 
      ~G4NRESP71M03(){;}; 

      void DKINMA(G4ReactionProduct *p1, G4ReactionProduct *p2, G4ReactionProduct *p3, G4ReactionProduct *p4, const G4double Q, const G4double costhcm3);

      G4int ApplyMechanismI_NBeA2A(G4ReactionProduct &neut, G4ReactionProduct &carb, G4ReactionProduct *theProds, const G4double QI);
      G4int ApplyMechanismII_ACN2A(G4ReactionProduct &neut, G4ReactionProduct &carb, G4ReactionProduct *theProds, const G4double QI);

      G4int ApplyMechanismABE(G4ReactionProduct &neut, G4ReactionProduct &carb, G4ReactionProduct *theProds);

   private:
      // Defining the arrays with the angular distribution data 
      static const G4int ndist = 32; // Number of angular distributions.
      static const G4int nrhos = 51; // Number of Rho values.
      // Energies for which angular distributions are given.
      static const G4double BEN2[ndist];
      // Angular distribution of alpha particles for each energy value in BEN2.
      static const G4double B2[ndist][nrhos];

};

#endif
