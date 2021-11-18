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
// Test program for G4GaussianQuadrature class. The function std::exp(-x)*std::cos(x) is
// integrated between zero and two pi. The true result is 0.499066278634
//

#include "G4ios.hh"
#include "globals.hh"
#include "G4PhysicalConstants.hh"
#include "G4GaussChebyshevQ.hh"
#include "G4GaussHermiteQ.hh"
#include "G4GaussJacobiQ.hh"
#include "G4GaussLaguerreQ.hh"
#include "G4GaussLegendreQ.hh"

G4double TestChebyshev(G4double x)
{
  return std::sqrt(1-x*x)*std::cos(x) ;
}

G4double TestFunction(G4double x)
{
  return std::exp(-x)*std::cos(x) ;
}

G4double TestHermite(G4double x)
{
  return x*x*std::cos(x) ;
}

G4double CosFunction(G4double x)
{
  return std::cos(x) ;
}

int main()
{
    G4int i, n = 20;
   // G4double pTolerance ;
   G4double a = -1.0 ;
   G4double b = 1.0 ;
   G4double truth = pi*0.4400505857 ;
   for(i=1;i<20;i++)
   {
      n = 1*i ;
      G4GaussChebyshevQ myChebyshev(TestChebyshev,n) ;
      G4cout<<"n = "<<n<<"\t"<<"true = "<<truth<<"  and n-point Gauss-Chebyshev =  "
	  <<myChebyshev.Integral(a,b)<<G4endl ;
   }
   for(i=1;i<20;i++)
   {
      n = 1*i ;
      G4GaussJacobiQ myJacobi(CosFunction,0.5,0.5,n) ;
      G4cout<<"n = "<<n<<"\t"<<"true = "<<truth<<"  and n-point Gauss-Jacobi =  "
	  <<myJacobi.Integral()<<G4endl ;
   }
   truth = 2*0.125*std::sqrt(pi)*std::exp(-0.25) ;
   for(i=1;i<20;i++)
   {
      n = 1*i ;
      G4GaussHermiteQ myHermite(TestHermite,n) ;
      G4cout<<"n = "<<n<<"\t"<<"true = "<<truth<<"  and n-point GaussHermite =  "
	  <<myHermite.Integral()<<G4endl ;
   }
   /* *******************
   G4GaussianQuadrature myHermite(n) ;
   for(i=0;i<(n+1)/2;i++)
   {
      G4cout<<i<<"\t"<<myHermite.GetAbscissa(i)<<"\t"
	  <<myHermite.GetWeight(i)<<G4endl ;
   }
   G4GaussianQuadrature myLaguerre(0.0,n) ;
   for(i=0;i<n;i++)
   {
      G4cout<<i<<"\t"<<myLaguerre.GetAbscissa(i)<<"\t"
	  <<myLaguerre.GetWeight(i)<<G4endl ;
   }
   */ /////////////////////////////////////
   for(i=1;i<20;i++)
   {
      n = 1*i ;
      G4GaussLaguerreQ myLaguerre(CosFunction,0.0,n) ;
      G4cout<<"n = "<<n<<"\t"<<"true = 0.5 "<<"  and n-point GaussLaguerre =  "
	  <<myLaguerre.Integral()<<G4endl ;
   }
   G4GaussLegendreQ myIntegrand(TestFunction) ;
   G4cout<<"true is 0.499066278634 "<<"  and QuickGaussLegendre is  "<<
      myIntegrand.QuickIntegral(0,2*pi)<<G4endl ;
   G4cout<<"true is 0.499066278634 "<<"  and AccurateGaussLegendre is  "<<
      myIntegrand.AccurateIntegral(0,2*pi)<<G4endl ;


   for(i=1;i<20;i++)
   {
      n = 8*i ;
      G4GaussLegendreQ myLegendre(TestFunction, n) ;
      G4cout<<myLegendre.GetNumber()<<
      "true is 0.5 "<<"  and n-point GaussLegendre is  "
	  <<myLegendre.Integral(0,200*pi)<<G4endl ;
   }
   /* **************************************************
   
   */ /////////////////////////////////////
  return 0;
}
