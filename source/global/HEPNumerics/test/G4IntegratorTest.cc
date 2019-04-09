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
// Test program for G4Integrator class. The function std::exp(-x)*std::cos(x) is
// integrated between zero and two pi. The exact result is 0.499066278634
//
// History:
//
// 05.04.97 V.Grichine, first implementation
// 04.09.99 V.Grichine, integration of member function from a class and main,
//                      as well as integration of global scope functions

#include "G4ios.hh"
#include "globals.hh"
#include "G4PhysicalConstants.hh"
#include "G4SimpleIntegration.hh"
#include "G4Integrator.hh"

G4double GlobalFunction( G4double x ){  return std::exp(-x)*std::cos(x) ; }

G4double GlobalCos( G4double x ){  return std::cos(x) ; }

G4double GlobalHermite(G4double x){  return x*x*std::cos(x) ; }


G4double fY;
G4int    fN = 100;

G4double X1 = -3.25/5.;
G4double X2 =  3.25/5.;

G4double Y1 = -7.5/5.;
G4double Y2 =  7.5/5.;

G4double Harp(G4double  x)
{

  G4double tmp = std::sqrt(1 +  x*x + fY*fY);
   tmp        *= 1 +  x*x + fY*fY ;
   return        1/tmp;
}


G4double HarpX(G4double y)
{
  fY = y;
  G4SimpleIntegration myIntegrand(Harp);
  return myIntegrand.Simpson(X1,X2,fN); 
  // return myIntegrand.Gauss(X1,X2,fN); 
}






G4double HarpY()
{
  // G4int i;
  G4SimpleIntegration myIntegrand(HarpX);
  // G4double sum =0.;
  // for(i = 0; i < fN; i++)
  G4double result = myIntegrand.Simpson(Y1,Y2,fN); 
  return result/twopi;

}

class B
{
  typedef G4double (B::* PBmem)(G4double);

public:

  B(){;}

 ~B(){;}

  G4double TestFunction(G4double x){  return std::exp(-x)*std::cos(x) ; }

  G4double CosFunction(G4double x) {  return std::cos(x) ; }

  G4double TestHermite(G4double x){  return x*x*std::cos(x) ; }

  G4double HarpX(G4double);
  G4double HarpY(G4double);


  void Integrand() ;

};


void B::Integrand()
{
     G4int i, n ;
     G4double pTolerance;
     G4double simpson1=0., simpson2=0. ;
     G4double legendre1=0., legendre2=0. ;
     G4double chebyshev1=0., chebyshev2=0. ;
     G4double adaptg1=0., adaptg2=0.;
     G4double a = 0.0 ;
     G4double b = twopi ;
      
     G4SimpleIntegration myIntegrand(GlobalFunction) ;

     G4Integrator<B,PBmem> integral;

     B bbb ;

     G4cout<<"Iter"
           <<"\t"<<"Simpson "<<"\t"<<"Simpson "
           <<"\t"<<"Legendre"<<"\t"<<"Legendre"
           <<"\t"<<"Chebyshev"<<"\t"<<"Chebyshev"<<G4endl ;

     for(i=1;i<=20;i++)
     {
        n = 2*i ;

        simpson1 = integral.Simpson(this,&B::TestFunction,a,b,n) ;
	legendre1 = integral.Legendre(this,&B::TestFunction,a,b,n) ;
	simpson2 = integral.Simpson(bbb,&B::TestFunction,a,b,n) ;
	legendre2 = integral.Legendre(bbb,&B::TestFunction,a,b,n) ;
	chebyshev1 = integral.Chebyshev(this,&B::TestFunction,a,b,n) ;
	chebyshev2 = integral.Chebyshev(bbb,&B::TestFunction,a,b,n) ;

        G4cout<<n
	      <<"\t"<<simpson1<<"\t"<<simpson2
              <<"\t"<<legendre1<<"\t"<<legendre2
              <<"\t"<<chebyshev1<<"\t"<<chebyshev2<<G4endl ;
     }
     G4cout<<G4endl ;

     for(i=0;i<8;i++)
     {
       pTolerance = std::pow(10.0,-i) ;
       adaptg1 = integral.AdaptiveGauss(bbb,&B::TestFunction,a,b,pTolerance) ; 
       adaptg2 = integral.AdaptiveGauss(this,&B::TestFunction,a,b,pTolerance) ; 
       G4cout<<pTolerance<<"\t"<<adaptg1<<"\t"<<adaptg2<<G4endl;
     }
   for(i=1;i<20;i++)
   {
      n = 1*i ;
      G4double laguerre1=0., laguerre2=0.;
      laguerre1 = integral.Laguerre(bbb,&B::CosFunction,0.0,n) ;
      laguerre2 = integral.Laguerre(this,&B::CosFunction,0.0,n) ;
      G4cout<<"n = "<<n<<"\t"<<"exact = 0.5 "
            <<"  and n-point GaussLaguerre =  "
	    <<laguerre1<<"\t"<<laguerre2<<G4endl ;
   }
   for(i=1;i<20;i++)
   {
      n = 1*i ;
      G4double hermite1=0., hermite2=0;
      G4double exactH = 2*0.125*std::sqrt(pi)*std::exp(-0.25) ;
      hermite1 = integral.Hermite(bbb,&B::TestHermite,n) ;
      hermite2 = integral.Hermite(this,&B::TestHermite,n) ;
      G4cout<<"n = "<<n<<"\t"<<"exact = "<<exactH
            <<"  and n-point GaussHermite =  "
	    <<hermite1<<"\t"<<hermite2<<G4endl ;
   }
   G4double exactJ = pi*0.4400505857 ;

   for(i=1;i<20;i++)
   {
      n = 1*i ;
      G4double jacobi1=0., jacobi2=0.;
      jacobi1 = integral.Jacobi(bbb,&B::CosFunction,0.5,0.5,n) ;
      jacobi2 = integral.Jacobi(this,&B::CosFunction,0.5,0.5,n) ;
      G4cout<<"n = "<<n<<"\t"<<"exact = "<<exactJ
            <<"  and n-point Gauss-Jacobi =  "
	    <<jacobi1<<"\t"<<jacobi2<<G4endl ;
   }
}


int main()
{
   B myIntegration ;

   myIntegration.Integrand() ;

   G4Integrator<B,function> iii ;
   G4int i, n ;
   G4double a = 0.0 ;
   G4double b = twopi ;
   G4double simpson3,legendre,legendre10,legendre96,chebyshev ;

   G4cout<<"Global function integration"<<G4endl ;
   G4cout<<"n = "<<"\t"<<"Simpson"<<"\t"
                 <<"\t"<<"Legendre""\t"<<"Chebyshev"<<G4endl ;
   for(i=1;i<=30;i++)
   {
     n = 2*i ;
     simpson3 = iii.Simpson(GlobalFunction,a,b,n) ;
     legendre = iii.Legendre(GlobalFunction,a,b,n) ;
     chebyshev = iii.Chebyshev(GlobalFunction,a,b,n) ;
     G4cout<<n<<"\t"<<simpson3<<"\t"<<legendre<<"\t"<<chebyshev<<G4endl ;
   }
   legendre10 = iii.Legendre10(GlobalFunction,a,b) ;
   legendre96 = iii.Legendre96(GlobalFunction,a,b) ;
   G4cout<<"Legendre 10 points = "<<"\t"<<legendre10<<G4endl ;
   G4cout<<"Legendre 96 points = "<<"\t"<<legendre96<<G4endl ;

   for(i=0;i<8;i++)
   {
     G4double  pTolerance = std::pow(10.0,-i) ;
     G4double  adaptg = iii.AdaptiveGauss(GlobalFunction,a,b,pTolerance) ;
     G4cout<<pTolerance<<"\t"<<adaptg<<G4endl;
   }
   for(i=1;i<20;i++)
   {
      n = 1*i ;
      G4double laguerre = iii.Laguerre(GlobalCos,0.0,n) ;
      G4cout<<"n = "<<n<<"\t"<<"exact = 0.5 "
            <<"  and n-point Laguerre =  "
	    <<laguerre<<G4endl ;
   }
   for(i=1;i<20;i++)
   {
      n = 1*i ;
      G4double exactH = 2*0.125*std::sqrt(pi)*std::exp(-0.25) ;
      G4double hermite = iii.Hermite(GlobalHermite,n) ;
      G4cout<<"n = "<<n<<"\t"<<"exact = "<<exactH
            <<"  and n-point Hermite =  "
	    <<hermite<<G4endl ;
   }
   G4double exactJ = pi*0.4400505857 ;

   for(i=1;i<20;i++)
   {
      n = 1*i ;
      G4double jacobi = iii.Jacobi(GlobalCos,0.5,0.5,n) ;
      G4cout<<"n = "<<n<<"\t"<<"exact = "<<exactJ
            <<"  and n-point Jacobi =  "
	    <<jacobi<<G4endl ;
   }

   G4cout<<"Integral = "<<HarpY()<<G4endl;

   return 0;
}
