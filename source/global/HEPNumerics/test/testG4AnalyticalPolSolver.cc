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
// Test program for G4AnalyticalPolSolver class. 
//

#include "G4ios.hh"
#include "globals.hh"
#include "Randomize.hh"
#include "G4SystemOfUnits.hh"
#include "G4AnalyticalPolSolver.hh"
#include "geomdefs.hh"

// #include "ApproxEqual.hh"

const G4double kApproxEqualTolerance = 1E-6;
// const G4double kApproxEqualTolerance = 1E-2;

// Return true if the double check is approximately equal to target
//
// Process:
//
// Return true is check if less than kApproxEqualTolerance from target

G4bool ApproxEqual(const G4double check,const G4double target)
{
  G4bool result;
  G4double mean, delta;
  mean  = 0.5*std::fabs(check + target);
  delta =     std::fabs(check - target);

  if(mean > 1.) delta /= mean;

  if(delta<kApproxEqualTolerance) result = true;
  else                            result = false;
  return result;
}


int main()
{
  G4int i, k, n, iMax = 10000;
  G4int iCheck = iMax/10;
  G4double p[5], r[3][5];
  G4double a, b, c, d, tmp, range = 10*mm;
  G4AnalyticalPolSolver solver;
  enum Eroot {k2, k3, k4};
  Eroot useCase = k4;

  G4cout.precision(20);


  a    = 14.511252641677856;
  b    = 14.7648024559021;
  c    = 14.82865571975708;
  d    = 14.437621831893921;


  p[0] =  1.;
  p[1] = -a - b - c - d;
  p[2] =  (a+b)*(c+d) + a*b + c*d;
  p[3] = -(a+b)*c*d - a*b*(c+d);
  p[4] =  a*b*c*d;

  // iRoot = solver.BiquadRoots(p,r);
  solver.QuarticRoots(p,r);

  for( k = 1; k <= 4; k++ )
  { 
    tmp = r[1][k];

    if ( ApproxEqual(tmp,a) || ApproxEqual(tmp,b) || 
         ApproxEqual(tmp,c) || ApproxEqual(tmp,d) )  continue;
    else
    {
      G4cout<<"k    = "<<k<<G4endl;
      G4cout<<"a    = "<<a<<G4endl;
      G4cout<<"b    = "<<b<<G4endl;
      G4cout<<"c    = "<<c<<G4endl;
      G4cout<<"d    = "<<d<<G4endl;
      G4cout<<"root = "<< r[1][k] << " " << r[2][k] <<" i" << G4endl<<G4endl;
    }
  }


  for ( n = 2; n <= 4; n++ )
  {
    // Various test cases 

    if( n == 4 )   // roots: 1,2,3,4 
    { 
      p[0] =  1.; 
      p[1] = -10.; 
      p[2] =  35.; 
      p[3] = -50.; 
      p[4] =  24.; 

      p[0] =  1.; 
      p[1] =  0.; 
      p[2] =  4.; 
      p[3] =  0.; 
      p[4] =  4.; // roots: +-i sqrt(2) 
    }
    if( n == 3 )  // roots:  1,2,3 
    { 
      p[0] =  1.; 
      p[1] = -6.; 
      p[2] = 11.; 
      p[3] = -6.; 
    }
    if(n==2)  // roots : 1 +- i 
    { 
      p[0] =  1.; 
      p[1] = -2.; 
      p[2] =  2.;
    }

    if( n == 2 )
    {
      // G4cout<<"Test QuadRoots(p,r):"<<G4endl;    
      i = solver.QuadRoots(p,r);
    }
    else if( n == 3 )
    {
      // G4cout<<"Test CUBICROOTS(p,r):"<<G4endl;    
      i = solver.CubicRoots(p,r);
    }
    else if( n == 4 )
    {
      // G4cout<<"Test BIQUADROOTS(p,r):"<<G4endl;    
      i = solver.BiquadRoots(p,r);
    }

    for( k = 1; k <= n; k++ )
    { 
      // G4cout << r[1][k] << " " << r[2][k] <<" i" << G4endl;
    }
  }

  G4cout << G4endl << G4endl;

  // Random test of quadratic, cubic, and biquadratic equations

  switch (useCase)
  {

    case k4:
    G4cout<<"Testing biquadratic:"<<G4endl<<G4endl;
  
    for( i = 0; i < iMax; i++ )
    {
      if(i%iCheck == 0) G4cout<<"i = "<<i<<G4endl<<G4endl;

      a = -range + 2*range*G4UniformRand();
      b = -range + 2*range*G4UniformRand();
      c = -range + 2*range*G4UniformRand();
      d = -range + 2*range*G4UniformRand();
      p[0] =  1.;
      p[1] = -a - b - c - d;
      // p[2] =  a*b + a*c + a*d + b*c + b*d + c*d;
      p[2] =  (a+b)*(c+d) + a*b + c*d;
      // p[3] = -a*b*c - b*c*d - a*b*d - a*c*d;
      p[3] = -(a+b)*c*d - a*b*(c+d);
      p[4] =  a*b*c*d;

      // iRoot = solver.BiquadRoots(p,r);
      solver.QuarticRoots(p,r);
      for( k = 1; k <= 4; k++ )
      { 
        tmp = r[1][k];

        if ( ApproxEqual(tmp,a) || ApproxEqual(tmp,b) || 
             ApproxEqual(tmp,c) || ApproxEqual(tmp,d) )  continue;
        else
        {
          G4cout<<"i = "<<i<<";    k = "<<k<<G4endl;
          G4cout<<"a    = "<<a<<G4endl;
          G4cout<<"b    = "<<b<<G4endl;
          G4cout<<"c    = "<<c<<G4endl;
          G4cout<<"d    = "<<d<<G4endl;
          G4cout<<"root = "<< r[1][k] << " " << r[2][k] <<" i"
                << G4endl << G4endl;
        }
      }
    }
    break;

    case k3:

    G4cout<<"Testing cubic:"<<G4endl<<G4endl;
  
    for( i = 0; i < iMax; i++ )
    {
      if(i%iCheck == 0) G4cout<<"i = "<<i<<G4endl;

      a = -range + 2*range*G4UniformRand();
      b = -range + 2*range*G4UniformRand();
      c = -range + 2*range*G4UniformRand();

      p[0] =  1.;
      p[1] = -a - b - c;
      p[2] =  (a+b)*c + a*b;
      p[3] =  -a*b*c;

      solver.CubicRoots(p,r);

      for( k = 1; k <= 3; k++ )
      { 
        tmp = r[1][k];

        if ( ApproxEqual(tmp,a) || ApproxEqual(tmp,b) || 
             ApproxEqual(tmp,c) )  continue;
        else
        {
          G4cout<<"i = "<<i<<G4endl;
          G4cout<<"k = "<<k<<G4endl;
          G4cout<<"a = "<<a<<G4endl;
          G4cout<<"b = "<<b<<G4endl;
          G4cout<<"c = "<<c<<G4endl;
          G4cout <<"root = "<< r[1][k] << " " << r[2][k] <<" i" << G4endl;
        }
      }
    }
    break;
  case k2:

    G4cout<<"Testing quadratic:"<<G4endl<<G4endl;
  
    for( i = 0; i < iMax; i++ )
    {
      if(i%iCheck == 0) G4cout<<"i = "<<i<<G4endl;

      a = -range + 2*range*G4UniformRand();
      b = -range + 2*range*G4UniformRand();

      p[0] =  1.;
      p[1] = -a - b ;
      p[2] =  a*b;

      solver.QuadRoots(p,r);

      for( k = 1; k <= 2; k++ )
      { 
        tmp = r[1][k];

        if ( ApproxEqual(tmp,a) || ApproxEqual(tmp,b) )  continue;
        else
        {
          G4cout<<"i = "<<i<<G4endl;
          G4cout<<"k = "<<k<<G4endl;
          G4cout<<"a = "<<a<<G4endl;
          G4cout<<"b = "<<b<<G4endl;
          G4cout <<"root = "<< r[1][k] << " " << r[2][k] <<" i" << G4endl;
        }
      }
    }
    break;

    default:
    break;   
  }
  return 0;
}
