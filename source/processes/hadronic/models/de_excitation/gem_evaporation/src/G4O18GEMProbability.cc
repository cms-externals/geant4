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
// $Id: G4O18GEMProbability.cc,v 1.5 2006/06/29 20:23:31 gunter Exp $
// GEANT4 tag $Name: geant4-09-01 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1999)
//


#include "G4O18GEMProbability.hh"

G4O18GEMProbability::G4O18GEMProbability() :
  G4GEMProbability(18,8,0.0) // A,Z,Spin
{

  ExcitEnergies.push_back(1982.*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(2.6*picosecond);

  ExcitEnergies.push_back(3552.9*keV);
  ExcitSpins.push_back(4.0);
  ExcitLifetimes.push_back(2.0*picosecond);

  ExcitEnergies.push_back(3631.7*keV);
  ExcitSpins.push_back(0.0);
  ExcitLifetimes.push_back(2.1*picosecond);

  ExcitEnergies.push_back(3919.1*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(0.12*picosecond);

  ExcitEnergies.push_back(4448.8*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(0.08*picosecond);

  ExcitEnergies.push_back(7620.0*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(2.5*keV));

  ExcitEnergies.push_back(8039.0*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(2.5*keV));

  ExcitEnergies.push_back(8213.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(1.0*keV));

  ExcitEnergies.push_back(8283.0*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(8.0*keV));

  ExcitEnergies.push_back(10119.0*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(16.0*keV));

}


G4O18GEMProbability::G4O18GEMProbability(const G4O18GEMProbability &) : G4GEMProbability()
{
  throw G4HadronicException(__FILE__, __LINE__, "G4O18GEMProbability::copy_constructor meant to not be accessable");
}




const G4O18GEMProbability & G4O18GEMProbability::
operator=(const G4O18GEMProbability &)
{
  throw G4HadronicException(__FILE__, __LINE__, "G4O18GEMProbability::operator= meant to not be accessable");
  return *this;
}


G4bool G4O18GEMProbability::operator==(const G4O18GEMProbability &) const
{
  return false;
}

G4bool G4O18GEMProbability::operator!=(const G4O18GEMProbability &) const
{
  return true;
}


