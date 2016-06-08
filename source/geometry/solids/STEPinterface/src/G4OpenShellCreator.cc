//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4OpenShellCreator.cc,v 1.4 2002/11/21 16:49:49 gcosmo Exp $
// GEANT4 tag $Name: geant4-05-00 $
//
// 
// ----------------------------------------------------------------------
// Class G4OpenShellCreator
//
// Authors: J.Sulkimo, P.Urban.
// Revisions by: L.Broglia, G.Cosmo.
//
// History:
//   18-Nov-1999: First step of re-engineering - G.Cosmo
// ----------------------------------------------------------------------

#include "G4OpenShellCreator.hh"
#include "G4GeometryTable.hh"

G4OpenShellCreator G4OpenShellCreator::csc;

G4OpenShellCreator::G4OpenShellCreator()
{
  G4GeometryTable::RegisterObject(this);
}

G4OpenShellCreator::~G4OpenShellCreator() {}

G4OpenShellCreator G4OpenShellCreator::GetInstance()
{
  return csc;
}

void G4OpenShellCreator::CreateG4Geometry(STEPentity& Ent)
{
}

void G4OpenShellCreator::CreateSTEPGeometry(void* G4obj)
{
}