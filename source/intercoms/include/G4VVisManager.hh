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
// $Id: G4VVisManager.hh,v 1.7 2002/11/20 14:46:00 gcosmo Exp $
// GEANT4 tag $Name: geant4-05-00 $
// John Allison 19/Oct/1996.
// 
// Class Description:
//
// G4VVisManager is an abstract interface for the GEANT4 Visualization Manager.
// The inheritance hierarchy is:
//   G4VVisManager <- G4VisManager <- YourVisManager
//
// See example/novice/N02 to see how to write YourVisManager and
// instantiate it.  You should *not* access it directly; instead you
// should obtain a pointer as follows:
// 
//   G4VVisManager* pVVMan =  G4VVisManager::GetConcreteInstance ();
//
// This ensures your code will link even if YourVisManager is not
// instantiated or even if not provided in a library.  Please protect
// your code by testing the pointer, for example, by:
//
//   if (pVVMan) pVVMan -> Draw (polyline);
//
// The Draw functions draw only "transient" objects.  This is useful
// for debugging, e.g., drawing the step in your UserSteppingAction,
// since G4Steps are not kept.
//
// Note: to draw "permanent" objects, i.e., objects which are always
// available, such as detector geometry components, or available in an
// event after tracking has finished, such as hits, digitisations and
// trajectories, can be drawn in a transient way if you wish but it is
// usually possible to draw them in a permanent way with /vis/
// commands.  The advantage is that permanent objects can be redrawn,
// e.g., when you change view or viewer; transient objects get
// forgotten.
//
// Note that the G4Transform3D argument refers to the transformation
// of the *object*, not the transformation of the coordinate syste.
//
// Note also that where a G4VisAttributes argument is specified, it
// overrides any attributes belonging to the object itself.

#ifndef G4VVISMANAGER_HH
#define G4VVISMANAGER_HH

#include "G4Transform3D.hh"
#include "G4ThreeVector.hh"       // Just a typedef Hep3Vector.
#include "G4RotationMatrix.hh"    // Just a typedef HepRotation.

class G4Polyline;
class G4Text;
class G4Circle;
class G4Scale;
class G4Square;
class G4Polymarker;
class G4Polyhedron;
class G4NURBS;
class G4VSolid;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4VisAttributes;

class G4VVisManager {

public: // With description

  static G4VVisManager* GetConcreteInstance ();
  // Returns pointer to actual visualization manager if a view is
  // available for drawing, else returns null.  Always check value.

public:

  virtual ~G4VVisManager ();

public: // With description

  ///////////////////////////////////////////////////////////////////////
  // Draw methods for Geant4 Visualization Primitives, useful
  // for representing hits, digis, etc.

  virtual void Draw (const G4Circle&,
    const G4Transform3D& objectTransformation = G4Transform3D::Identity) = 0;

  virtual void Draw (const G4NURBS&,
    const G4Transform3D& objectTransformation = G4Transform3D::Identity) = 0;

  virtual void Draw (const G4Polyhedron&,
    const G4Transform3D& objectTransformation = G4Transform3D::Identity) = 0;

  virtual void Draw (const G4Polyline&,
    const G4Transform3D& objectTransformation = G4Transform3D::Identity) = 0;

  virtual void Draw (const G4Polymarker&,
    const G4Transform3D& objectTransformation = G4Transform3D::Identity) = 0;

  virtual void Draw (const G4Scale&,
    const G4Transform3D& objectTransformation = G4Transform3D::Identity) = 0;

  virtual void Draw (const G4Square&,
    const G4Transform3D& objectTransformation = G4Transform3D::Identity) = 0;

  virtual void Draw (const G4Text&,
    const G4Transform3D& objectTransformation = G4Transform3D::Identity) = 0;

  ///////////////////////////////////////////////////////////////////////
  // Draw methods for Geant4 Geometry Objects as if they were
  // Visualization Primitives.  Useful for representing hits.

  virtual void Draw (const G4LogicalVolume&, const G4VisAttributes&,
    const G4Transform3D& objectTransformation = G4Transform3D::Identity) = 0;

  virtual void Draw (const G4VPhysicalVolume&, const G4VisAttributes&,
    const G4Transform3D& objectTransformation = G4Transform3D::Identity) = 0;

  virtual void Draw (const G4VSolid&, const G4VisAttributes&,
    const G4Transform3D& objectTransformation = G4Transform3D::Identity) = 0;

  //////////////////////////////////////////////////////////////////////
  // Other methods...

  virtual void GeometryHasChanged () = 0;
  // This is used by the run manager to notify a change of geometry.

protected:

  static void SetConcreteInstance (G4VVisManager*);

  static G4VVisManager* fpConcreteInstance;  // Pointer to real G4VisManager.

};

#endif