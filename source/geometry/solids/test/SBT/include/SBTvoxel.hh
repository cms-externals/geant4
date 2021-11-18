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
// SBTvoxel.hh
//
// Definition of a batch based voxel test
//
#ifndef SBTvoxel_hh
#define SBTvoxel_hh

#include "globals.hh"
#include "G4ThreeVector.hh"
#include "geomdefs.hh"
#include <iostream>

class G4VSolid;
class G4VoxelLimits;
class G4AffineTransform;

class SBTVisManager;

class SBTvoxel {
	public:
	SBTvoxel();
	~SBTvoxel();
	
	void SetDefaults();
	
	void RunTest( const G4VSolid *testVolume, std::ostream &logger );
	
	void Draw( const G4VSolid *testVolume, 
		   const G4VoxelLimits &voxel, const G4AffineTransform &transform,
		   const G4ThreeVector *points, const G4int numPoints,
		   const EAxis	axis, const G4double limits[2],
	 	   SBTVisManager *visManager ) const;

	void Debug( const G4VSolid *testVolume, const EAxis axis,
		    const G4VoxelLimits &voxel, const G4AffineTransform &transform,
		    const G4ThreeVector *point ) const;
	


	void SetTarget( const G4ThreeVector &newTarget ) { target = newTarget; }
	G4ThreeVector GetTarget() const { return target; }

	void SetMaxVoxels( const G4int newMaxVoxels ) { maxVoxels = newMaxVoxels; }
	G4int GetMaxVoxels() const { return maxVoxels; }
	
	void SetMaxErrors( const G4int newMaxErrors ) { maxErrors = newMaxErrors; }
	G4int GetMaxErrors() const { return maxErrors; }

	inline void SetWidths( const G4ThreeVector &newWidths ) { widths = newWidths; }
	inline G4ThreeVector GetWidths() const { return widths; }


	protected:
	G4ThreeVector	target,
			widths;
	G4int		maxVoxels,
			maxErrors;
	
	G4bool TestOneVoxel( const G4VSolid *testVolume, 
			     const G4VoxelLimits &voxel,
			     const G4AffineTransform &transform,
			     const G4ThreeVector inside[], const G4int numInside,
			     std::ostream &logger ) const;

	void GetInsidePoints( const G4VSolid *testVolume,
			      G4ThreeVector inside[], G4int *numInside,
			      const G4int numPoints,
			      const G4int numAttempts ) const;

	G4ThreeVector GetRandomPoint() const;

	G4double GaussianRandom(const G4double cutoff) const;

	G4bool GetRandomLimit( G4double x[2], const G4double range ) const;

	G4VoxelLimits *NewRandomVoxel( const G4ThreeVector &theWidths ) const;
	
	void DumpVoxel( const G4VoxelLimits &voxel, std::ostream &logger ) const;
	
	void DumpTransform( const G4AffineTransform &transform, std::ostream &logger ) const;
	
	void MakeVoxelTestPoints( const G4VoxelLimits &voxel,
				  const EAxis axis, const G4double value,
				  G4ThreeVector testPoints[9] ) const;
};

#endif
