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
// FredHit.hh
//
// Definition of FredHit: intersection of a track with a test volume
//

#ifndef FredHit_HH
#define FredHit_HH

//
// Geez... this stuff would be very hard without examples...
//
#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4Color.hh"
#include "G4VisAttributes.hh"

class G4Track;

class FredHit : public G4VHit
{
	public:
	FredHit() {;}
	~FredHit() {;}
	FredHit( const G4ThreeVector pos, G4bool enters, const G4Track *track );
	
	//
	// Tsk Tsk. c++ is such hell.
	// Define copy and comparison. Note that hits are never equal.
	FredHit( const FredHit &right );
	const FredHit& operator=( const FredHit &right );
	int operator==(const FredHit &) const { return 0; }
	
	//
	// See below. We use G4's special paging allocator
	inline void * operator new(size_t);
	inline void operator delete( void *hit );
	
	//
	// These G4VHit virtual functions must be overridden, otherwise
	// c++ gets upset. We wouldn't want that.
	void Draw();
	void Print();
	
	//
	// Special draw modes
	void DrawShadowHit();
	void DrawErrorHit();
	void DrawMaskHit();
	
	//
	// As tempting as it is, we should keep data structures private,
	// and write public routines to access them. Otherwise,
	// flexibility suffers.
	private:
	G4ThreeVector pos;
	G4bool enters;
	const G4Track *track;
	
	// Draw primitive
	void DrawPrim(G4int imode);

        public:
	void SetPos( const G4ThreeVector newPos ) { pos = newPos; }
	void SetEnters( G4bool newEnters ) { enters = newEnters; }
	void SetTrack( const G4Track *newTrack ) { track = newTrack; }
	
	G4ThreeVector GetPos() const { return pos; }
	G4bool GetEnters() const { return enters; }
	const G4Track *GetTrack() const { return track; }
};

//
// G4 likes to play with "collections" of hits, in a standard manner.
// The supplied template class allows us to build this standard collection.
// Not that I would *dare* challenge G4's designers, but couldn't they
// have defined a hit collection as a collection of the G4VHit base
// class??
typedef G4THitsCollection<FredHit> FredHitCollection;

//
// G4 has it's own paging allocator. Good idea for a program of this
// complexity. It also allows all FredHit instances to be erased with
// one line of code. Kinda dangerous.
//
// Normally you'd put these inline codes inside the class declaration,
// but we can hardly do that, since we need to specify the allocator
// template. But, was "inline" really necessary?
extern G4Allocator<FredHit> FredHitAllocator;

inline void *FredHit::operator new(size_t) 
{
	void *hit = (void *)FredHitAllocator.MallocSingle();
	return hit;
}

inline void FredHit::operator delete( void *hit )
{
	FredHitAllocator.FreeSingle( (FredHit *)hit );
}

#endif
