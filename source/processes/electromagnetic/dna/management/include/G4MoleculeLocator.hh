#ifndef G4MOLECULELOCATOR_HH
#define G4MOLECULELOCATOR_HH 1
#pragma once

#include "G4ThreeVector.hh"
#include "G4TouchableHandle.hh"

#include <memory>

class G4Navigator;
class G4TouchableHistory;
class G4VPhysicalVolume;

class G4MoleculeLocator final
{
public:
  static G4MoleculeLocator* Instance();

private:
  G4MoleculeLocator();
  ~G4MoleculeLocator() = default;

  G4MoleculeLocator(const G4MoleculeLocator&) = delete;
  G4MoleculeLocator(G4MoleculeLocator&&) = delete;
  G4MoleculeLocator& operator=(const G4MoleculeLocator&) = delete;
  G4MoleculeLocator& operator=(G4MoleculeLocator&&) = delete;

  static G4ThreadLocal G4MoleculeLocator* fpInstance;

  G4bool fIsInitialized;
  std::unique_ptr<G4Navigator> fNavigator;

  void Initialize();

public:
  G4TouchableHistory* CreateTouchableHistory() const;

  G4VPhysicalVolume* LocateGlobalPointAndSetup(const G4ThreeVector& point,
                                               const G4ThreeVector* direction = nullptr);

  G4VTouchable* LocateGlobalPointAndReturnNewTouchable(const G4ThreeVector& position);
  G4VTouchable* LocateGlobalPointAndReturnNewTouchable(const G4ThreeVector& position,
                                                       const G4ThreeVector& direction);

  G4TouchableHandle LocateGlobalPointAndReturnNewTouchableHandle(const G4ThreeVector& position);
  G4TouchableHandle LocateGlobalPointAndReturnNewTouchableHandle(const G4ThreeVector& position,
                                                                 const G4ThreeVector& direction);

  void LocateGlobalPointAndUpdateTouchable(const G4ThreeVector& position, G4VTouchable* touchable);
  void LocateGlobalPointAndUpdateTouchable(const G4ThreeVector& position,
                                           const G4ThreeVector& direction, G4VTouchable* touchable);
};

#endif
