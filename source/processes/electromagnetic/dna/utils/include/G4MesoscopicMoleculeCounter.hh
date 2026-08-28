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
// Author: Christian Velten (2026)

#ifndef G4MESOSCOPICMOLECULECOUNTER_HH
#define G4MESOSCOPICMOLECULECOUNTER_HH

#include "G4DNAMesh.hh"
#include "G4VUserMoleculeCounter.hh"

//------------------------------------------------------------------------------

struct G4MesoscopicMoleculeCounterIndex : public G4VMoleculeCounter::G4VMoleculeCounterIndex
{
    const G4MolecularConfiguration* Molecule{nullptr};
    G4float resolution{0.};  // mesh.GetResolution() at snapshot time; disambiguates epochs
    G4int linear{0};  // ix + nx*iy + nx*ny*iz

    G4MesoscopicMoleculeCounterIndex() = default;
    G4MesoscopicMoleculeCounterIndex(const G4MolecularConfiguration* mol) : Molecule(mol) {}
    G4MesoscopicMoleculeCounterIndex(const G4MolecularConfiguration* mol, G4float res, G4int lin)
      : Molecule(mol), resolution(res), linear(lin)
    {}
    ~G4MesoscopicMoleculeCounterIndex() override = default;

    G4bool operator<(G4VMoleculeCounterIndex const& rhs) const override
    {
      const auto& o = static_cast<const G4MesoscopicMoleculeCounterIndex&>(rhs);
      if (Molecule != o.Molecule) return std::less<>{}(Molecule, o.Molecule);
      if (resolution != o.resolution) return resolution < o.resolution;
      return linear < o.linear;
    }
    G4bool operator==(G4VMoleculeCounterIndex const& rhs) const override
    {
      const auto& o = static_cast<const G4MesoscopicMoleculeCounterIndex&>(rhs);
      return Molecule == o.Molecule && resolution == o.resolution && linear == o.linear;
    }
    G4String GetInfo() const override;
    const G4MolecularConfiguration* GetMolecule() const override { return Molecule; }
};

//------------------------------------------------------------------------------

class G4MesoscopicMoleculeCounter : public G4VUserMoleculeCounter<G4MesoscopicMoleculeCounterIndex>,
                                    public G4VMoleculeCounterMeshSupport
{
  public:

    G4MesoscopicMoleculeCounter();
    explicit G4MesoscopicMoleculeCounter(G4String);
    ~G4MesoscopicMoleculeCounter() override = default;

    void InitializeUser() override;
    void SetMeshSnapshot(const G4DNAMesh* const, G4double) override;

    G4int GetNbMoleculesAtTime(const G4MesoscopicMoleculeCounterIndex&, G4double) const override;
    G4int GetNbMoleculesAtTime(Search&, const G4MesoscopicMoleculeCounterIndex&,
                               G4double) const override;

  public:

    // Build{*}Index -> null; Add/Remove{*} are no-ops -> {}
    std::unique_ptr<G4VMoleculeCounterIndex> BuildIndex(const G4Track*) const override
    {
      return nullptr;
    }
    std::unique_ptr<G4VMoleculeCounterIndex> BuildIndex(const G4Track*,
                                                        const G4StepPoint*) const override
    {
      return nullptr;
    }
    std::unique_ptr<G4VMoleculeCounterIndex>
    BuildSimpleIndex(const G4MolecularConfiguration*) const override
    {
      return nullptr;
    }
    void AddMolecule(std::unique_ptr<G4VMoleculeCounterIndex>, G4double, G4int = 1) override {}
    void RemoveMolecule(std::unique_ptr<G4VMoleculeCounterIndex>, G4double, G4int = 1) override {}

  private:

    G4bool fRecordMeshData{false};

  public:

    G4bool GetRecordMeshData() const { return fRecordMeshData; }
    void SetRecordMeshData(G4bool flag = true) { fRecordMeshData = flag; }
};

#endif
