#ifndef G4ScoringRealWorld_h
#define G4ScoringRealWorld_h 1

#include "globals.hh"
#include "G4VScoringMesh.hh"
class G4VPhysicalVolume;
class G4LogicalVolume;

class G4ScoringRealWorld : public G4VScoringMesh
{
public:
  G4ScoringRealWorld(G4String lvName);
  ~G4ScoringRealWorld();

protected:
  // construct this mesh
  virtual void SetupGeometry(G4VPhysicalVolume* );

protected:
  G4String logVolName;

public:
  virtual void List() const;

public:
    //++++++++++ visualization method not yet implemented
    virtual void Draw(RunScore * /*map*/, G4VScoreColorMap* /*colorMap*/, G4int /*axflg=111*/)
    {;}
    virtual void DrawColumn(RunScore * /*map*/, G4VScoreColorMap* /*colorMap*/,
                          G4int /*idxProj*/, G4int /*idxColumn*/)
    {;}
};




#endif
