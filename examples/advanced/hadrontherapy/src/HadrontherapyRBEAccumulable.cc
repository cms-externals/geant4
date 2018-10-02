#include "HadrontherapyRBEAccumulable.hh"
#include "HadrontherapyRBE.hh"

#include <tuple>
#include <G4SystemOfUnits.hh>

using namespace std;

HadrontherapyRBEAccumulable::HadrontherapyRBEAccumulable()
    : G4VAccumulable("RBE")
{

}

void HadrontherapyRBEAccumulable::Merge(const G4VAccumulable& rhs)
{
    if (GetVerboseLevel() > 1)
    {
        G4cout << "HadrontherapyRBEAccumulable::Merge()" << G4endl;
    }
    const HadrontherapyRBEAccumulable& other = dynamic_cast<const HadrontherapyRBEAccumulable&>(rhs);
    fAlphaNumerator += other.fAlphaNumerator;
    fDenominator += other.fDenominator;
    fBetaNumerator += other.fBetaNumerator;
    fEnergyDeposit += other.fEnergyDeposit;
}

void HadrontherapyRBEAccumulable::Reset()
{
    if (GetVerboseLevel() > 1)
    {
        G4cout << "HadrontherapyRBEAccumulable::Reset()" << G4endl;
    }
    if (fInitialized)
    {
        fAlphaNumerator = 0.0;
        fBetaNumerator = 0.0;
        fDenominator = 0.0;
        fEnergyDeposit = 0.0;
    }
    else
    {
        Initialize();
    }
}

void HadrontherapyRBEAccumulable::Accumulate(G4double E, G4double energyDeposit, G4double dX, G4int Z, G4int i, G4int j, G4int k)
{
    if (!fInitialized)
    {
        G4Exception("HadrontherapyRBEAccumulable::Accumulate", "NotInitialized", FatalException, "Accumulable not initialized. Must be a programming error.");
    }
    if (GetVerboseLevel() > 2)
    {
        G4cout << "HadrontherapyRBEAccumulable::Accumulate() in " << i << ", " << j << ", " << k << G4endl;
    }
    if (energyDeposit <= 0)
    {
        return;
    }
    size_t n = GetIndex(i, j, k);
    fEnergyDeposit[n] += energyDeposit;

    if ((Z >= 1) && (dX > 0) && (E > 0)) // TODO: Verify this condition
    {
        tuple<G4double, G4double> alpha_beta = HadrontherapyRBE::GetInstance()->GetHitAlphaAndBeta(E, Z);
        fDenominator[n] += energyDeposit;
        fAlphaNumerator[n] += get<0>(alpha_beta) * energyDeposit;
        fBetaNumerator[n] += sqrt(get<1>(alpha_beta)) * energyDeposit;
    }
}

const HadrontherapyRBEAccumulable::array_type HadrontherapyRBEAccumulable::GetEnergyDeposit() const
{
    return fEnergyDeposit;
}

G4int HadrontherapyRBEAccumulable::GetVerboseLevel() const
{
    return HadrontherapyRBE::GetInstance()->GetVerboseLevel();
}

void HadrontherapyRBEAccumulable::Initialize()
{
    if (GetVerboseLevel() > 1)
    {
        G4cout << "HadrontherapyRBEAccumulable::Initialize(): ";
    }
    auto rbe = HadrontherapyRBE::GetInstance();

    fVoxelsAlongX = rbe->GetNumberOfVoxelsAlongX();
    fVoxelsAlongY = rbe->GetNumberOfVoxelsAlongY();
    fVoxelsAlongZ = rbe->GetNumberOfVoxelsAlongZ();
    fVoxels = fVoxelsAlongX * fVoxelsAlongY * fVoxelsAlongZ;

    if (GetVerboseLevel() > 1)
    {
        G4cout << fVoxels << " voxels." << G4endl;
    }

    fAlphaNumerator = array_type(0.0, fVoxels);
    fBetaNumerator = array_type(0.0, fVoxels);
    fDenominator = array_type(0.0, fVoxels);
    fEnergyDeposit = array_type(0.0, fVoxels);
    fInitialized = true;
}
