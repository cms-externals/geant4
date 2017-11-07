#include "G4LENDCombinedModel.hh"
#include "G4LENDCombinedCrossSection.hh"
#include "G4LENDElastic.hh"
#include "G4LENDInelastic.hh"
#include "G4LENDCapture.hh"
#include "G4LENDFission.hh"
#include "G4DynamicParticle.hh"

G4LENDCombinedModel::G4LENDCombinedModel( G4ParticleDefinition* pd )
:G4LENDModel( "LENDCombinedModel" ) {
   proj = pd;
   crossSection = new G4LENDCombinedCrossSection( pd );
   elastic = new G4LENDElastic( pd );
   inelastic = new G4LENDInelastic( pd );
   capture = new G4LENDCapture( pd );
   fission = new G4LENDFission( pd );
   channels[0] = elastic;
   channels[1] = inelastic;
   channels[2] = capture;
   channels[3] = fission;
}

void G4LENDCombinedModel::BuildPhysicsTable(const G4ParticleDefinition& projectile) {
   crossSection->BuildPhysicsTable( projectile );
   create_used_target_map();
}

G4bool G4LENDCombinedModel::HasData( const G4DynamicParticle* , G4int iZ, G4int iA , G4int iM,
                                     const G4Isotope* , const G4Element* , const G4Material* ) {

   G4bool result = false;
   if ( get_target_from_map ( lend_manager->GetNucleusEncoding ( iZ, iA , iM ) ) != NULL ) result=true;
   return result;
}


G4HadFinalState * G4LENDCombinedModel::ApplyYourself(const G4HadProjectile& aTrack, G4Nucleus& aTarg ){

   G4LENDModel* channel = NULL;

   G4int iZ = aTarg.GetZ_asInt();
   G4int iA = aTarg.GetA_asInt();
   //To pass kinetic energy, need to generate dynamic particle  
   G4DynamicParticle* dp = new G4DynamicParticle( proj , G4ThreeVector(0.,0.,1.) , aTrack.GetKineticEnergy() );
   G4int ichannel = crossSection->SelectChannel( dp , iZ , iA , aTarg.GetIsotope(), NULL , aTrack.GetMaterial() );
   delete dp;
   //ichannel=1;
   channel = channels[ichannel];
   return channel->ApplyYourself(aTrack,aTarg);
}
