#include "G4LENDorBERTModel.hh"
#include "G4LENDCombinedModel.hh"
#include "G4CascadeInterface.hh"
#include "G4DynamicParticle.hh"

G4LENDorBERTModel::G4LENDorBERTModel( G4ParticleDefinition* pd )
:G4LENDModel( "LENDorBERTModel" ) {
   proj = pd;
   lend = new G4LENDCombinedModel( proj ); 
   bert = new G4CascadeInterface;
}

void G4LENDorBERTModel::BuildPhysicsTable(const G4ParticleDefinition& projectile) {
   lend->BuildPhysicsTable( projectile );
   create_used_target_map();
}

G4HadFinalState * G4LENDorBERTModel::ApplyYourself(const G4HadProjectile& aTrack, G4Nucleus& aTarg ){
   //G4cout << "Hello from G4LENDorBERTModel::ApplyYourself" << G4endl;

   G4int iZ = aTarg.GetZ_asInt();
   G4int iA = aTarg.GetA_asInt();
   G4int iM = 0;
   if ( aTarg.GetIsotope() != NULL ) iM = aTarg.GetIsotope()->Getm();

   G4DynamicParticle* dp = new G4DynamicParticle( aTrack.GetDefinition() , G4ThreeVector(0.,0.,1.) , aTrack.GetKineticEnergy() );
   G4bool lendIsOK = lend->HasData( dp , iZ , iA , iM , aTarg.GetIsotope() , NULL , aTrack.GetMaterial() );
   delete dp;

   G4HadronicInteraction* model=NULL;
   if ( lendIsOK ) { 
      //G4cout << "LEND is selected" << G4endl;
      model = lend;
   } else { 
      //G4cout << "BERT is selected" << G4endl;
      model = bert;
   }

   return model->ApplyYourself(aTrack,aTarg);
}
