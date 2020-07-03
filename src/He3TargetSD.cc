#include "He3TargetSD.hh"
//______________________________________________________________________________
He3TargetSD::He3TargetSD(G4String name,G4String hcName)
: G4VSensitiveDetector(name)
{
   collectionName.insert(hcName);
   fHitCollection = NULL;  
}
//______________________________________________________________________________
He3TargetSD::~He3TargetSD(){

}
//______________________________________________________________________________
void He3TargetSD::Initialize(G4HCofThisEvent *hce){
   fHitCollection = new He3TargetHitCollection(SensitiveDetectorName,collectionName[0]); 
}
//______________________________________________________________________________
G4bool He3TargetSD::ProcessHits(G4Step *aStep,G4TouchableHistory *){
   // get particle type and Edep 
   G4int pid   = aStep->GetTrack()->GetParticleDefinition()->GetPDGEncoding();
   double edep = aStep->GetTotalEnergyDeposit(); 

   if(edep<=0.0||pid==0) return false; 

   He3TargetHit *hit = new He3TargetHit(); 

   G4ThreeVector pos = aStep->GetPreStepPoint()->GetPosition();
   G4ThreeVector mom = aStep->GetPreStepPoint()->GetMomentum();

   G4double E        = aStep->GetPreStepPoint()->GetTotalEnergy();

   // G4TouchableHistory *hist = (G4TouchableHistory *)( aStep->GetPreStepPoint()->GetTouchable() );  
   // G4AffineTransform aTrans = hist->GetHistory()->GetTopTransform();

   hit->SetLabPos(pos); // global coordinates 
   hit->SetVertex( aStep->GetTrack()->GetVertexPosition() ); 
   hit->SetTime( aStep->GetPreStepPoint()->GetGlobalTime() ); 
   
   hit->SetEdep(edep); 
   hit->SetEnergy(E); 
   hit->SetLStep( aStep->GetStepLength() ); 
   hit->SetMomentum(mom); 
   hit->SetPID( aStep->GetTrack()->GetParticleDefinition()->GetPDGEncoding() ); 
   hit->SetTrID( aStep->GetTrack()->GetTrackID() ); 
   hit->SetMID( aStep->GetTrack()->GetParentID() );

   // G4Track *aTrack = aStep->GetTrack();

   // // hit->SetOriginalTrackInformation(aTrack);
   // // hit->SetPrimaryTrackInformation(aTrack);
   // // hit->SetSDTrackInformation( aTrack, fullPathName );
   // hit->SetOTrIdx( SDtracks.InsertOriginalTrackInformation( aTrack ) );
   // hit->SetPTrIdx( SDtracks.InsertPrimaryTrackInformation( aTrack ) );
   // hit->SetSDTrIdx( SDtracks.InsertSDTrackInformation( aTrack ) ); 
 
   fHitCollection->insert(hit);

   return true;   
}
//______________________________________________________________________________
void He3TargetSD::EndOfEvent(G4HCofThisEvent *hce){
   G4int hcid = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]); 
   hce->AddHitsCollection(hcid,fHitCollection); 
}
//______________________________________________________________________________
void He3TargetSD::Clear(){

}
//______________________________________________________________________________
void He3TargetSD::DrawAll(){

}
//______________________________________________________________________________
void He3TargetSD::PrintAll(){

}
