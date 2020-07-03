#include "He3TargetHit.hh"
//______________________________________________________________________________
G4Allocator<He3TargetHit> He3TargetHitAllocator;
//______________________________________________________________________________
He3TargetHit::He3TargetHit(){
  fPos       = G4ThreeVector(0,0,0);
  fVertex    = G4ThreeVector(0,0,0);
  fLabPos    = G4ThreeVector(0,0,0);
  fMom       = G4ThreeVector(0,0,0);
  fEnergy    = 0;
  fEdep      = 0.0;
  fMID       = -1;
  fPID       = -1e9;
  fVerbosity = 0; 
}
//______________________________________________________________________________
He3TargetHit::~He3TargetHit(){

}
//______________________________________________________________________________
He3TargetHit::He3TargetHit(const He3TargetHit &rhs)
  : G4VHit()
{ 
   // copy constructor 
   fPos    = rhs.fPos;
   fVertex = rhs.fVertex;
   fLabPos = rhs.fLabPos;
   fMom    = rhs.fMom;

   fHitTime = rhs.fHitTime;
   fEdep    = rhs.fEdep;
   fEnergy  = rhs.fEnergy;

   fPID     = rhs.fPID;
   fMID     = rhs.fMID;
   fTRID    = rhs.fTRID;

   //Original track info:
   fOTrIdx = rhs.fOTrIdx;

   //Primary track info:
   fPTrIdx = rhs.fPTrIdx;

   //SD boundary crossing track info:
   fSDTrIdx = rhs.fSDTrIdx;

}
//______________________________________________________________________________
const He3TargetHit& He3TargetHit::operator=(const He3TargetHit &rhs){
   // assignment operator
   fPos    = rhs.fPos;
   fVertex = rhs.fVertex;
   fLabPos = rhs.fLabPos;
   fMom    = rhs.fMom;

   fHitTime = rhs.fHitTime;
   fEdep    = rhs.fEdep;
   fEnergy  = rhs.fEnergy;

   fPID     = rhs.fPID;
   fMID     = rhs.fMID;
   fTRID    = rhs.fTRID;

   //Original track info:
   fOTrIdx = rhs.fOTrIdx;

   //Primary track info:
   fPTrIdx = rhs.fPTrIdx;

   //SD boundary crossing track info:
   fSDTrIdx = rhs.fSDTrIdx;

   return *this;

}
//______________________________________________________________________________
G4int He3TargetHit::operator==(const He3TargetHit &rhs) const{ 
   return (this==&rhs) ? 1 : 0; 
}
//______________________________________________________________________________
void He3TargetHit::Draw(){
   G4VVisManager *pVVisManager = G4VVisManager::GetConcreteInstance(); 
   if( pVVisManager ){ 
      G4Circle circle(fLabPos); 
      circle.SetScreenSize(5); 
      circle.SetFillStyle( G4Circle::filled ); 
      G4VisAttributes attr( G4Colour::Red() );
      circle.SetVisAttributes(attr); 
      pVVisManager->Draw(circle);  
   }
}
//______________________________________________________________________________
void He3TargetHit::Print(){
   if(fVerbosity>0){
      std::cout << "[He3TargetHit]: PID = " << fPID << std::endl;

      std::cout << "[He3TargetHit]: Pos (x,y,z) = (" << fPos.getX()/cm << " cm, " 
	                                             << fPos.getY()/cm << " cm, " 
	                                             << fPos.getZ()/cm << " cm)" << std::endl; 

      std::cout << "[He3TargetHit]: Vertex (x,y,z) = (" << fVertex.getX()/cm << " cm, " 
	                                                << fVertex.getY()/cm << " cm, " 
	                                                << fVertex.getZ()/cm << " cm)" << std::endl;

      std::cout << "[He3TargetHit]: Lab Pos (x,y,z) = (" << fLabPos.getX()/cm << " cm, " 
	                                                 << fLabPos.getY()/cm << " cm, " 
	                                                 << fLabPos.getZ()/cm << " cm)" << std::endl;

      std::cout << "[He3TargetHit]: Momentum (x,y,z) = (" << fMom.getX()/MeV << " MeV, " 
	                                                  << fMom.getY()/MeV << " MeV, " 
	                                                  << fMom.getZ()/MeV << " MeV)" << std::endl;

      std::cout << "[He3TargetHit]: Energy (init) = " << fEnergy/MeV << " MeV" << std::endl;
      std::cout << "[He3TargetHit]: Edep          = " << fEdep/MeV   << " MeV" << std::endl;
   }
}
