#ifndef He3TargetHit_hh 
#define He3TargetHit_hh

// a hit class for the 3He target  

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"
#include "G4Track.hh"
#include "G4VVisManager.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4LogicalVolume.hh"
#include "G4ios.hh"
#include "G4Circle.hh"

class He3TargetHit : public G4VHit { 

   private: 
      G4ThreeVector fPos;      // local coordinates 
      G4ThreeVector fVertex;   // vertex coordinates 
      G4ThreeVector fLabPos;   // global coordinates 
      G4ThreeVector fMom;      // moentum of the particle *before* the step 

      G4double fHitTime;       // time of the hit (relative to start of event)  
      G4double fEdep;          // energy deposition     
      G4double fEnergy;        // initial energy of the particle prior to the hit 
      G4double fLStep;         // spatial length of the step (magnitude) 

      G4int fPID,fMID,fTRID;   // PID of particle causing the hit, mother track ID, track ID 

      G4int fOTrIdx;           // original track info 
      G4int fPTrIdx;           // primary track info 
      G4int fSDTrIdx;          // SD boundary crossing info  

      G4int fVerbosity;        // how verbose is the Print function  

   public:

      He3TargetHit();  
      ~He3TargetHit(); 

      He3TargetHit(const He3TargetHit &rhs); 
      const He3TargetHit& operator=(const He3TargetHit &rhs); 
      G4int operator==(const He3TargetHit &rhs) const;

      inline void *operator new(size_t); 
      inline void operator delete(void *aHit); 
      void *operator new(size_t,void *p){ return p; }  

#ifndef G4NOT_ISO_DELETES
      void operator delete(void *,void *){}
#endif 

      void Draw(); 
      void Print();

      // Setter methods 
      inline void SetPos(G4ThreeVector v)      { fPos     = v; }
      inline void SetVertex(G4ThreeVector v)   { fVertex  = v; } 
      inline void SetLabPos(G4ThreeVector v)   { fLabPos  = v; }   
      inline void SetMomentum(G4ThreeVector v) { fMom     = v; }  

      inline void SetTime(G4double t)          { fHitTime = t; }  
      inline void SetEnergy(G4double e)        { fEnergy  = e; }  
      inline void SetLStep(G4double l)         { fLStep   = l; } 

      inline void SetPID(G4int p)              { fPID     = p; }  
      inline void SetMID(G4int m)              { fMID     = m; }  
      inline void SetTrID(G4int t)             { fTRID    = t; }  

      inline void SetOTrIdx(G4int idx)         { fOTrIdx  = idx; }
      inline void SetPTrIdx(G4int idx)         { fPTrIdx  = idx; }
      inline void SetSDTrIdx(G4int idx)        { fSDTrIdx = idx; }

      inline void SetVerbosity(G4int v)        { fVerbosity = v; } 

      // Getter methods
      inline G4ThreeVector GetPos()            { return fPos;     }  
      inline G4ThreeVector GetVertex()         { return fVertex;  }  
      inline G4ThreeVector GetLabPos()         { return fLabPos;  }  

      inline G4double GetTime()                { return fHitTime; }  
      inline G4double GetEdep()                { return fEdep;    } 

      inline G4int GetPID()                    { return fPID;     }  
      inline G4int GetMID()                    { return fMID;     }  
      inline G4int GetTrID()                   { return fTRID;    } 
 
      inline G4int GetOTrIdx()                 { return fOTrIdx;  }
      inline G4int GetPTrIdx()                 { return fPTrIdx;  }
      inline G4int GetSDTrIdx()                { return fSDTrIdx; }

};

typedef G4THitsCollection<He3TargetHit> He3TargetHitCollection;  

extern G4Allocator<He3TargetHit> He3TargetHitAllocator; 

inline void *He3TargetHit::operator new(size_t) { 
   void *aHit;
   aHit = (void *)He3TargetHitAllocator.MallocSingle(); 
   return aHit;  
} 

inline void He3TargetHit::operator delete(void *aHit) {
   He3TargetHitAllocator.FreeSingle( (He3TargetHit *)aHit ); 
}

#endif 
