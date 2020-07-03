#ifndef He3TargetSD_hh
#define He3TargetSD_hh

// Sensitive detector for 3He
#include <cstdlib>
#include <vector>  

#include "He3TargetHit.hh"

#include "G4Step.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VSensitiveDetector.hh"
#include "G4LogicalVolume.hh"
#include "G4Track.hh"
#include "G4ParticleDefinition.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4Box.hh"
#include "G4SystemOfUnits.hh"
#include "G4String.hh"

class G4Step; 
class G4HCofThisEvent; 
class G4TouchableHistory; 

class He3TargetSD: public G4VSensitiveDetector { 

   private: 
      He3TargetHitCollection *fHitCollection;

      // may need these later... 
      G4int fNTimeBins;          // number of time bins to record energy deposition after the start of the hit 
      G4double fHitTimeWindow;   // maximum time *after* the first tracking step of the hit efore starting a new hit
      G4double fEnergyThreshold; // threshold on the minimum (summed) energy deposition of a hit to record to the output     

   public:
      He3TargetSD(G4String name,G4String hcName); 
      ~He3TargetSD(); 

      // methods from the base class that we need to define 
      virtual void Initialize(G4HCofThisEvent *hce); 
      virtual void EndOfEvent(G4HCofThisEvent *hce); 
      virtual G4bool ProcessHits(G4Step *step,G4TouchableHistory *hist);  

      // our own methods  
      void Clear();
      void DrawAll();
      void PrintAll();
  
      inline void SetNTimeBins(G4int N)          { fNTimeBins       = N; }  
      inline void SetHitTimeWindow(G4double t)   { fHitTimeWindow   = t; }  
      inline void SetEnergyThreshold(G4double e) { fEnergyThreshold = e; } 

}; 

#endif 
