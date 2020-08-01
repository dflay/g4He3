#ifndef He3TargetDetectorConstruction_h
#define He3TargetDetectorConstruction_h 1

#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>
#include <map>  

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

#include "He3TargetPartParameters.hh"
#include "He3TargetEnumerate.hh"
#include "He3TargetSD.hh"

#include "G4Transform3D.hh"
#include "G4PhysicalConstants.hh"
#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4AssemblyVolume.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Torus.hh"
#include "G4Trd.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4UnionSolid.hh"
#include "G4SDManager.hh"
#include "G4SubtractionSolid.hh"
#include "G4SystemOfUnits.hh"
#include "G4String.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;

/// Detector construction class to define materials and geometry.

class He3TargetDetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    He3TargetDetectorConstruction(int config=kSBS_GEN_146);
    virtual ~He3TargetDetectorConstruction();

    virtual G4VPhysicalVolume* Construct();
    
    G4LogicalVolume* GetScoringVolume() const { return fScoringVolume; }
   
    void SetEXPConfig(int x){ fEXPConfig = x; }  // see He3TargetEnumerate for options 

  protected:
    G4LogicalVolume*  fScoringVolume;

  private:
   // data members 
   int fEXPConfig;  // switch for helmholz coil angle and possibly shield panels  
   bool fDebug,fCheckOverlaps;

   // logical volumes 
   G4LogicalVolume **fLogicShield;
   G4LogicalVolume **fLogicEndWindow;
   G4LogicalVolume **fLogicPUCoil;  
   G4LogicalVolume **fLogicPUCoilMount; 
   G4LogicalVolume **fLogicHelmholtzSMaj;    
   G4LogicalVolume **fLogicHelmholtzSMin;    
   G4LogicalVolume **fLogicHelmholtzSRFY;    
   G4LogicalVolume **fLogicHelmholtzMaj;    
   G4LogicalVolume **fLogicHelmholtzMin;    
   G4LogicalVolume **fLogicHelmholtzRFY;  
   G4LogicalVolume *fLogicLadder; 
   G4LogicalVolume *fLogicGlassCell;
   G4LogicalVolume *fLogicHe3;  

   // lists
   std::vector<partParameters_t> fPartData;
   std::map<G4String,G4Material *> fMaterialsMap;  

   // methods 
   G4Material *GetMaterial(G4String name);

   void BuildGlassCell(G4LogicalVolume *logicMother);
   void BuildPolarizedHe3(G4LogicalVolume *logicMother);    
   void BuildBeam(G4LogicalVolume *logicMother); 
   void BuildPickupCoils(G4LogicalVolume *logicMother);    
   void BuildLadderPlate(G4LogicalVolume *logicMother); 
   void BuildShield(int config,G4LogicalVolume *logicMother); 
   void BuildHelmholtzCoils(int config,const std::string type,G4LogicalVolume *logicMother); 
   void BuildEndWindow(const std::string type,G4LogicalVolume *logicMother);

   int ConstructMaterials(); 
   int ReadData(const char *inpath);
   int SplitString(const char delim,const std::string inStr,std::vector<std::string> &out); 
   int GetPart(const std::string partName,partParameters_t &data); 
   int PrintPart(const partParameters_t data,std::string LEN="mm",std::string ANG="deg");  

};

#endif

