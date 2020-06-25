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
//
/// \file He3TargetDetectorConstruction.hh
/// \brief Definition of the He3TargetDetectorConstruction class

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
#include "G4SubtractionSolid.hh"
#include "G4SystemOfUnits.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;

/// Detector construction class to define materials and geometry.

class He3TargetDetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    He3TargetDetectorConstruction();
    virtual ~He3TargetDetectorConstruction();

    virtual G4VPhysicalVolume* Construct();
    
    G4LogicalVolume* GetScoringVolume() const { return fScoringVolume; }
   
    void SetEXPConfig(int x){ fEXPConfig = x; }  // see He3TargetEnumerate for options 

  protected:
    G4LogicalVolume*  fScoringVolume;

  private:
   int fEXPConfig;  // switch for helmholz coil angle and possibly shield panels  
   bool fDebug,fCheckOverlaps;
   std::vector<partParameters_t> fPartData;
   std::map<G4String,G4Material *> fMaterialsMap;  

   G4Material *GetMaterial(G4String name);

   void BuildBeam(G4LogicalVolume *logicMother); 
   void BuildPickupCoils(G4LogicalVolume *logicMother);    
   void BuildLadderPlate(G4LogicalVolume *logicMother); 
   void BuildShield(int config,G4LogicalVolume *logicMother); 
   void BuildPolarizedHe3(G4LogicalVolume *logicMother);    
   void BuildHelmholtzCoils(int config,const std::string type,G4LogicalVolume *logicMother); 
   void BuildEndWindow(const std::string type,G4LogicalVolume *logicMother);

   G4LogicalVolume *BuildGlassCell();

   int ConstructMaterials(); 
   int ReadData(const char *inpath);
   int SplitString(const char delim,const std::string inStr,std::vector<std::string> &out); 
   int GetPart(const std::string partName,partParameters_t &data);  

};

#endif

