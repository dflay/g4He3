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
/// \file He3TargetDetectorConstruction.cc
/// \brief Implementation of the He3TargetDetectorConstruction class

#include "He3TargetDetectorConstruction.hh"

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
#include "G4SystemOfUnits.hh"

//______________________________________________________________________________
He3TargetDetectorConstruction::He3TargetDetectorConstruction()
: G4VUserDetectorConstruction(),
  fScoringVolume(0)
{ }
//______________________________________________________________________________
He3TargetDetectorConstruction::~He3TargetDetectorConstruction()
{ }
//______________________________________________________________________________
G4VPhysicalVolume* He3TargetDetectorConstruction::Construct()
{  
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();
  
   
  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = true;

  //     
  // World
  //
  G4double world_sizeXY = 5*m;
  G4double world_sizeZ  = 5*m;
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");
  
  G4Box* solidWorld =    
    new G4Box("World",                       //its name
       0.5*world_sizeXY, 0.5*world_sizeXY, 0.5*world_sizeZ);     //its size
      
  G4LogicalVolume* logicWorld =                         
    new G4LogicalVolume(solidWorld,          //its solid
                        world_mat,           //its material
                        "World");            //its name
                                   
  G4VPhysicalVolume* physWorld = 
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      logicWorld,            //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking

  // // Envelope parameters
  // //
  // G4double env_sizeXY = 20*cm, env_sizeZ = 30*cm;
  // G4Material* env_mat = nist->FindOrBuildMaterial("G4_WATER");
  //                    
  // //     
  // // Envelope
  // //  
  // G4Box* solidEnv =    
  //   new G4Box("Envelope",                    //its name
  //       0.5*env_sizeXY, 0.5*env_sizeXY, 0.5*env_sizeZ); //its size
  //     
  // G4LogicalVolume* logicEnv =                         
  //   new G4LogicalVolume(solidEnv,            //its solid
  //                       env_mat,             //its material
  //                       "Envelope");         //its name
  //              
  // new G4PVPlacement(0,                       //no rotation
  //                   G4ThreeVector(),         //at (0,0,0)
  //                   logicEnv,                //its logical volume
  //                   "Envelope",              //its name
  //                   logicWorld,              //its mother  volume
  //                   false,                   //no boolean operation
  //                   0,                       //copy number
  //                   checkOverlaps);          //overlaps checking
 
  // //     
  // // Shape 1
  // //  
  // G4Material* shape1_mat = nist->FindOrBuildMaterial("G4_A-150_TISSUE");
  // G4ThreeVector pos1 = G4ThreeVector(0, 2*cm, -7*cm);
  //       
  // // Conical section shape       
  // G4double shape1_rmina =  0.*cm, shape1_rmaxa = 2.*cm;
  // G4double shape1_rminb =  0.*cm, shape1_rmaxb = 4.*cm;
  // G4double shape1_hz = 3.*cm;
  // G4double shape1_phimin = 0.*deg, shape1_phimax = 360.*deg;
  // G4Cons* solidShape1 =    
  //   new G4Cons("Shape1", 
  //   shape1_rmina, shape1_rmaxa, shape1_rminb, shape1_rmaxb, shape1_hz,
  //   shape1_phimin, shape1_phimax);
  //                     
  // G4LogicalVolume* logicShape1 =                         
  //   new G4LogicalVolume(solidShape1,         //its solid
  //                       shape1_mat,          //its material
  //                       "Shape1");           //its name
  //              
  // new G4PVPlacement(0,                       //no rotation
  //                   pos1,                    //at position
  //                   logicShape1,             //its logical volume
  //                   "Shape1",                //its name
  //                   logicEnv,                //its mother  volume
  //                   false,                   //no boolean operation
  //                   0,                       //copy number
  //                   checkOverlaps);          //overlaps checking

  // //     
  // // Shape 2
  // //
  // G4Material* shape2_mat = nist->FindOrBuildMaterial("G4_BONE_COMPACT_ICRU");
  // G4ThreeVector pos2 = G4ThreeVector(0, -1*cm, 7*cm);

  // // Trapezoid shape       
  // G4double shape2_dxa = 12*cm, shape2_dxb = 12*cm;
  // G4double shape2_dya = 10*cm, shape2_dyb = 16*cm;
  // G4double shape2_dz  = 6*cm;      
  // G4Trd* solidShape2 =    
  //   new G4Trd("Shape2",                      //its name
  //             0.5*shape2_dxa, 0.5*shape2_dxb, 
  //             0.5*shape2_dya, 0.5*shape2_dyb, 0.5*shape2_dz); //its size
  //               
  // G4LogicalVolume* logicShape2 =                         
  //   new G4LogicalVolume(solidShape2,         //its solid
  //                       shape2_mat,          //its material
  //                       "Shape2");           //its name
  //              
  // new G4PVPlacement(0,                       //no rotation
  //                   pos2,                    //at position
  //                   logicShape2,             //its logical volume
  //                   "Shape2",                //its name
  //                   logicEnv,                //its mother  volume
  //                   false,                   //no boolean operation
  //                   0,                       //copy number
  //                   checkOverlaps);          //overlaps checking
  //               
  // // Set Shape2 as scoring volume
  // //
  // fScoringVolume = logicShape2;

  // Define elements, compounds here that we'll need 
  // We don't need to supply the molar mass, because G4NistManager does it for us!
  G4int Z,N,ncomponents; 
  G4double abundance;

  G4Isotope *iso_3He = new G4Isotope( "He3", Z=2, N=3 );

  G4Element *el3He = new G4Element("Helium3","3He",ncomponents=1); //Define isotopically pure Helium-3 
  el3He->AddIsotope( iso_3He, abundance=100.0*perCent );

  G4Element *elO  = nist->FindOrBuildElement("O");
  G4Element *elAl = nist->FindOrBuildElement("Al");
  G4Element *elSi = nist->FindOrBuildElement("Si");
  G4Element *elCa = nist->FindOrBuildElement("Ca");
  G4Element *elSr = nist->FindOrBuildElement("Sr");
  G4Element *elBa = nist->FindOrBuildElement("Ba");

  // GE180   
  G4double bigden = 1e9*g/cm3; // why so big? To make these materials not weigh as much in the physics?  
  // gather necessary molecules and compounds  
  // SiO2 60.3%
  G4Material* SiO2 = new G4Material("GE180_SiO2", 2.2*g/cm3, 2 );
  SiO2->AddElement(elSi, 1);
  SiO2->AddElement(elO, 2);
  // fMaterialsMap["GE180_SiO2"] = SiO2;
  // BaO  18.2%
  G4Material* BaO = new G4Material("GE180_BaO", bigden, 2 );
  BaO->AddElement(elBa, 1);
  BaO->AddElement(elO, 1);
  // fMaterialsMap["GE180_BaO"] = BaO;
  // Al2O3 14.3%
  G4Material* Al2O3 = new G4Material("GE180_Al2O3", bigden, 2 );
  Al2O3->AddElement(elAl, 2);
  Al2O3->AddElement(elO, 3);
  // fMaterialsMap["GE180_Al2O3"] = Al2O3;
  // CaO   6.5%
  G4Material* CaO = new G4Material("GE180_CaO", bigden, 2 );
  CaO->AddElement(elCa, 1);
  CaO->AddElement(elO, 1);
  // fMaterialsMap["GE180_CaO"] = CaO;
  // SrO   0.25%
  G4Material* SrO = new G4Material("GE180_SrO", bigden, 2 );
  SrO->AddElement(elSr, 1);
  SrO->AddElement(elO, 1);
  // fMaterialsMap["GE180_SrO"] = SrO;
 
  // Density 2.76 g/cm^3
  // Index of Refraction 1.536
  G4Material* GE180 = new G4Material("GE180", 2.76*g/cm3, 5);
  GE180->AddMaterial(SiO2 , 0.6039);
  GE180->AddMaterial(BaO  , 0.1829);
  GE180->AddMaterial(Al2O3, 0.1439);
  GE180->AddMaterial(CaO  , 0.0659);
  GE180->AddMaterial(SrO  , 0.0034);

  //===================================== END Materials =====================================

  G4double GLASS_WALL_THICKNESS = 0.2*cm;

  // Try a different approach.  Use an assembly volume... 
  G4AssemblyVolume *targetAssembly = new G4AssemblyVolume();

  // use this to set the physical location
  G4ThreeVector P; 

  // use this to rotate the geometry 
  G4RotationMatrix *rm = new G4RotationMatrix(); 
 
  //---- pumping chamber ----
  G4double pcWall       = GLASS_WALL_THICKNESS;            
  G4double pcR_max      = (10.8/2.)*cm;     // 4.25 inches OD  
  G4double pcR_min      = pcR_max - pcWall; // derive inner radius from wall thickness and max radius   
  G4double pcStartTheta = 0.*deg;           // polar angle (relative to z axis) 
  G4double pcStopTheta  = 360.0*deg;        // polar angle
  G4double pcStartPhi   = 0.*deg;           // azimuthal angle (about z axis in xy plane)    
  G4double pcStopPhi    = 360.0*deg;        // azimuthal angle
  G4double r_pc[3]      = {0.*cm,-33.02*cm,0.*cm};  // its position in the logical world 
  G4double rot_pc[3]    = {0.*deg,0.*deg,0.*deg}; 

  G4Sphere *pumpChamberShape        = new G4Sphere("pumpChamberShape",pcR_min,pcR_max,pcStartPhi,pcStopPhi,pcStartTheta,pcStopTheta); 
  G4LogicalVolume *logicPumpChamber = new G4LogicalVolume(pumpChamberShape,GE180,"logicPumpChamber"); 

  // set its position and rotation 
  P.setX(r_pc[0]);        P.setY(r_pc[1]);        P.setZ(r_pc[2]); 
  rm->rotateX(rot_pc[0]); rm->rotateY(rot_pc[1]); rm->rotateZ(rot_pc[2]);  
  // place in the target assembly  
  targetAssembly->AddPlacedVolume(logicPumpChamber,P,rm);  

  //---- end window on target chamber, downstream ----  
  // FIXME: What is the correct material?  Using GE180 for now  
  G4double ewR_max      = 1.05*cm; 
  G4double ewR_min      = ewR_max-GLASS_WALL_THICKNESS; 
  G4double ewStartTheta = 0.*deg;      // polar angle (relative to z axis) 
  G4double ewStopTheta  = 90.0*deg;   // polar angle
  G4double ewStartPhi   = 0.*deg;      // azimuthal angle (about z axis in xy plane)    
  G4double ewStopPhi    = 360.0*deg;   // azimuthal angle
  G4double r_ewd[3]     = {0.*cm,0.*cm,(61.2/2.)*cm}; 
  G4double rot_ewd[3]   = {0.*deg,0.*deg,0.*deg}; 

  // shape 
  G4Sphere *endWindowShapeDn        = new G4Sphere("endWindowShapeDn",ewR_min,ewR_max,ewStartPhi,ewStopPhi,ewStartTheta,ewStopTheta);

  // logical volume 
  G4LogicalVolume *logicEndWindowDn = new G4LogicalVolume(endWindowShapeDn,GE180,"logicEndWindowDn"); 

  // set its position and rotation 
  P.setX(r_ewd[0]);        P.setY(r_ewd[1]);        P.setZ(r_ewd[2]); 
  rm->rotateX(rot_ewd[0]); rm->rotateY(rot_ewd[1]); rm->rotateZ(rot_ewd[2]);  
  // place in the target assembly  
  targetAssembly->AddPlacedVolume(logicEndWindowDn,P,rm);  

  //---- end window on target chamber, upstream ----  
  // FIXME: What is the correct material?  Using GE180 for now  
  ewR_max      = 1.05*cm; 
  ewR_min      = ewR_max-GLASS_WALL_THICKNESS; 
  ewStartTheta = 0.*deg;      // polar angle (relative to z axis) 
  ewStopTheta  = 90.0*deg;    // polar angle
  ewStartPhi   = 0.*deg;      // azimuthal angle (about z axis in xy plane)    
  ewStopPhi    = 360.0*deg;   // azimuthal angle
  G4double r_ewu[3]   = {0.*cm,0.*cm,(-61./2.)*cm};
  G4double rot_ewu[3] = {180.*deg,0.*deg,0.*deg}; 

  // shape 
  G4Sphere *endWindowShapeUp      = new G4Sphere("endWindowShapeUp",ewR_min,ewR_max,ewStartPhi,ewStopPhi,ewStartTheta,ewStopTheta);

  // logical volume 
  G4LogicalVolume *logicEndWindowUp = new G4LogicalVolume(endWindowShapeUp,GE180,"logicEndWindowUp"); 

  // set its position and rotation 
  P.setX(r_ewu[0]);       P.setY(r_ewu[1]);       P.setZ(r_ewu[2]); 
  rm->rotateX(rot_ewu[0]); rm->rotateY(rot_ewu[1]); rm->rotateZ(rot_ewu[2]);  
  // place in the target assembly  
  targetAssembly->AddPlacedVolume(logicEndWindowUp,P,rm);  

  //---- target chamber ----
  // along z axis 

  // shape  
  G4double tcR_max       = 1.05*cm;       
  G4double tcR_min       = tcR_max-GLASS_WALL_THICKNESS; 
  G4double tcLength      = (57.9/2.)*cm; // half length! 
  G4double tcStartAngle  = 0.*deg;  
  G4double tcStopAngle   = 360.*deg;     // full circle  
  G4double r_tc[3]       = {0.*cm,0.*cm,0.*cm};
  G4double rot_tc[3]     = {0.*deg,0.*deg,0.*deg};

  G4Tubs *targetChamberShape = new G4Tubs("targetChamberShape",tcR_min,tcR_max,tcLength,tcStartAngle,tcStopAngle);

  // logical volume 
  G4LogicalVolume *logicTargetChamber = new G4LogicalVolume(targetChamberShape,GE180,"logicTargetChamber");  

  // set its position and rotation 
  P.setX(r_tc[0]);       P.setY(r_tc[1]);       P.setZ(r_tc[2]); 
  rm->rotateX(rot_tc[0]); rm->rotateY(rot_tc[1]); rm->rotateZ(rot_tc[2]);  
  // place in the target assembly  
  targetAssembly->AddPlacedVolume(logicTargetChamber,P,rm);  

  //---- transfer tube elbow, downstream 
  G4double tteR_max     = 0.45*cm;                           // max radius of tube
  G4double tteR_min     = tteR_max - GLASS_WALL_THICKNESS;  // min radius of tube 
  G4double tteRtor      = 0.9*cm;   // major radius of torus   
  G4double tteStartPhi  = 0.*deg;   // probably the span angle about z   
  G4double tteStopPhi   = 90.*deg;   
  G4double r_tted[3]    = {0.*cm,-3.35*cm,25.85*cm}; 
  G4double rot_tted[3]  = {-90.*deg,90.*deg,0.*deg};  
 
  G4Torus *transTubeElDnShape = new G4Torus("transTubeElDnShape",tteR_min,tteR_max,tteRtor,tteStartPhi,tteStopPhi);

  G4LogicalVolume *logicTransTubeElDn = new G4LogicalVolume(transTubeElDnShape,GE180,"logicTransTubeElDn"); 

  // set its position and rotation 
  P.setX(r_tted[0]);        P.setY(r_tted[1]);        P.setZ(r_tted[2]); 
  rm->rotateX(rot_tted[0]); rm->rotateY(rot_tted[1]); rm->rotateZ(rot_tted[2]);  
  // place in the target assembly  
  targetAssembly->AddPlacedVolume(logicTransTubeElDn,P,rm);  

  //---- transfer tube elbow, upstream 
  tteR_max     = 0.45*cm;                           // max radius of tube
  tteR_min     = tteR_max - GLASS_WALL_THICKNESS;  // min radius of tube 
  tteRtor      = 0.9*cm;   // major radius of torus   
  tteStartPhi  = 0.*deg;   // probably the span angle about z   
  tteStopPhi   = 90.*deg;    
  G4double r_tteu[3]    = {0.*cm,-3.35*cm,-25.85*cm}; 
  G4double rot_tteu[3]  = {180.*deg,90.*deg,0.*deg};  

  G4Torus *transTubeElUpShape = new G4Torus("transTubeElUpShape",tteR_min,tteR_max,tteRtor,tteStartPhi,tteStopPhi);

  G4LogicalVolume *logicTransTubeElUp = new G4LogicalVolume(transTubeElUpShape,GE180,"logicTransTubeElUp");

  // set its position and rotation 
  P.setX(r_tteu[0]);        P.setY(r_tteu[1]);        P.setZ(r_tteu[2]); 
  rm->rotateX(rot_tteu[0]); rm->rotateY(rot_tteu[1]); rm->rotateZ(rot_tteu[2]);  
  // place in the target assembly  
  targetAssembly->AddPlacedVolume(logicTransTubeElUp,P,rm);  
 
  //---- transfer tube elbow, downstream lower  
  tteR_max     = 0.45*cm;                           // max radius of tube
  tteR_min     = tteR_max - GLASS_WALL_THICKNESS;  // min radius of tube 
  tteRtor      = 0.9*cm;   // major radius of torus   
  tteStartPhi  = 0.*deg;   // probably the span angle about z   
  tteStopPhi   = 90.*deg;   
  G4double r_ttedl[3]   = {0.*cm,-5.15*cm,1.9*cm};  
  G4double rot_ttedl[3] = {90.*deg,90.*deg,0.*deg};  
 
  G4Torus *transTubeElDnLoShape = new G4Torus("transTubeElDnLoShape",tteR_min,tteR_max,tteRtor,tteStartPhi,tteStopPhi);

  G4LogicalVolume *logicTransTubeElDnLo = new G4LogicalVolume(transTubeElDnLoShape,GE180,"logicTransTubeElDnLo"); 

  // set its position and rotation 
  P.setX(r_ttedl[0]);        P.setY(r_ttedl[1]);        P.setZ(r_ttedl[2]); 
  rm->rotateX(rot_ttedl[0]); rm->rotateY(rot_ttedl[1]); rm->rotateZ(rot_ttedl[2]);  
  // place in the target assembly  
  targetAssembly->AddPlacedVolume(logicTransTubeElDnLo,P,rm);  

  //---- transfer tube elbow, upstream lower  
  tteR_max     = 0.45*cm;                           // max radius of tube
  tteR_min     = tteR_max - GLASS_WALL_THICKNESS;  // min radius of tube 
  tteRtor      = 0.9*cm;   // major radius of torus   
  tteStartPhi  = 0.*deg;   // probably the span angle about z   
  tteStopPhi   = 90.*deg;    
  G4double r_tteul[3]   = {0.*cm,-5.15*cm,-1.9*cm};  
  G4double rot_tteul[3] = {0.*deg,90.*deg,0.*deg};  
 
  G4Torus *transTubeElUpLoShape = new G4Torus("transTubeElUpLoShape",tteR_min,tteR_max,tteRtor,tteStartPhi,tteStopPhi);

  G4LogicalVolume *logicTransTubeElUpLo = new G4LogicalVolume(transTubeElUpLoShape,GE180,"logicTransTubeElUpLo");

  // set its position and rotation 
  G4ThreeVector P_tteul      = G4ThreeVector(r_tteul[0],r_tteul[1],r_tteul[2]);  
  G4RotationMatrix *rm_tteul = new G4RotationMatrix();
  rm_tteul->rotateX(rot_tteul[0]);  
  rm_tteul->rotateY(rot_tteul[1]);  
  rm_tteul->rotateZ(rot_tteul[2]);  
  // place in the target assembly  
  targetAssembly->AddPlacedVolume(logicTransTubeElUpLo,P_tteul,rm_tteul);  

  //---- transfer tube sphere --- 
  G4double ttsR_max      = 1.39*cm; 
  G4double ttsR_min      = ttsR_max - GLASS_WALL_THICKNESS; 
  G4double ttsStartTheta = 0.*deg;      // polar angle (relative to z axis) 
  G4double ttsStopTheta  = 180.0*deg;   // polar angle
  G4double ttsStartPhi   = 0.*deg;      // azimuthal angle (about z axis in xy plane)    
  G4double ttsStopPhi    = 360.0*deg;   // azimuthal angle
  G4double r_tts[3]      = {0.*cm,-8.33*cm,1.0*cm};  
  G4double rot_tts[3]    = {0.*deg,0.*deg,0.*deg};  

  // shape 
  G4Sphere *transTubeSphere       = new G4Sphere("transTubeSphere",ttsR_min,ttsR_max,ttsStartPhi,ttsStopPhi,ttsStartTheta,ttsStopTheta);

  // logical volume 
  G4LogicalVolume *logicTransTubeSphere = new G4LogicalVolume(transTubeSphere,GE180,"logicTransTubeSphere"); 

  // set its position and rotation 
  G4ThreeVector P_tts      = G4ThreeVector(r_tts[0],r_tts[1],r_tts[2]);  
  G4RotationMatrix *rm_tts = new G4RotationMatrix();
  rm_tts->rotateX(rot_tts[0]);  
  rm_tts->rotateY(rot_tts[1]);  
  rm_tts->rotateZ(rot_tts[2]);  
  // place in the target assembly  
  targetAssembly->AddPlacedVolume(logicTransTubeSphere,P_tts,rm_tts);  

  //---- transfer tube: upstream, along z ---

  G4double TT_MAX_RADIUS  = 0.45*cm; 
  G4double TT_MIN_RADIUS  = TT_MAX_RADIUS-GLASS_WALL_THICKNESS; 

  G4double TT_LONG_LENGTH = 23.0*cm; // calculated: 21.51*cm; 
  G4double TT_VERT_LENGTH = 23.0*cm; // calculated: 22.45*cm;

  // shape  
  G4double ttuR_min        = TT_MIN_RADIUS; 
  G4double ttuR_max        = TT_MAX_RADIUS;            
  G4double ttuLength       = TT_LONG_LENGTH/2.; // half length! 
  G4double ttuStartAngle   = 0.*deg;  
  G4double ttuStopAngle    = 360.*deg;    // full circle 
  G4double r_ttuz[3]       = {0.*cm,-4.25*cm,-14.19*cm}; 
  G4double rot_ttuz[3]     = {0.*deg,0.*deg,0.*deg}; 

  G4Tubs *transTubeUpZShape = new G4Tubs("transTubeUpZShape",ttuR_min,ttuR_max,ttuLength,ttuStartAngle,ttuStopAngle);
 
  // logical volume 
  G4LogicalVolume *logicTransTubeUpZ = new G4LogicalVolume(transTubeUpZShape,GE180,"logicTransTubeUpZ"); 

  // set its position and rotation 
  G4ThreeVector P_ttuz      = G4ThreeVector(r_ttuz[0],r_ttuz[1],r_ttuz[2]);  
  G4RotationMatrix *rm_ttuz = new G4RotationMatrix();
  rm_ttuz->rotateX(rot_ttuz[0]);  
  rm_ttuz->rotateY(rot_ttuz[1]);  
  rm_ttuz->rotateZ(rot_ttuz[2]);  
  // place in the target assembly  
  targetAssembly->AddPlacedVolume(logicTransTubeUpZ,P_ttuz,rm_ttuz);  

  //---- transfer tube: upstream, along y ----  

  // shape  
  ttuR_min        = TT_MIN_RADIUS; 
  ttuR_max        = TT_MAX_RADIUS;  
  ttuLength       = TT_VERT_LENGTH/2.; // half length! 
  ttuStartAngle   = 0.*deg;  
  ttuStopAngle    = 360.*deg;    // full circle
  G4double r_ttuy[3]   = {0.*cm,-16.70*cm,-1.0*cm}; 
  G4double rot_ttuy[3] = {90.*deg,0.*deg,0.*deg}; 
 
  G4Tubs *transTubeUpYShape = new G4Tubs("transTubeUpYShape",ttuR_min,ttuR_max,ttuLength,ttuStartAngle,ttuStopAngle);

  // logical volume 
  G4LogicalVolume *logicTransTubeUpY = new G4LogicalVolume(transTubeUpYShape,GE180,"logicTransTubeUpY");  

  // set its position and rotation 
  G4ThreeVector P_ttuy      = G4ThreeVector(r_ttuy[0],r_ttuy[1],r_ttuy[2]);  
  G4RotationMatrix *rm_ttuy = new G4RotationMatrix();
  rm_ttuy->rotateX(rot_ttuy[0]);  
  rm_ttuy->rotateY(rot_ttuy[1]);  
  rm_ttuy->rotateZ(rot_ttuy[2]);  
  // place in the target assembly  
  targetAssembly->AddPlacedVolume(logicTransTubeUpY,P_ttuy,rm_ttuy);  
 
  //---- transfer tube: downstream, along z ---

  // shape  
  G4double ttdR_min        = TT_MIN_RADIUS;  
  G4double ttdR_max        = TT_MAX_RADIUS; 
  G4double ttdLength       = TT_LONG_LENGTH/2.; // half length! 
  G4double ttdStartAngle   = 0.*deg;  
  G4double ttdStopAngle    = 360.*deg;    // full circle
  G4double r_ttdz[3]       = {0.*cm,-4.25*cm,14.19*cm}; 
  G4double rot_ttdz[3]     = {0.*deg,0.*deg,0.*deg}; 
 
  G4Tubs *transTubeDnZShape = new G4Tubs("transTubeDnZShape",ttdR_min,ttdR_max,ttdLength,ttdStartAngle,ttdStopAngle);

  // logical volume 
  G4LogicalVolume *logicTransTubeDnZ = new G4LogicalVolume(transTubeDnZShape,GE180,"logicTransTubeDnZ");  

  // set its position and rotation 
  G4ThreeVector P_ttdz      = G4ThreeVector(r_ttdz[0],r_ttdz[1],r_ttdz[2]);  
  G4RotationMatrix *rm_ttdz = new G4RotationMatrix();
  rm_ttdz->rotateX(rot_ttdz[0]);  
  rm_ttdz->rotateY(rot_ttdz[1]);  
  rm_ttdz->rotateZ(rot_ttdz[2]);  
  // place in the target assembly  
  targetAssembly->AddPlacedVolume(logicTransTubeDnZ,P_ttdz,rm_ttdz);  
 
  //---- transfer tube: upstream, along y ----  
  ttdR_min        = TT_MIN_RADIUS; 
  ttdR_max        = TT_MAX_RADIUS;  
  ttdLength       = TT_VERT_LENGTH/2.; // half length! 
  ttdStartAngle   = 0.*deg;  
  ttdStopAngle    = 360.*deg;    // full circle
  G4double r_ttdy[3]   = {0.*cm,-16.70*cm,1.0*cm}; 
  G4double rot_ttdy[3] = {90.*deg,0.*deg,0.*deg}; 
 
  G4Tubs *transTubeDnYShape = new G4Tubs("transTubeDnYShape",ttdR_min,ttdR_max,ttdLength,ttdStartAngle,ttdStopAngle);

  // logical volume 
  G4LogicalVolume *logicTransTubeDnY = new G4LogicalVolume(transTubeDnYShape,GE180,"logicTransTubeDnY"); 

  // set its position and rotation 
  G4ThreeVector P_ttdy      = G4ThreeVector(r_ttdy[0],r_ttdy[1],r_ttdy[2]);  
  G4RotationMatrix *rm_ttdy = new G4RotationMatrix();
  rm_ttdy->rotateX(rot_ttdy[0]);  
  rm_ttdy->rotateY(rot_ttdy[1]);  
  rm_ttdy->rotateZ(rot_ttdy[2]);  
  // place in the target assembly  
  targetAssembly->AddPlacedVolume(logicTransTubeDnY,P_ttdy,rm_ttdy);  
 
  //---- transfer tube post: downstream, along y ----  

  // shape  
  G4double ttpR_min        = TT_MIN_RADIUS; 
  G4double ttpR_max        = TT_MAX_RADIUS;  
  G4double ttpLength       = 2.3*cm/2.; // half length! 
  G4double ttpStartAngle   = 0.*deg;  
  G4double ttpStopAngle    = 360.*deg;    // full circle 
  G4double r_ttpdy[3]      = {0.*cm,-2.2*cm,26.75*cm}; 
  G4double rot_ttpdy[3]    = {90.*deg,0.*deg,0.*deg}; 

  G4Tubs *transTubePostDnShape = new G4Tubs("transTubePostDnShape",ttpR_min,ttuR_max,ttpLength,ttpStartAngle,ttpStopAngle);
 
  // logical volume 
  G4LogicalVolume *logicTransTubePostDn = new G4LogicalVolume(transTubePostDnShape,GE180,"logicTransTubePostDn");  

  // set its position and rotation 
  G4ThreeVector P_ttpdy      = G4ThreeVector(r_ttpdy[0],r_ttpdy[1],r_ttpdy[2]);  
  G4RotationMatrix *rm_ttpdy = new G4RotationMatrix();
  rm_ttpdy->rotateX(rot_ttpdy[0]);  
  rm_ttpdy->rotateY(rot_ttpdy[1]);  
  rm_ttpdy->rotateZ(rot_ttpdy[2]);  
  // place in the target assembly  
  targetAssembly->AddPlacedVolume(logicTransTubePostDn,P_ttdy,rm_ttpdy);  

  //---- transfer tube post: upstream, along y ----  

  // shape  
  ttpR_min        = TT_MIN_RADIUS; 
  ttpR_max        = TT_MAX_RADIUS;  
  ttpLength       = 2.3*cm/2.; // half length! 
  ttpStartAngle   = 0.*deg;  
  ttpStopAngle    = 360.*deg;    // full circle 
  G4double r_ttpuy[3]   = {0.*cm,-2.2*cm,-26.75*cm}; 
  G4double rot_ttpuy[3] = {90.*deg,0.*deg,0.*deg}; 

  G4Tubs *transTubePostUpShape = new G4Tubs("transTubePostUpShape",ttpR_min,ttuR_max,ttpLength,ttpStartAngle,ttpStopAngle);

  // logical volume 
  G4LogicalVolume *logicTransTubePostUp = new G4LogicalVolume(transTubePostUpShape,GE180,"logicTransTubePostUp");  

  // set its position and rotation 
  G4ThreeVector P_ttpuy      = G4ThreeVector(r_ttpuy[0],r_ttpuy[1],r_ttpuy[2]);  
  G4RotationMatrix *rm_ttpuy = new G4RotationMatrix();
  rm_ttpuy->rotateX(rot_ttpuy[0]);  
  rm_ttpuy->rotateY(rot_ttpuy[1]);  
  rm_ttpuy->rotateZ(rot_ttpuy[2]);  
  // place in the target assembly  
  targetAssembly->AddPlacedVolume(logicTransTubePostUp,P_ttpuy,rm_ttpuy);  

  // now build and orient in mother volume (logicWorld)
  G4ThreeVector PA = G4ThreeVector(0,0,0);  
  G4RotationMatrix *RA = new G4RotationMatrix(); 
  RA->rotateX(0); RA->rotateY(0); RA->rotateZ(0);  
  targetAssembly->MakeImprint(logicWorld,PA,RA,0,checkOverlaps); 

  //---- polarized 3He ----
  
  // G4double gasden = 10.77*atmosphere*(3.016*g/Avogadro)/(300*kelvin*k_Boltzmann);
  // G4Material *pol3He = new G4Material("pol3He", gasden, 1 );
  // pol3He->AddElement(el3He, 1); 

  // // cylinder of polarized 3He 
  // G4double innerRadius = 0.0*cm; 
  // G4double outerRadius = 1.0*cm;
  // G4double length      = (25./2.)*cm; // half length!      
  // G4double startAngle  = 0.*deg;  
  // G4double stopAngle   = 360.*deg;   // full circle 
  // G4Tubs *he3Tube = new G4Tubs("He3_tube",innerRadius,outerRadius,length,startAngle,stopAngle);

  // // logical volume of He3
  // G4LogicalVolume *logicHe3 = new G4LogicalVolume(he3Tube,pol3He,"logicHe3");  
  //   
  // // placement of He3 -- at origin  
  // G4double xhe = 0.*cm; 
  // G4double yhe = 0.*cm; 
  // G4double zhe = 0.*cm; 
  // G4ThreeVector he3pos = G4ThreeVector(xhe,yhe,zhe); 
  // new G4PVPlacement(0,                 // rotation
  //                  he3pos,             // position 
  //                  logicHe3,           // logical volume 
  //                  "physHe3",          // name 
  //                  logicTargetChamber, // logical mother volume  
  //                  false,              // no boolean operations 
  //                  0,                  // copy number 
  //                  checkOverlaps);     // check overlaps

  
  // always return the physical World
  return physWorld;
}

