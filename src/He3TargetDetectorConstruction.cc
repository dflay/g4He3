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

//______________________________________________________________________________
He3TargetDetectorConstruction::He3TargetDetectorConstruction()
: G4VUserDetectorConstruction(),
  fScoringVolume(0),fDebug(false)
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

  // load all part parameters 
  ReadData("./input/He3-parts.csv");   
 
  //---- pumping chamber ----

  partParameters_t pumpCh;
  GetPart("pumpingChamber",pumpCh); 

  G4Sphere *pumpChamberShape = new G4Sphere(pumpCh.name,
                                            pumpCh.r_min     ,pumpCh.r_max,
                                            pumpCh.startPhi  ,pumpCh.stopPhi, 
                                            pumpCh.startTheta,pumpCh.stopTheta); 

  G4ThreeVector P_pc = G4ThreeVector(pumpCh.x,pumpCh.y,pumpCh.z); 
  G4RotationMatrix *rm_pc = new G4RotationMatrix();
  rm_pc->rotateX(pumpCh.rx); 
  rm_pc->rotateY(pumpCh.ry); 
  rm_pc->rotateZ(pumpCh.rz); 
  
  //---- target chamber ----
  // along z axis 

  partParameters_t tgtCh;
  GetPart("targetChamber",tgtCh); 

  G4Tubs *targetChamberShape = new G4Tubs(tgtCh.name,
                                          tgtCh.r_min    ,tgtCh.r_max,
                                          tgtCh.length/2.,
                                          tgtCh.startPhi ,tgtCh.stopPhi); 

  G4ThreeVector P_tc = G4ThreeVector(tgtCh.x,tgtCh.y,tgtCh.z); 
  G4RotationMatrix *rm_tc = new G4RotationMatrix();
  rm_tc->rotateX(tgtCh.rx); 
  rm_tc->rotateY(tgtCh.ry); 
  rm_tc->rotateZ(tgtCh.rz); 

  //---- end window on target chamber, downstream ----  
  // FIXME: What is the correct material?  Using GE180 for now  

  partParameters_t ewDn;
  GetPart("endWindow_dn",ewDn); 

  G4Sphere *endWindowShapeDn = new G4Sphere(ewDn.name,
                                            ewDn.r_min     ,ewDn.r_max,
                                            ewDn.startPhi  ,ewDn.stopPhi, 
                                            ewDn.startTheta,ewDn.stopTheta); 

  G4ThreeVector P_ewd = G4ThreeVector(ewDn.x,ewDn.y,ewDn.z); 
  G4RotationMatrix *rm_ewd = new G4RotationMatrix();
  rm_ewd->rotateX(ewDn.rx); 
  rm_ewd->rotateY(ewDn.ry); 
  rm_ewd->rotateZ(ewDn.rz); 

  //---- end window on target chamber, upstream ----  
  // FIXME: What is the correct material?  Using GE180 for now  

  partParameters_t ewUp;
  GetPart("endWindow_up",ewUp); 

  G4Sphere *endWindowShapeUp = new G4Sphere(ewUp.name,
                                            ewUp.r_min     ,ewUp.r_max,
                                            ewUp.startPhi  ,ewUp.stopPhi, 
                                            ewUp.startTheta,ewUp.stopTheta); 

  G4ThreeVector P_ewu = G4ThreeVector(ewUp.x,ewUp.y,ewUp.z); 
  G4RotationMatrix *rm_ewu = new G4RotationMatrix();
  rm_ewu->rotateX(ewUp.rx); 
  rm_ewu->rotateY(ewUp.ry); 
  rm_ewu->rotateZ(ewUp.rz); 

  //---- transfer tube elbow, downstream 

  partParameters_t tted; 
  GetPart("transTubeEl_dn",tted); 

  G4Torus *transTubeElDnShape = new G4Torus(tted.name,
                                            tted.r_min     ,tted.r_max,tted.r_tor,
                                            tted.startPhi  ,tted.stopPhi); 

  G4ThreeVector P_tted = G4ThreeVector(tted.x,tted.y,tted.z); 
  G4RotationMatrix *rm_tted = new G4RotationMatrix();
  rm_tted->rotateX(tted.rx); 
  rm_tted->rotateY(tted.ry); 
  rm_tted->rotateZ(tted.rz); 

  //---- transfer tube elbow, upstream 
  partParameters_t tteu; 
  GetPart("transTubeEl_up",tteu); 

  G4Torus *transTubeElUpShape = new G4Torus(tteu.name,
                                            tteu.r_min     ,tteu.r_max,tteu.r_tor,
                                            tteu.startPhi  ,tteu.stopPhi); 

  G4ThreeVector P_tteu = G4ThreeVector(tteu.x,tteu.y,tteu.z); 
  G4RotationMatrix *rm_tteu = new G4RotationMatrix();
  rm_tteu->rotateX(tteu.rx); 
  rm_tteu->rotateY(tteu.ry); 
  rm_tteu->rotateZ(tteu.rz); 
 
  //---- transfer tube elbow, downstream lower  
  partParameters_t ttedl; 
  GetPart("transTubeElLo_dn",ttedl); 

  G4Torus *transTubeElDnLoShape = new G4Torus(ttedl.name,
                                              ttedl.r_min     ,ttedl.r_max,ttedl.r_tor,
                                              ttedl.startPhi  ,ttedl.stopPhi); 

  G4ThreeVector P_ttedl = G4ThreeVector(ttedl.x,ttedl.y,ttedl.z); 
  G4RotationMatrix *rm_ttedl = new G4RotationMatrix();
  rm_ttedl->rotateX(ttedl.rx); 
  rm_ttedl->rotateY(ttedl.ry); 
  rm_ttedl->rotateZ(ttedl.rz); 

  //---- transfer tube elbow, upstream lower  
  partParameters_t tteul; 
  GetPart("transTubeElLo_up",tteul); 

  G4Torus *transTubeElUpLoShape = new G4Torus(tteul.name,
                                              tteul.r_min     ,tteul.r_max,tteul.r_tor,
                                              tteul.startPhi  ,tteul.stopPhi); 

  G4ThreeVector P_tteul = G4ThreeVector(tteul.x,tteul.y,tteul.z); 
  G4RotationMatrix *rm_tteul = new G4RotationMatrix();
  rm_tteul->rotateX(tteul.rx); 
  rm_tteul->rotateY(tteul.ry); 
  rm_tteul->rotateZ(tteul.rz); 

  //---- transfer tube sphere --- 
  partParameters_t tts; 
  GetPart("transTubeSphere",tts);

  G4Sphere *transTubeSphere = new G4Sphere(tts.name,
                                           tts.r_min     ,tts.r_max,
                                           tts.startPhi  ,tts.stopPhi,
                                           tts.startTheta,tts.stopTheta); 

  G4ThreeVector P_tts = G4ThreeVector(tts.x,tts.y,tts.z); 
  G4RotationMatrix *rm_tts = new G4RotationMatrix();
  rm_tts->rotateX(tts.rx); 
  rm_tts->rotateY(tts.ry); 
  rm_tts->rotateZ(tts.rz); 

  //---- transfer tube: upstream, along z ---
  partParameters_t ttuz; 
  GetPart("transTubeZ_up",ttuz);

  G4Tubs *transTubeUpZShape = new G4Tubs(ttuz.name,
                                         ttuz.r_min    ,ttuz.r_max,
                                         ttuz.length/2.,
                                         ttuz.startPhi ,ttuz.stopPhi); 

  G4ThreeVector P_ttuz = G4ThreeVector(ttuz.x,ttuz.y,ttuz.z); 
  G4RotationMatrix *rm_ttuz = new G4RotationMatrix();
  rm_ttuz->rotateX(ttuz.rx); 
  rm_ttuz->rotateY(ttuz.ry); 
  rm_ttuz->rotateZ(ttuz.rz); 

  //---- transfer tube: upstream, along y ----  
  partParameters_t ttuy; 
  GetPart("transTubeY_up",ttuy);

  G4Tubs *transTubeUpYShape = new G4Tubs(ttuy.name,
                                         ttuy.r_min    ,ttuy.r_max,
                                         ttuy.length/2.,
                                         ttuy.startPhi ,ttuy.stopPhi); 

  G4ThreeVector P_ttuy = G4ThreeVector(ttuy.x,ttuy.y,ttuy.z); 
  G4RotationMatrix *rm_ttuy = new G4RotationMatrix();
  rm_ttuy->rotateX(ttuy.rx); 
  rm_ttuy->rotateY(ttuy.ry); 
  rm_ttuy->rotateZ(ttuy.rz); 
 
  //---- transfer tube: downstream, along z ---
  partParameters_t ttdz; 
  GetPart("transTubeZ_dn",ttdz);

  G4Tubs *transTubeDnZShape = new G4Tubs(ttdz.name,
                                         ttdz.r_min    ,ttdz.r_max,
                                         ttdz.length/2.,
                                         ttdz.startPhi ,ttdz.stopPhi); 

  G4ThreeVector P_ttdz = G4ThreeVector(ttdz.x,ttdz.y,ttdz.z); 
  G4RotationMatrix *rm_ttdz = new G4RotationMatrix();
  rm_ttdz->rotateX(ttdz.rx); 
  rm_ttdz->rotateY(ttdz.ry); 
  rm_ttdz->rotateZ(ttdz.rz); 
 
  //---- transfer tube: downstream, along y ----  
  // this has two components since it joins with a sphere 
  // - short (above sphere) 
  // - long  (below sphere)  

  // long (below) 
  partParameters_t ttdby; 
  GetPart("transTubeYB_dn",ttdby);

  G4Tubs *transTubeDnBYShape = new G4Tubs(ttdby.name,
                                          ttdby.r_min    ,ttdby.r_max,
                                          ttdby.length/2.,
                                          ttdby.startPhi ,ttdby.stopPhi); 

  G4ThreeVector P_ttdby = G4ThreeVector(ttdby.x,ttdby.y,ttdby.z); 
  G4RotationMatrix *rm_ttdby = new G4RotationMatrix();
  rm_ttdby->rotateX(ttdby.rx); 
  rm_ttdby->rotateY(ttdby.ry); 
  rm_ttdby->rotateZ(ttdby.rz); 

  // short (above) 
  partParameters_t ttday; 
  GetPart("transTubeYA_dn",ttday);

  G4Tubs *transTubeDnAYShape = new G4Tubs(ttday.name,
                                          ttday.r_min    ,ttday.r_max,
                                          ttday.length/2.,
                                          ttday.startPhi ,ttday.stopPhi); 

  G4ThreeVector P_ttday = G4ThreeVector(ttday.x,ttday.y,ttday.z); 
  G4RotationMatrix *rm_ttday = new G4RotationMatrix();
  rm_ttday->rotateX(ttday.rx); 
  rm_ttday->rotateY(ttday.ry); 
  rm_ttday->rotateZ(ttday.rz); 
 
  //---- transfer tube post: downstream, along y ---- 

  // shape  
  // G4double ttpR_min        = TT_MIN_RADIUS; 
  // G4double ttpR_max        = TT_MAX_RADIUS;  
  // G4double ttpLength       = 2.3*cm/2.; // half length! 
  // G4double ttpStartAngle   = 0.*deg;  
  // G4double ttpStopAngle    = 360.*deg;    // full circle 
  // G4double r_ttpdy[3]      = {0.*cm,-2.2*cm,26.75*cm}; 
  // G4double rot_ttpdy[3]    = {90.*deg,0.*deg,0.*deg}; 

  // G4Tubs *transTubePostDnShape = new G4Tubs("transTubePostDnShape",ttpR_min,ttuR_max,ttpLength,ttpStartAngle,ttpStopAngle);

  // // set its position and rotation 
  // G4ThreeVector P_ttpdy      = G4ThreeVector(r_ttpdy[0],r_ttpdy[1],r_ttpdy[2]);  
  // G4RotationMatrix *rm_ttpdy = new G4RotationMatrix();
  // rm_ttpdy->rotateX(rot_ttpdy[0]);  
  // rm_ttpdy->rotateY(rot_ttpdy[1]);  
  // rm_ttpdy->rotateZ(rot_ttpdy[2]);  

  partParameters_t ttpdy; 
  GetPart("transTubePost_dn",ttpdy);

  G4Tubs *transTubePostDnShape = new G4Tubs(ttpdy.name,
                                            ttpdy.r_min    ,ttpdy.r_max,
                                            ttpdy.length/2.,
                                            ttpdy.startPhi ,ttpdy.stopPhi); 

  G4ThreeVector P_ttpdy = G4ThreeVector(ttpdy.x,ttpdy.y,ttpdy.z); 
  G4RotationMatrix *rm_ttpdy = new G4RotationMatrix();
  rm_ttpdy->rotateX(ttpdy.rx); 
  rm_ttpdy->rotateY(ttpdy.ry); 
  rm_ttpdy->rotateZ(ttpdy.rz); 

  //---- transfer tube post: upstream, along y ----  
  partParameters_t ttpuy; 
  GetPart("transTubePost_up",ttpuy);

  G4Tubs *transTubePostUpShape = new G4Tubs(ttpuy.name,
                                            ttpuy.r_min    ,ttpuy.r_max,
                                            ttpuy.length/2.,
                                            ttpuy.startPhi ,ttpuy.stopPhi); 

  G4ThreeVector P_ttpuy = G4ThreeVector(ttpuy.x,ttpuy.y,ttpuy.z); 
  G4RotationMatrix *rm_ttpuy = new G4RotationMatrix();
  rm_ttpuy->rotateX(ttpuy.rx); 
  rm_ttpuy->rotateY(ttpuy.ry); 
  rm_ttpuy->rotateZ(ttpuy.rz); 

  // Create unions to make this a single continuous piece of material  
  G4UnionSolid *glassCell; // this is the top level object everything becomes 

  // build the target chamber + endcaps
  // target chamber + upstream window 
  glassCell = new G4UnionSolid("gc_tc_ewu",targetChamberShape,endWindowShapeUp,rm_ewu,P_ewu);
  // downstream window  
  glassCell = new G4UnionSolid("gc_tc_ewud",glassCell,endWindowShapeDn,rm_ewd,P_ewd);  

  // transfer tube posts
  // upstream  
  glassCell = new G4UnionSolid("gc_tc_ewud_pu" ,glassCell,transTubePostUpShape,rm_ttpuy,P_ttpuy);  
  // downstream 
  glassCell = new G4UnionSolid("gc_tc_ewud_pud",glassCell,transTubePostDnShape,rm_ttpdy,P_ttpdy); 

  // transfer tube elbows 
  // upstream
  glassCell = new G4UnionSolid("gc_tc_ewud_pud_eu" ,glassCell,transTubeElUpShape,rm_tteu,P_tteu);  
  // downstream  
  glassCell = new G4UnionSolid("gc_tc_ewud_pud_eud",glassCell,transTubeElDnShape,rm_tted,P_tted); 

  // transfer tubes along z 
  // upstream
  glassCell = new G4UnionSolid("gc_tc_ewud_pud_eud_tzu" ,glassCell,transTubeUpZShape,rm_ttuz,P_ttuz); 
  // downstream  
  glassCell = new G4UnionSolid("gc_tc_ewud_pud_eud_tzud",glassCell,transTubeDnZShape,rm_ttdz,P_ttdz); 

  // transfer tube elbows [lower]  
  // upstream
  glassCell = new G4UnionSolid("gc_tc_ewud_pud_eud_tzud_elu" ,glassCell,transTubeElUpLoShape,rm_tteul,P_tteul);  
  // downstream  
  glassCell = new G4UnionSolid("gc_tc_ewud_pud_eud_tzud_elud",glassCell,transTubeElDnLoShape,rm_ttedl,P_ttedl); 

  // transfer tubes along y 
  // upstream
  glassCell = new G4UnionSolid("gc_tc_ewud_pud_eud_tzud_tyu" ,glassCell,transTubeUpYShape,rm_ttuy,P_ttuy); 
  // downstream: above  
  glassCell = new G4UnionSolid("gc_tc_ewud_pud_eud_tzud_tyua",glassCell,transTubeDnAYShape,rm_ttday,P_ttday); 
  // downstream: sphere   
  glassCell = new G4UnionSolid("gc_tc_ewud_pud_eud_tzud_tyuas",glassCell,transTubeSphere,rm_tts,P_tts); 
  // downstream: below  
  glassCell = new G4UnionSolid("gc_tc_ewud_pud_eud_tzud_tyuasb",glassCell,transTubeDnBYShape,rm_ttdby,P_ttdby);

  // pumping chamber.  also change the name since everything is connected now 
  glassCell = new G4UnionSolid("glassCell",glassCell,pumpChamberShape,rm_pc,P_pc);

  G4LogicalVolume *logicGlassCell = new G4LogicalVolume(glassCell,GE180,"logicGlassCell");

  G4VisAttributes *visGC = new G4VisAttributes(); 
  visGC->SetColour( G4Colour::White() ); 

  logicGlassCell->SetVisAttributes(visGC); 

  G4ThreeVector P_tgt_o = G4ThreeVector(0,20.*cm,0); 
  new G4PVPlacement(0,P_tgt_o,logicGlassCell,"physGC",logicWorld,false,0,checkOverlaps);       

  //---- polarized 3He ----
  // G4double gasden = 10.77*atmosphere*(3.016*g/Avogadro)/(300*kelvin*k_Boltzmann);
  // G4Material *pol3He = new G4Material("pol3He", gasden, 1 );
  // pol3He->AddElement(el3He, 1); 

  // // cylinder of polarized 3He 
  // G4double innerRadius = 0.0*cm; 
  // G4double outerRadius = tgtCh.r_min;     // fill the target chamber  
  // G4double length      = tgtCh.length/2.; // half-length of target chamber       
  // G4double startAngle  = 0.*deg;  
  // G4double stopAngle   = 360.*deg;   // full circle 
  // G4Tubs *he3Tube = new G4Tubs("He3_tube",innerRadius,outerRadius,length,startAngle,stopAngle);

  // // logical volume of He3
  // G4LogicalVolume *logicHe3 = new G4LogicalVolume(he3Tube,pol3He,"logicHe3");  
  //  
  // // set the color of He3 
  // G4VisAttributes *He3VisAtt = new G4VisAttributes(); 
  // He3VisAtt->SetColour( G4Colour::Yellow() );

  // logicHe3->SetVisAttributes(He3VisAtt);  
 
  // // placement of He3 -- at origin  
  // G4ThreeVector he3pos = G4ThreeVector(0,0,0); 
  // new G4PVPlacement(0,                 // rotation
  //                  he3pos,             // position 
  //                  logicHe3,           // logical volume 
  //                  "physHe3",          // name 
  //                  logicGlassCell,     // logical mother volume is the target chamber 
  //                  false,              // no boolean operations 
  //                  0,                  // copy number 
  //                  checkOverlaps);     // check overlaps

  // always return the physical World
  return physWorld;
}
//______________________________________________________________________________
int He3TargetDetectorConstruction::GetPart(const char *partName,partParameters_t &data){
   // get the part data based on name
   bool found=false; 
   const int N = fPartData.size();
   for(int i=0;i<N;i++){
      if( fPartData[i].name.compare(partName)==0 ){
	 data = fPartData[i];
	 found = true;
      } 
   }

   if(found && fDebug){
      std::cout << "[He3TargetDetectorConstruction::GetPart]: FOUND part: " << std::endl;
      std::cout << " name: " << data.name << " shape: " << data.shape << std::endl;
      std::cout << " len unit: " << data.len_unit << " ang unit: " << data.ang_unit << std::endl;
      std::cout << " major radius (torus): " << data.r_tor << std::endl;
      std::cout << " min radius: " << data.r_min << " max radius: " << data.r_max << std::endl; 
      std::cout << " length: " << data.length << std::endl; 
      std::cout << " THETA start theta: " << data.startTheta << " stop theta: " << data.stopTheta << std::endl; 
      std::cout << " PHI start phi: " << data.startPhi   << " stop phi: "   << data.stopPhi << std::endl; 
      std::cout << " POS x: " << data.x << " y: " << data.y << " z: " << data.z << std::endl; 
      std::cout << " ROT rx: " << data.rx << " ry: " << data.ry << " rz: " << data.rz << std::endl; 
      std::cout << "------------------------------------------------------" << std::endl; 
   }

   if(!found){
      std::cout << "[He3TargetDetectorConstruction::GetPart]: ERROR! Cannot find part " 
                << partName << std::endl;
      exit(1);
   }
   return 0; 
}
//______________________________________________________________________________
int He3TargetDetectorConstruction::ReadData(const char *inpath){

   std::string aLine;
   std::vector<std::string> row;

   std::ifstream infile;
   infile.open(inpath);

   int k=0;

   if( infile.fail() ){
      G4cout << "[He3TargetDetectorConstruction::ReadData]: Cannot open the file: " << inpath << G4endl;
      return 1;
   }else{
      G4cout << "[He3TargetDetectorConstruction::ReadData]: Opened the file: " << inpath << G4endl;
      while( !infile.eof() ){
         std::getline(infile,aLine);
         if(k>0) row.push_back(aLine); // skip the header row 
	 k++;
      }
      // row.pop_back();  // since we skip the first row, don't need to worry about this 
      infile.close();
   }

   int NROW = row.size();
   std::vector<std::string> col;

   partParameters_t dataPt; 
   double r_min=0,r_max=0,r_tor=0;
   double length=0,startTheta=0,stopTheta=0,startPhi=0,stopPhi=0;
   double x=0,y=0,z=0,rx=0,ry=0,rz=0;
   double LEN_UNIT=1.;
   double DEG_UNIT=1.;

   // now parse the data
   int rc=0;
   for(int i=0;i<NROW;i++){
      // split the row into a vector which represents the columns 
      rc   = SplitString(',',row[i],col);
      if(rc!=0){
         G4cout << "[He3TargetDetectorConstruction::ReadData]: Cannot parse string " << row[i] << G4endl;
	 return 1;
      }
      // fill the data vector
      dataPt.name     = col[0];  
      dataPt.shape    = col[1];  
      dataPt.len_unit = col[2];  
      dataPt.ang_unit = col[3];  
      r_tor           = std::atof( col[4].c_str()  );  
      r_min           = std::atof( col[5].c_str()  );  
      r_max           = std::atof( col[6].c_str()  );  
      length          = std::atof( col[7].c_str()  );  
      startTheta      = std::atof( col[8].c_str()  );  
      stopTheta       = std::atof( col[9].c_str()  );  
      startPhi        = std::atof( col[10].c_str() );  
      stopPhi         = std::atof( col[11].c_str() );  
      x               = std::atof( col[12].c_str() );  
      y               = std::atof( col[13].c_str() );  
      z               = std::atof( col[14].c_str() );  
      rx              = std::atof( col[15].c_str() );  
      ry              = std::atof( col[16].c_str() );  
      rz              = std::atof( col[17].c_str() );
      // convert to units
      if( dataPt.len_unit.compare("mm")==0 ) LEN_UNIT = mm;
      if( dataPt.len_unit.compare("cm")==0)  LEN_UNIT = cm; 
      if( dataPt.len_unit.compare("m")==0)   LEN_UNIT = m; 
      if( dataPt.len_unit.compare("in")==0)  LEN_UNIT = cm/2.54; 
      if( dataPt.ang_unit.compare("deg")==0 ) DEG_UNIT = deg;
      if( dataPt.ang_unit.compare("rad")==0)  DEG_UNIT = rad; 
      // store data 
      dataPt.r_tor      = LEN_UNIT*r_tor;  
      dataPt.r_min      = LEN_UNIT*r_min;  
      dataPt.r_max      = LEN_UNIT*r_max;  
      dataPt.length     = LEN_UNIT*length;  
      dataPt.startTheta = DEG_UNIT*startTheta;  
      dataPt.stopTheta  = DEG_UNIT*stopTheta;  
      dataPt.startPhi   = DEG_UNIT*startPhi;  
      dataPt.stopPhi    = DEG_UNIT*stopPhi;  
      dataPt.x          = LEN_UNIT*x;  
      dataPt.y          = LEN_UNIT*y;  
      dataPt.z          = LEN_UNIT*z;  
      dataPt.rx         = DEG_UNIT*rx;  
      dataPt.ry         = DEG_UNIT*ry;  
      dataPt.rz         = DEG_UNIT*rz;
      if(i==0) G4cout << "Loading 3He target part parameters: " << G4endl;
      G4cout << "name: "                  << dataPt.name       << " shape: "      << dataPt.shape 
	     << " major radius (torus): " << dataPt.r_tor 
	     << " min radius: "           << dataPt.r_min      << " max radius: " << dataPt.r_max 
	     << " length: "               << dataPt.length 
	     << " start theta: "          << dataPt.startTheta << " stop theta: " << dataPt.stopTheta 
	     << " start phi: "            << dataPt.startPhi   << " stop phi: "   << dataPt.stopPhi 
	     << " x: "                    << dataPt.x          << " y: "          << dataPt.y 
	     << " z: "                    << dataPt.z 
	     << " rx: "                   << dataPt.rx         << " ry: " << dataPt.ry 
	     << " rz: "                   << dataPt.rz << G4endl;  
      fPartData.push_back(dataPt);  
      // clean up
      col.clear();
   }
   return 0;
}
//______________________________________________________________________________
int He3TargetDetectorConstruction::SplitString(const char delim,const std::string inStr,std::vector<std::string> &out){
   // split a string by a delimiter
   std::stringstream ss(inStr);
   while( ss.good() ){
      std::string substr;
      std::getline(ss,substr,delim);
      out.push_back(substr);
   }
   return 0;
}
