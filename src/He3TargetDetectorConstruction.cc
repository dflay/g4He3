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
  fScoringVolume(0),fEXPConfig(kSBS_GEN_368),fDebug(false),
  fCheckOverlaps(true)
{ }
//______________________________________________________________________________
He3TargetDetectorConstruction::~He3TargetDetectorConstruction()
{ 
   fPartData.clear(); 
}
//______________________________________________________________________________
G4VPhysicalVolume* He3TargetDetectorConstruction::Construct()
{  

  // build materials 
  if( fMaterialsMap.empty() ) ConstructMaterials();

  // load all part parameters 
  ReadData("./input/He3-parts.csv");   

  // Option to switch on/off checking of volumes overlaps
  // G4bool fCheckOverlaps = true;

  // World
  G4double world_sizeXY = 5*m;
  G4double world_sizeZ  = 5*m;
  G4Material* world_mat = GetMaterial("G4air"); 
  
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
                      fCheckOverlaps);        //overlaps checking

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
  //                   fCheckOverlaps);          //overlaps checking
 
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
  //                   fCheckOverlaps);          //overlaps checking

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
  //                   fCheckOverlaps);          //overlaps checking
  //               
  // // Set Shape2 as scoring volume
  // //
  // fScoringVolume = logicShape2;

  //---- Cu end windows for the target cell 
  // - main shaft [tube] 
  // - lip [tube]  
  // - rounded lip [hemisphere/tube] 
  // - endcap [hemisphere] 
 
  BuildEndWindow("upstream"  ,logicWorld); 
  BuildEndWindow("downstream",logicWorld); 

  //---- glass cell ---- 
  // we return the logical volume since we want to place the 3He inside of it.  
  // could consider building the 3He material inside the BuildGlassCell function  
  G4LogicalVolume *logicGlassCell = BuildGlassCell();

  G4VisAttributes *visGC = new G4VisAttributes(); 
  visGC->SetColour( G4Colour::White() );
  visGC->SetForceWireframe(true);  
  logicGlassCell->SetVisAttributes(visGC); 

  // place the volume
  // - note that this is relative to the *target chamber* as that is the first object in the union 
  // - rotation puts the cell oriented such that the pumping chamber is vertically above
  //   and the beam enters from the side where the small sphere on the transfer tube is 
  //   closest to the upstream side 
  G4ThreeVector P_tgt_o = G4ThreeVector(0.*cm,0.*cm,0.*cm);   
  G4RotationMatrix *rm_gc = new G4RotationMatrix(); 
  rm_gc->rotateX(0.*deg);  
  rm_gc->rotateY(180.*deg);  
  rm_gc->rotateZ(180.*deg); 
 
  new G4PVPlacement(rm_gc,P_tgt_o,logicGlassCell,"physGC",logicWorld,false,0,fCheckOverlaps);       

  // cylinder of polarized 3He
  BuildPolarizedHe3(logicGlassCell);

  // helmholtz coils
  fEXPConfig = kSBS_GEN_368;  
  BuildHelmholtzCoils(fEXPConfig,"maj",logicWorld);  
  BuildHelmholtzCoils(fEXPConfig,"rfy",logicWorld);  
  BuildHelmholtzCoils(fEXPConfig,"min",logicWorld); 

  // shield 
  BuildShield(fEXPConfig,logicWorld); 

  // support ladder 
  BuildLadderPlate(logicWorld);

  // pickup coils 
  BuildPickupCoils(logicWorld); 

  // beam path (for reference only) 
  BuildBeam(logicWorld);  

  // always return the physical World
  return physWorld;
}
//______________________________________________________________________________
void He3TargetDetectorConstruction::BuildBeam(G4LogicalVolume *logicMother){
   // build the nominal trajectory of the beam to use as a reference

   partParameters_t beam; 
   beam.r_min  = 0.*mm; 
   beam.r_max  = 2.*mm; 
   beam.length = 5.*m;
   beam.dPhi   = 360.*deg;

   G4Tubs *beamShape = new G4Tubs("beam",
                                  beam.r_min    ,beam.r_max,
                                  beam.length/2., 
                                  beam.startPhi ,beam.dPhi); 
  
   G4VisAttributes *vis = new G4VisAttributes();
   vis->SetColour( G4Colour::Blue() ); 

   G4LogicalVolume *logicBeam = new G4LogicalVolume(beamShape,GetMaterial("G4air"),"logicBeam");
   logicBeam->SetVisAttributes(vis); 

   new G4PVPlacement(0,                     // rotation relative to mother         
                     G4ThreeVector(0,0,0),  // position relative to mother           
                     logicBeam,             // logical object     
                     "physBeam",            // name of physical placement     
                     logicMother,           // logical mother       
                     false,                 // boolean object? (true or false)    
                     0,                     // copy no   
                     fCheckOverlaps);       // check overlaps       

}
//______________________________________________________________________________
void He3TargetDetectorConstruction::BuildPolarizedHe3(G4LogicalVolume *logicMother){
  // build the polarized 3He logical volume and place it in the logicMother 
  // FIXME: - Check volume overlap with glass cell 
  //        - Check volume overlap on endcap; seems like rlip could use a + 0.5*mm offset 

  //---- target chamber component  
  partParameters_t tc; 
  GetPart("targetChamber",tc); 
  tc.r_max = tc.r_min; 
  tc.r_min = 0.*cm;
 
  G4Tubs *tcShape = new G4Tubs("He3_tc",
                               tc.r_min    ,tc.r_max,
                               tc.length/2.,
                               tc.startPhi ,tc.dPhi);

  //---- end window (upstream, new)

  // build the main shaft 
  partParameters_t mshu;
  GetPart("ew_mainShaft_up",mshu);
  mshu.r_max = mshu.r_min;   
  mshu.r_min = 0.*cm;  
  
  G4double z0 = 0.5*(tc.length + mshu.length);  
  mshu.z     = -z0;  

  G4Tubs *mainShaft_up = new G4Tubs(mshu.name,
	                            mshu.r_min    ,mshu.r_max,
	                            mshu.length/2.,
	                            mshu.startPhi ,mshu.dPhi);

  G4ThreeVector P_mshu      = G4ThreeVector(mshu.x,mshu.y,mshu.z);  
  G4RotationMatrix *rm_mshu = new G4RotationMatrix();
  rm_mshu->rotateX(mshu.rx);   
  rm_mshu->rotateY(mshu.ry);   
  rm_mshu->rotateZ(mshu.rz);   

  // lip  
  partParameters_t lipu;
  GetPart("ew_lip_up",lipu);
  lipu.r_max = lipu.r_min; 
  lipu.r_min = 0.*cm;  
  lipu.z     = -z0 - lipu.z;  

  G4Tubs *lip_up = new G4Tubs(lipu.name,
	                      lipu.r_min    ,lipu.r_max,
	                      lipu.length/2.,
	                      lipu.startPhi ,lipu.dPhi);

  G4ThreeVector P_lipu      = G4ThreeVector(lipu.x,lipu.y,lipu.z);  
  G4RotationMatrix *rm_lipu = new G4RotationMatrix();
  rm_lipu->rotateX(lipu.rx);   
  rm_lipu->rotateY(lipu.ry);   
  rm_lipu->rotateZ(lipu.rz);   

  // rounded lip 
  partParameters_t rlipu;
  GetPart("ew_rlip_up",rlipu); 
  rlipu.r_max = rlipu.r_min; 
  rlipu.r_min = 0.*cm; 
  rlipu.z     = -z0 - rlipu.z; // subtract an additional 0.5*mm?  
  rlipu.ry    = 180.*deg;  

  G4Sphere *roundLip_up = new G4Sphere(rlipu.name,
	                               rlipu.r_min     ,rlipu.r_max,
	                               rlipu.startPhi  ,rlipu.dPhi,
	                               rlipu.startTheta,rlipu.dTheta);

  G4ThreeVector P_rlipu      = G4ThreeVector(rlipu.x,rlipu.y,rlipu.z);  
  G4RotationMatrix *rm_rlipu = new G4RotationMatrix();
  rm_rlipu->rotateX(rlipu.rx);   
  rm_rlipu->rotateY(rlipu.ry);   
  rm_rlipu->rotateZ(rlipu.rz);   

  // endcap 
  partParameters_t ecu;
  GetPart("ew_cap_up",ecu);
  ecu.r_max = ecu.r_min;
  ecu.r_min = 0.*cm;  
  ecu.z     = -z0 - ecu.z;  
  ecu.ry    = 180.*deg;  

  G4Sphere *endcap_up = new G4Sphere(ecu.name,
	                             ecu.r_min     ,ecu.r_max,
	                             ecu.startPhi  ,ecu.dPhi,
	                             ecu.startTheta,ecu.dTheta);

  G4ThreeVector P_ecu      = G4ThreeVector(ecu.x,ecu.y,ecu.z);  
  G4RotationMatrix *rm_ecu = new G4RotationMatrix();
  rm_ecu->rotateX(ecu.rx);   
  rm_ecu->rotateY(ecu.ry);   
  rm_ecu->rotateZ(ecu.rz);  

  //---- end window (downstream, new)

  // build the main shaft 
  partParameters_t mshd = mshu;
  mshd.name = "mshd_dn"; 
  mshd.z *= -1.; 

  G4Tubs *mainShaft_dn = new G4Tubs(mshd.name,
	                            mshd.r_min    ,mshd.r_max,
	                            mshd.length/2.,
	                            mshd.startPhi ,mshd.dPhi);

  G4ThreeVector P_mshd      = G4ThreeVector(mshd.x,mshd.y,mshd.z);  
  G4RotationMatrix *rm_mshd = new G4RotationMatrix();
  rm_mshd->rotateX(mshd.rx);   
  rm_mshd->rotateY(mshd.ry);   
  rm_mshd->rotateZ(mshd.rz);   

  // lip  
  partParameters_t lipd = lipu;
  lipd.name = "lipd_dn"; 
  lipd.z *= -1.;

  G4Tubs *lip_dn = new G4Tubs(lipd.name,
	                      lipd.r_min    ,lipd.r_max,
	                      lipd.length/2.,
	                      lipd.startPhi ,lipd.dPhi);

  G4ThreeVector P_lipd      = G4ThreeVector(lipd.x,lipd.y,lipd.z);  
  G4RotationMatrix *rm_lipd = new G4RotationMatrix();
  rm_lipd->rotateX(lipd.rx);   
  rm_lipd->rotateY(lipd.ry);   
  rm_lipd->rotateZ(lipd.rz);   

  // rounded lip 
  partParameters_t rlipd = rlipu;
  rlipd.name = "rlipd_dn"; 
  rlipd.z *= -1.; 
  rlipd.ry = 0.*deg; 

  G4Sphere *roundLip_dn = new G4Sphere(rlipd.name,
	                               rlipd.r_min     ,rlipd.r_max,
	                               rlipd.startPhi  ,rlipd.dPhi,
	                               rlipd.startTheta,rlipd.dTheta);

  G4ThreeVector P_rlipd      = G4ThreeVector(rlipd.x,rlipd.y,rlipd.z);  
  G4RotationMatrix *rm_rlipd = new G4RotationMatrix();
  rm_rlipd->rotateX(rlipd.rx);   
  rm_rlipd->rotateY(rlipd.ry);   
  rm_rlipd->rotateZ(rlipd.rz);   

  // endcap 
  partParameters_t ecd = ecu;
  ecd.name = "ecd_dn"; 
  ecd.z *= -1.; 
  ecd.ry = 0.*deg; 

  G4Sphere *endcap_dn = new G4Sphere(ecd.name,
	                             ecd.r_min     ,ecd.r_max,
	                             ecd.startPhi  ,ecd.dPhi,
	                             ecd.startTheta,ecd.dTheta);

  G4ThreeVector P_ecd      = G4ThreeVector(ecd.x,ecd.y,ecd.z);  
  G4RotationMatrix *rm_ecd = new G4RotationMatrix();
  rm_ecd->rotateX(ecd.rx);   
  rm_ecd->rotateY(ecd.ry);   
  rm_ecd->rotateZ(ecd.rz);  

  // end window (upstream)
  // combine info from lip, rlip, end cap to one object here
  // try to make this a simple union of a cylinder and a hemisphere.  
  // TODO: may need something more involved later... 

  // // build the "base" -- combines end window main shaft + lip   
  // partParameters_t mshu;
  // GetPart("ew_mainShaft_up",mshu);

  // partParameters_t lipu;
  // GetPart("ew_lip_up",lipu);

  // partParameters_t baseu; 
  // baseu.name       = "He3_base_u"; 
  // baseu.r_max      = tc.r_max; // mshu.r_min; 
  // baseu.r_min      = 0.*cm;
  // baseu.length     = mshu.length + lipu.length;  
  // baseu.x          = 0.*cm; 
  // baseu.y          = 0.*cm; 
  // baseu.z          = (-1.)*( 0.5*tc.length + 0.5*baseu.length ); 
  // baseu.startTheta = 0.*deg; 
  // baseu.dTheta     = 0.*deg; 
  // baseu.startPhi   = 0.*deg; 
  // baseu.dPhi       = 360.*deg; 

  // G4Tubs *base_u = new G4Tubs(baseu.name,
  //                             baseu.r_min    ,baseu.r_max,
  //                             baseu.length/2.,
  //                             baseu.startPhi ,baseu.dPhi);

  // G4ThreeVector P_baseu      = G4ThreeVector(baseu.x,baseu.y,baseu.z);  
  // G4RotationMatrix *rm_baseu = new G4RotationMatrix();
  // rm_baseu->rotateX(baseu.rx);   
  // rm_baseu->rotateY(baseu.ry);   
  // rm_baseu->rotateZ(baseu.rz);   
  // 
  // // hemisphere
  // partParameters_t hspu; 
  // hspu.name       = "He3_hspu"; 
  // hspu.r_max      = tc.r_max;  
  // hspu.r_min      = 0.*cm;
  // hspu.length     = 0.*cm;  
  // hspu.x          = 0.*cm; 
  // hspu.y          = 0.*cm; 
  // hspu.z          = baseu.z - 0.5*hspu.r_max; 
  // hspu.rx         = 0.*deg; 
  // hspu.ry         = 180.*deg; 
  // hspu.rz         = 0.*deg; 
  // hspu.startTheta = 0.*deg; 
  // hspu.dTheta     = 90.*deg; 
  // hspu.startPhi   = 0.*deg; 
  // hspu.dPhi       = 360.*deg; 

  // G4Sphere *hSphere_u = new G4Sphere(hspu.name,
  //                                    hspu.r_min      ,hspu.r_max,
  //                                    hspu.startPhi   ,hspu.dPhi,
  //                                    hspu.startTheta ,hspu.dTheta);

  // G4ThreeVector P_hspu      = G4ThreeVector(hspu.x,hspu.y,hspu.z);  
  // G4RotationMatrix *rm_hspu = new G4RotationMatrix();
  // rm_hspu->rotateX(hspu.rx);   
  // rm_hspu->rotateY(hspu.ry);   
  // rm_hspu->rotateZ(hspu.rz);   

  // downstream 
  // partParameters_t based = baseu; 
  // based.name = "He3_base_d"; 
  // based.z *= (-1.); 

  // G4Tubs *base_d = new G4Tubs(based.name,
  //                             based.r_min    ,based.r_max,
  //                             based.length/2.,
  //                             based.startPhi ,based.dPhi);

  // G4ThreeVector P_based      = G4ThreeVector(based.x,based.y,based.z);  
  // G4RotationMatrix *rm_based = new G4RotationMatrix();
  // rm_based->rotateX(based.rx);   
  // rm_based->rotateY(based.ry);   
  // rm_based->rotateZ(based.rz);   

  // partParameters_t hspd = hspu;
  // hspd.name  = "He3_hspd";  
  // hspd.z    *= (-1.); 
  // hspd.ry    = 0.*deg;

  // G4Sphere *hSphere_d = new G4Sphere(hspd.name,
  //                                    hspd.r_min      ,hspd.r_max,
  //                                    hspd.startPhi   ,hspd.dPhi,
  //                                    hspd.startTheta ,hspd.dTheta);

  // G4ThreeVector P_hspd      = G4ThreeVector(hspd.x,hspd.y,hspd.z);  
  // G4RotationMatrix *rm_hspd = new G4RotationMatrix();
  // rm_hspd->rotateX(hspd.rx);   
  // rm_hspd->rotateY(hspd.ry);   
  // rm_hspd->rotateZ(hspd.rz);   

  // Union solid 
  // use same rotation and positional vectors as glass shell!  
  G4UnionSolid *he3Tube;
  // // target chamber + base upstream  
  // he3Tube = new G4UnionSolid("tc_bu"     ,tcShape,base_u   ,rm_baseu,P_baseu);
  // // assembly + hsphere upstream 
  // he3Tube = new G4UnionSolid("tc_bhu"    ,he3Tube,hSphere_u,rm_hspu,P_hspu);  
  // // assembly + base downstream 
  // he3Tube = new G4UnionSolid("tcbhu_hspd",he3Tube,base_d   ,rm_based,P_based);  
  // // assembly + hsphere downstream 
  // he3Tube = new G4UnionSolid("heTube"    ,he3Tube,hSphere_d,rm_hspd,P_hspd);  

  // more complicated geometry 
  he3Tube = new G4UnionSolid("tc_um"        ,tcShape,mainShaft_up,rm_mshu ,P_mshu );
  he3Tube = new G4UnionSolid("tc_uml"       ,he3Tube,lip_up      ,rm_lipu ,P_lipu );
  he3Tube = new G4UnionSolid("tc_umlr"      ,he3Tube,roundLip_up ,rm_rlipu,P_rlipu);
  he3Tube = new G4UnionSolid("tc_umlre"     ,he3Tube,endcap_up   ,rm_ecu  ,P_ecu  );
  // upstream end window  
  he3Tube = new G4UnionSolid("tc_umlre_dm"  ,he3Tube,mainShaft_dn,rm_mshd ,P_mshd );
  he3Tube = new G4UnionSolid("tc_umlre_dml" ,he3Tube,lip_dn      ,rm_lipd ,P_lipd );
  he3Tube = new G4UnionSolid("tc_umlre_dmlr",he3Tube,roundLip_dn ,rm_rlipd,P_rlipd);
  he3Tube = new G4UnionSolid("he3Tube"      ,he3Tube,endcap_dn   ,rm_ecd  ,P_ecd  );

  // set the color of He3 
  G4VisAttributes *visHe3 = new G4VisAttributes(); 
  visHe3->SetColour( G4Colour::Yellow() );
  // visHe3->SetForceWireframe(true);  

  // logical volume of He3
  G4LogicalVolume *logicHe3 = new G4LogicalVolume(he3Tube,GetMaterial("He3"),"logicHe3");  
  logicHe3->SetVisAttributes(visHe3);  
 
  // placement of He3 is *inside target chamber*  
  G4ThreeVector posHe3 = G4ThreeVector(0.*cm,0.*cm,0.*cm); 
  new G4PVPlacement(0,                 // rotation
                    posHe3,            // position 
                    logicHe3,          // logical volume 
                    "physHe3",         // name 
                    logicMother,       // logical mother volume is the target chamber 
                    false,             // no boolean operations 
                    0,                 // copy number 
                    fCheckOverlaps);    // check overlaps

}
//______________________________________________________________________________
void He3TargetDetectorConstruction::BuildEndWindow(const std::string type,G4LogicalVolume *logicMother){
   // build the metal end window based on type = upstream or downstream 

   std::string suffix; 
   if( type.compare("upstream")==0 ){
      suffix = "up"; 
   }else if( type.compare("downstream")==0 ){
      suffix = "dn"; 
   }else{
      std::cout << "[He3TargetDetectorConstruction::BuildEndWindow]: ERROR! Invalid type = " << type << std::endl;
      exit(1);
   }  

   std::string ms_str = "ew_mainShaft_" + suffix; 
   std::string l_str  = "ew_lip_"       + suffix; 
   std::string rl_str = "ew_rlip_"      + suffix; 
   std::string c_str  = "ew_cap_"       + suffix; 

   // build the main shaft 
   partParameters_t msh;
   GetPart(ms_str.c_str(),msh); 

   G4Tubs *mainShaft = new G4Tubs(msh.name,
	                          msh.r_min    ,msh.r_max,
	                          msh.length/2.,
	                          msh.startPhi ,msh.dPhi);

   G4ThreeVector P_msh      = G4ThreeVector(msh.x,msh.y,msh.z);  
   G4RotationMatrix *rm_msh = new G4RotationMatrix();
   rm_msh->rotateX(msh.rx);   
   rm_msh->rotateY(msh.ry);   
   rm_msh->rotateZ(msh.rz);   

   // lip  
   partParameters_t lip;
   GetPart(l_str.c_str(),lip); 

   G4Tubs *lipTube = new G4Tubs(lip.name,
	                        lip.r_min    ,lip.r_max,
	                        lip.length/2.,
	                        lip.startPhi ,lip.dPhi);

   G4ThreeVector P_lip      = G4ThreeVector(lip.x,lip.y,lip.z);  
   G4RotationMatrix *rm_lip = new G4RotationMatrix();
   rm_lip->rotateX(lip.rx);   
   rm_lip->rotateY(lip.ry);   
   rm_lip->rotateZ(lip.rz);   

   // rounded lip 
   partParameters_t rlip;
   GetPart(rl_str.c_str(),rlip); 

   G4Sphere *roundLip = new G4Sphere(rlip.name,
	                             rlip.r_min     ,rlip.r_max,
	                             rlip.startPhi  ,rlip.dPhi,
	                             rlip.startTheta,rlip.dTheta);

   G4ThreeVector P_rlip = G4ThreeVector(rlip.x,rlip.y,rlip.z);  
   G4RotationMatrix *rm_rlip = new G4RotationMatrix();
   rm_rlip->rotateX(rlip.rx);   
   rm_rlip->rotateY(rlip.ry);   
   rm_rlip->rotateZ(rlip.rz);   

   // endcap 
   partParameters_t ec;
   GetPart(c_str.c_str(),ec); 

   G4Sphere *endcap = new G4Sphere(ec.name,
	                           ec.r_min     ,ec.r_max,
	                           ec.startPhi  ,ec.dPhi,
	                           ec.startTheta,ec.dTheta);

   G4ThreeVector P_ec      = G4ThreeVector(ec.x,ec.y,ec.z);  
   G4RotationMatrix *rm_ec = new G4RotationMatrix();
   rm_ec->rotateX(ec.rx);   
   rm_ec->rotateY(ec.ry);   
   rm_ec->rotateZ(ec.rz);  

   // create the union solid 
   G4UnionSolid *endWindow; 
   // main shaft + lip 
   G4String label1 = "ew_ms_l_"    + suffix; 
   G4String label2 = "ew_ms_l_rl_" + suffix; 
   G4String label3 = "endWindow_"  + suffix; 
   endWindow = new G4UnionSolid(label1,mainShaft,lipTube ,rm_lip ,P_lip); 
   // add rounded lip 
   endWindow = new G4UnionSolid(label2,endWindow,roundLip,rm_rlip,P_rlip);
   // endcap  
   endWindow = new G4UnionSolid(label3,endWindow,endcap  ,rm_ec  ,P_ec); 

   // create the logical volume
   G4String logicLabel = "logicEndWindow_" + suffix;  
   G4LogicalVolume *logicEndWindow = new G4LogicalVolume(endWindow,GetMaterial("Copper"),logicLabel);

   // visualization
   G4VisAttributes *vis = new G4VisAttributes();
   vis->SetColour( G4Colour::Red() ); 
   vis->SetForceWireframe(true);

   logicEndWindow->SetVisAttributes(vis); 

   // place it at the end of the target chamber 
   partParameters_t tc; 
   GetPart("targetChamber",tc); 
   G4double x_ew = 0.*cm; 
   G4double y_ew = 0.*cm;   
   G4double z_ew = (tc.length/2. + 0.5*msh.length);  // this will make the large lip flush with the glass
   G4double ra[3] = {0.*deg,0.*deg,0.*deg}; 

   // adjustments for upstream part 
   if( type.compare("upstream")==0 ){
      z_ew *= -1.;
      ra[1] = 180.*deg;  
   }

   G4ThreeVector P_ew = G4ThreeVector(x_ew,y_ew,z_ew); // location  
   G4RotationMatrix *rm_ew = new G4RotationMatrix(); 
   rm_ew->rotateX(ra[0]); 
   rm_ew->rotateY(ra[1]); 
   rm_ew->rotateZ(ra[2]);
   
   G4String physLabel = "physCuEndWindow_" + suffix;  

   new G4PVPlacement(rm_ew,
	             P_ew, 
	             logicEndWindow,      // logical volume 
	             physLabel,           // name 
	             logicMother,         // logical mother volume is the target chamber 
	             false,               // no boolean operations 
	             0,                   // copy number 
	             fCheckOverlaps);     // check overlaps
 
}
//______________________________________________________________________________
void He3TargetDetectorConstruction::BuildShield(int config,G4LogicalVolume *logicMother){
   // shield box for the target magnetic field 
   // - Material: Carbon steel 1008
   // - config: kSBS_GEN_146, kSBS_GEN_368, kSBS_GEN_677, kSBS_GEN_1018
   // - The shield is actually two layers
   //   - each layer is 0.25" thick
   //   - outer surfaces of layers separated by 1.29" 
   //     => 1.29 - 0.25 = 1.04" center-to-center distance 
   // - Drawing number: A09016-03-05-0000, A09016-03-05-0800

   // constants from drawings  
   G4double wall = 0.25*2.54*cm;
   G4double sp   = 0.79*2.54*cm;  // inner-spacing: |<-s->| 

   // get outer box dimensions  
   partParameters_t sh; 
   GetPart("shield",sh); 

   //---- window cuts
   // FIXME: these will change!     
   // downstream, along beam.  sizes are estimates! 
   G4double xw = 9.48*2.54*cm; 
   G4double yw = 9.48*2.54*cm; 
   G4double zw = 9.48*2.54*cm; 
   G4Box *windowCut_dn = new G4Box("windowCut_dn",xw/2.,yw/2.,zw/2.); 
   G4ThreeVector Pw_dn = G4ThreeVector(sh.x_len/2.,0.,sh.z_len/3.);  // position of cut 
   
   G4double door = 40.5*2.54*cm; // for OWU2A door panel (from JT model)  

   // downstream, beam left 
   // x = 4.62 + ??, y = 28.44", z = 4.62 + ?? 
   // drawings: A09016-03-05-0851  
   // cut as of 6/27/20 (new design from Bert); comments are old design 

   // old 
   // G4double xw_bl = 15.*2.54*cm; 
   // G4double yw_bl = 28.44*2.54*cm; 
   // G4double zw_bl = 15.*2.54*cm; 
   // G4ThreeVector Pw_bl_dn = G4ThreeVector(sh.x_len/2.,0.,0.);  // position of cut 

   // new cut as of 6/27/20 (new design from Bert)
   G4double xw_bl = 10.*cm; // arbitrary size
   G4double yw_bl = 0.8*sh.y_len; 
   G4double zw_bl = door; 

   // coordinates 
   G4double xs = sh.x_len/2.; // distance to midpoint of door, aligns door close to edge 
   G4double ys = (-1.)*(sh.y_len/2.-yw_bl/2.);          // this should center the cut properly  
   G4double zs = -zw_bl/2.; // right on top of the right side wall (before rotation)  
   G4ThreeVector Pw_bl_dn = G4ThreeVector(xs,ys,zs);  // position of cut 
  
   G4Box *windowCut_beamLeft_dn = new G4Box("windowCut_beamLeft_dn",xw_bl/2.,yw_bl/2.,zw_bl/2.); 

   // downstream, beam right
   // FIXME: might be affected by SBS configuration
   config +=0;  // to get rid of compiler warnings for now 
   // for now consider a single large cut 
   // large opening; distance along x and z = panel 1 + door w/handles + panel 3
   // panel 1
   // - x = 8.62", y = 21.41", z = 8.62" 
   // - drawing: A09016-03-05-0811
   // panel 2 
   // - x = 6.62", y = 21.41", z = 6.62" 
   // door with handles  
   // - x = 16.70", y = 21.41", z = 16.70" 
   // - drawing: A09016-03-05-0841
   // panel 3
   // - x = 5.12", y = 21.41", z = 5.12" 
   // - drawing: A09016-03-05-0831
   // G4double p1 = 8.62*2.54*cm; 
   // G4double p2 = 6.62*2.54*cm; 
   // G4double p3 = 5.12*2.54*cm; 
   // G4double dh = 16.70*2.54*cm; 
   // G4double T  = p1+p2+p3+dh;  // total cutaway size 

   // old 
   // G4double xw_br = T; 
   // G4double yw_br = 21.41*2.54*cm; 
   // G4double zw_br = 5.*cm; // arbitray size
   // // coordinates 
   // G4double xs = sh.x_len/2. - xw_br/2. - 2.*cm; // distance to midpoint of door, aligns door close to edge 
   // G4double ys = (-1.)*(sh.y_len/2.-yw_br/2.);  // this should center the cut properly  
   // G4double zs = sh.z_len/2.; // right on top of the right side wall (before rotation)  
 
   // cut as of 6/27/20 (new design from Bert)  
   G4double xw_br = door; 
   G4double yw_br = 0.8*sh.y_len; 
   G4double zw_br = 10.*cm; 
   // coordinates 
   xs = sh.x_len/2. - door/2. - 5.*cm; // distance to midpoint of door, aligns door close to edge 
   ys = (-1.)*(sh.y_len/2.-yw_br/2.);  // this should center the cut properly  
   zs = sh.z_len/2.; // right on top of the right side wall (before rotation)  
 
   G4Box *windowCut_beamRight_dn = new G4Box("windowCut_beamRight_dn",xw_br/2.,yw_br/2.,zw_br/2.); 
   G4ThreeVector Pw_br_dn        = G4ThreeVector(xs,ys,zs);  // position of cut 
  
   // upstream, along beam  
   G4Box *windowCut_up = new G4Box("windowCut_up",xw/2.,yw/2.,zw/2.); 
   G4ThreeVector Pw_up = G4ThreeVector(-sh.x_len/2.,0.,-sh.z_len/3.);  // position of cut 
   
   //---- shield  

   // outer box 
   G4Box *outer = new G4Box("outer"   ,sh.x_len/2.,sh.y_len/2.,sh.z_len/2.);
   // cut away inner material  
   G4double xc = sh.x_len - wall*2.; 
   G4double yc = sh.y_len - wall*2.; 
   G4double zc = sh.z_len - wall*2.; 
   G4Box *outerCut = new G4Box("outerCut",xc/2.,yc/2.,zc/2.);

   // subtract the parts 
   G4SubtractionSolid *outerShield = new G4SubtractionSolid("outerShield_1",outer,outerCut,0,G4ThreeVector(0,0,0));
   outerShield = new G4SubtractionSolid("outerShield_2",outerShield,windowCut_dn,0,Pw_dn); 
   outerShield = new G4SubtractionSolid("outerShield_3",outerShield,windowCut_beamLeft_dn,0,Pw_bl_dn); 
   outerShield = new G4SubtractionSolid("outerShield_4",outerShield,windowCut_beamRight_dn,0,Pw_br_dn); 
   outerShield = new G4SubtractionSolid("outerShield"  ,outerShield,windowCut_up,0,Pw_up); 

   // logical volume
   G4VisAttributes *visOut = new G4VisAttributes();
   // visOut->SetForceWireframe(true); 
   visOut->SetColour( G4Colour::Magenta() ); 
   
   G4LogicalVolume *logicOuter     = new G4LogicalVolume(outerShield,GetMaterial("Carbon_Steel_1008"),"logicOuter");
   logicOuter->SetVisAttributes(visOut);  

   // inner box 
   partParameters_t sh_inner = sh; 
   sh_inner.x_len = sh.x_len - 2.*wall - 2.*sp; 
   sh_inner.y_len = sh.y_len - 2.*wall - 2.*sp; 
   sh_inner.z_len = sh.z_len - 2.*wall - 2.*sp; 

   G4Box *inner = new G4Box("inner",sh_inner.x_len/2.,sh_inner.y_len/2.,sh_inner.z_len/2.);
   // cut away inner material  
   xc = sh_inner.x_len - wall*2.; 
   yc = sh_inner.y_len - wall*2.; 
   zc = sh_inner.z_len - wall*2.; 
   G4Box *innerCut = new G4Box("innerCut",xc/2.,yc/2.,zc/2.);

   // subtract the parts 
   G4SubtractionSolid *innerShield = new G4SubtractionSolid("innerShield_1",inner,innerCut,0,G4ThreeVector(0,0,0)); 
   innerShield = new G4SubtractionSolid("innerShield_2",innerShield,windowCut_dn,0,Pw_dn); 
   innerShield = new G4SubtractionSolid("innerShield_3",innerShield,windowCut_beamLeft_dn,0,Pw_bl_dn); 
   innerShield = new G4SubtractionSolid("innerShield_4",innerShield,windowCut_beamRight_dn,0,Pw_br_dn); 
   innerShield = new G4SubtractionSolid("innerShield",innerShield,windowCut_up,0,Pw_up); 

   // logical volume
   G4VisAttributes *visIn = new G4VisAttributes();
   // visIn->SetForceWireframe(true); 
   visIn->SetColour( G4Colour::Magenta() ); 

   G4LogicalVolume *logicInner = new G4LogicalVolume(innerShield,GetMaterial("Carbon_Steel_1008"),"logicInner");
   logicInner->SetVisAttributes(visIn);  

   // place the parts 
   G4double RY = 55.0*deg;  // FIXME: This angle is still an estimate!  

   G4RotationMatrix *rm = new G4RotationMatrix();
   rm->rotateX(0.*deg); rm->rotateY(RY); rm->rotateZ(0.*deg);

   G4ThreeVector P = G4ThreeVector(0.*cm,0.*cm,0.*cm); 

   bool isBoolean = true; 

   // outer
   new G4PVPlacement(rm,                // rotation relative to mother       
                     P,                 // position relative to mother         
                     logicOuter,        // logical volume        
                     "physShieldOuter", // physical volume name           
                     logicMother,       // logical mother     
                     isBoolean,         // is boolean device? (true or false)    
                     0,                 // copy number    
                     fCheckOverlaps);   // check overlaps  

   // inner
   new G4PVPlacement(rm,                // rotation relative to mother       
                     P,                 // position relative to mother         
                     logicInner,        // logical volume        
                     "physShieldInner", // physical volume name           
                     logicMother,       // logical mother     
                     isBoolean,         // is boolean device? (true or false)    
                     0,                 // copy number    
                     fCheckOverlaps);   // check overlaps

}
//______________________________________________________________________________
void He3TargetDetectorConstruction::BuildLadderPlate(G4LogicalVolume *logicMother){
   // ladder plate  
   // - This comes close to the beam path
   // - Drawing number: A09016-03-04-0601
   G4double x0 = -0.438*2.54*cm; // beam right from JT model  
   G4double y0 = -5.33*cm;       // low according to JT model
   G4double z0 = 1.5*2.54*cm;   // downstream according to JT model   

   //---- vertical posts along y axis 
   // upstream 
   partParameters_t lvu; 
   GetPart("ladder_vert_up",lvu); 

   G4Box *ladder_vert_up    = new G4Box("lvu",lvu.x_len/2.,lvu.y_len/2.,lvu.z_len/2.);
   G4ThreeVector P_lvu      = G4ThreeVector(lvu.x,lvu.y,lvu.z);  
   G4RotationMatrix *rm_lvu = new G4RotationMatrix(); 
   rm_lvu->rotateX(lvu.rx); rm_lvu->rotateY(lvu.ry); rm_lvu->rotateZ(lvu.rz);   

   // downstream 
   partParameters_t lvd; 
   GetPart("ladder_vert_dn",lvd); 

   G4Box *ladder_vert_dn    = new G4Box("lvd",lvd.x_len/2.,lvd.y_len/2.,lvd.z_len/2.);
   G4ThreeVector P_lvd      = G4ThreeVector(lvd.x,lvd.y,lvd.z);  
   G4RotationMatrix *rm_lvd = new G4RotationMatrix(); 
   rm_lvd->rotateX(lvd.rx); rm_lvd->rotateY(lvd.ry); rm_lvd->rotateZ(lvd.rz);   

   //---- horizontal posts along z axis 
   partParameters_t la;  
   GetPart("ladder_above",la); 

   G4Box *ladder_above     = new G4Box("la",la.x_len/2.,la.y_len/2.,la.z_len/2.);
   G4ThreeVector P_la      = G4ThreeVector(la.x,la.y,la.z);  
   G4RotationMatrix *rm_la = new G4RotationMatrix(); 
   rm_la->rotateX(la.rx); rm_la->rotateY(la.ry); rm_la->rotateZ(la.rz);   

   //---- horizontal part below, along z axis 
   // FIXME: This is an approximated length in z!
   partParameters_t lb; 
   GetPart("ladder_below",lb); 

   G4Box *ladder_below     = new G4Box("lb",lb.x_len/2.,lb.y_len/2.,lb.z_len/2.);
   G4ThreeVector P_lb      = G4ThreeVector(lb.x,lb.y,lb.z);  
   G4RotationMatrix *rm_lb = new G4RotationMatrix(); 
   rm_lb->rotateX(lb.rx); rm_lvd->rotateY(lb.ry); rm_lb->rotateZ(lb.rz);   

   //---- angular part above, along z axis, upstream 
   partParameters_t aau; 
   GetPart("ladder_ang_above_up",aau); 

   G4Box *ladder_aau        = new G4Box("aau",aau.x_len/2.,aau.y_len/2.,aau.z_len/2.);
   G4ThreeVector P_aau      = G4ThreeVector(aau.x,aau.y,aau.z);  
   G4RotationMatrix *rm_aau = new G4RotationMatrix(); 
   rm_aau->rotateX(aau.rx); rm_aau->rotateY(aau.ry); rm_aau->rotateZ(aau.rz);   

   //---- angular part above, along z axis, downstream 
   partParameters_t aad; 
   GetPart("ladder_ang_above_dn",aad); 

   G4Box *ladder_aad        = new G4Box("aad",aad.x_len/2.,aad.y_len/2.,aad.z_len/2.);
   G4ThreeVector P_aad      = G4ThreeVector(aad.x,aad.y,aad.z);  
   G4RotationMatrix *rm_aad = new G4RotationMatrix(); 
   rm_aad->rotateX(aad.rx); rm_aad->rotateY(aad.ry); rm_aad->rotateZ(aad.rz);   

   //---- angular part below, along z axis, upstream 
   partParameters_t abu; 
   GetPart("ladder_ang_below_up",abu); 

   G4Box *ladder_abu        = new G4Box("abu",abu.x_len/2.,abu.y_len/2.,abu.z_len/2.);
   G4ThreeVector P_abu      = G4ThreeVector(abu.x,abu.y,abu.z);  
   G4RotationMatrix *rm_abu = new G4RotationMatrix(); 
   rm_abu->rotateX(abu.rx); rm_abu->rotateY(abu.ry); rm_abu->rotateZ(abu.rz);   

   //---- angular part below, along z axis, downstream 
   partParameters_t abd; 
   GetPart("ladder_ang_below_dn",abd); 

   G4Box *ladder_abd        = new G4Box("abd",abd.x_len/2.,abd.y_len/2.,abd.z_len/2.);
   G4ThreeVector P_abd      = G4ThreeVector(abd.x,abd.y,abd.z);  
   G4RotationMatrix *rm_abd = new G4RotationMatrix(); 
   rm_abd->rotateX(abd.rx); rm_abd->rotateY(abd.ry); rm_abd->rotateZ(abd.rz);   

   // unions -- vectors are relative to the upstream post  
   G4UnionSolid *ladder; 
   // [downstream] pane inner and outer  
   G4ThreeVector P     = G4ThreeVector(0.,0.,lvd.z*2.);     
   ladder = new G4UnionSolid("lud",ladder_vert_up,ladder_vert_dn,0,P); 
   // [above] horizontal 
   G4ThreeVector P_a   = G4ThreeVector(0.,la.y,la.z);     
   ladder = new G4UnionSolid("lud_a" ,ladder,ladder_above,0,P_a); 
   // [below] horizontal 
   G4ThreeVector P_b   = G4ThreeVector(0.,lb.y,lb.z);     
   ladder = new G4UnionSolid("lud_ab",ladder,ladder_below,0,P_b); 
   // [above] angled above, upstream  
   ladder = new G4UnionSolid("lud_ab_aau",ladder,ladder_aau,rm_aau,P_aau); 
   // [above] angled below, downstream  
   ladder = new G4UnionSolid("lud_ab_aaud",ladder,ladder_aad,rm_aad,P_aad); 
   // [below] angled below, upstream  
   ladder = new G4UnionSolid("lud_ab_aaud_abu",ladder,ladder_abu,rm_abu,P_abu); 
   // [below] angled, downstream  
   ladder = new G4UnionSolid("ladder",ladder,ladder_abd,rm_abd,P_abd); 

   G4VisAttributes *vis = new G4VisAttributes(); 
   vis->SetColour( G4Colour::Red() ); 
   // vis->SetForceWireframe(true);

   G4LogicalVolume *logicLadder = new G4LogicalVolume(ladder,GetMaterial("Aluminum"),"logicLadder"); 
   logicLadder->SetVisAttributes(vis); 

   G4double lx = x0; 
   G4double ly = y0;  
   G4double lz = lvu.z + z0; 
    
   G4ThreeVector P_l      = G4ThreeVector(lx,ly,lz);   
   G4RotationMatrix *rm_l = new G4RotationMatrix();
   rm_l->rotateX(0.*deg); rm_l->rotateY(0.*deg); rm_l->rotateZ(0.*deg); 
 
   new G4PVPlacement(rm_l,                // rotation [relative to mother]    
                     P_l,                 // position [relative to mother] 
                     logicLadder,         // logical volume    
                     "physLadder",        // name                          
                     logicMother,         // logical mother volume            
                     false,               // no boolean operations          
                     0,                   // copy number                   
                     fCheckOverlaps);     // check overlaps               
   
}
//______________________________________________________________________________
void He3TargetDetectorConstruction::BuildPickupCoils(G4LogicalVolume *logicMother){
   // Pickup coils used to determine polarization of 3He (NMR, EPR) 
   // - Dimensions based on JT model
   // - global y offset: pickup coils sit 1.1" below the target cell (measured from top of coil mount)
   G4double y0 = -1.1*2.54*cm - 1.5*cm;  // 1.5 cm accounts for center of coil mount  

   //---- coil mount [upstream, beam left] 
   // note: the subtraction dimensions are estimates at best  
   partParameters_t cbul; 
   GetPart("pu_coil_base",cbul);

   G4Box *coilB = new G4Box("coilB",cbul.x_len/2.,cbul.y_len/2.,cbul.z_len/2.); 


   // U-cut: create a *subtraction* with a component for a U shape  
   G4double xlen = 2.*cbul.x_len;  // far exceed the thickness to make sure it's a cutaway 
   G4double ylen = 2.*1.3*cm; // cmul.y_len/2.;
   G4double zlen = 3.2*cm; // cmul.z_len/3.;

   G4Box *Ucut = new G4Box("Ucut",xlen/2.,ylen/2.,zlen/2.);

   // do the subtraction 
   G4ThreeVector Psub = G4ThreeVector(0,1.1*cm+ylen/2.,0.);  // centers the second volume relative to the first 
   G4SubtractionSolid *coilBase = new G4SubtractionSolid("coilMount",coilB,Ucut,0,Psub); 
 
   // now the coil itself 
   partParameters_t pu;
   GetPart("pu_coil",pu); 

   G4Box *coilOuter = new G4Box("coilOuter",pu.x_len/2.,pu.y_len/2.,pu.z_len/2.);

   std::cout << "PU COIL THICKNESS = " << pu.x_len/cm << " cm" << std::endl;

   // create a cutaway that actually defines the coil since the initial dimensions are the OUTER values 
   G4double xc = 2.*pu.x_len; // cut straight through in this dimension  
   G4double yc = pu.y_len - 2.*0.5*cm;  
   G4double zc = pu.z_len - 2.*0.5*cm; 
   
   G4Box *coilCut = new G4Box("coilCut",xc,yc/2.,zc/2.); 
   
   // subtraction; no coordinates needed since we're centered on the outer coil  
   G4SubtractionSolid *coil = new G4SubtractionSolid("coil",coilOuter,coilCut);

   // ok, now we use the *coil* to do a cut on the coil mount (this is where the coil will sit)
   partParameters_t cmul; 
   GetPart("pu_coil_mnt",cmul);
   
   G4Box *coilMnt = new G4Box("cmnt",cmul.x_len/2.,cmul.y_len/2.,cmul.z_len/2.); 

   // need to locate the coil properly
   xc = pu.x_len/2.;  
   yc = 0.*cm; 
   zc = 0.*cm; 
   G4ThreeVector Pc = G4ThreeVector(xc,yc,zc);
   G4SubtractionSolid *coilMount = new G4SubtractionSolid("coilMount",coilMnt,coil,0,Pc); 

   // now create a union of the coil base and coil mount
   G4double xm = cbul.x_len/2. + cmul.x_len/2.; 
   G4double ym = 0.9*cm; 
   G4double zm = 0.*cm;  
   G4ThreeVector Pm = G4ThreeVector(xm,ym,zm);  
   G4UnionSolid *coilBMNT = new G4UnionSolid("coilBMNT",coilBase,coilMount,0,Pm); 
  
   G4VisAttributes *vis = new G4VisAttributes();
   vis->SetColour( G4Colour::Magenta() );  
   // vis->SetForceWireframe(true);  

   G4LogicalVolume *logicCoilBMNT = new G4LogicalVolume(coilBMNT,GetMaterial("Aluminum"),"logicCoilBMNT");
   logicCoilBMNT->SetVisAttributes(vis); 

   bool isBoolean = true;  // not sure what this really does...  

   // FIXME: constraints come from JT file; is the 2.5 cm correct?
   G4double xbm = (2.5*cm)/2. + cbul.x_len/2. +  cmul.x_len + pu.x_len; 
   G4double ybm = -0.9*cm + y0; // placement relative to BEAM 
   G4double zbm = 7.2*cm; 
   
   G4ThreeVector P_ul = G4ThreeVector(xbm,ybm,zbm);
   G4RotationMatrix *rm = new G4RotationMatrix(); 
   rm->rotateX(0); rm->rotateY(0); rm->rotateZ(0);   

   G4RotationMatrix *rm180 = new G4RotationMatrix(); 
   rm180->rotateX(0); rm180->rotateY(180*deg); rm180->rotateZ(0);   

   // beam left
   new G4PVPlacement(rm180,               // rotation [relative to mother]             
                     P_ul,                // position [relative to mother]          
                     logicCoilBMNT,       // logical volume                            
                     "physPUCoilMNT_ubl", // name                                       
                     logicMother,         // logical mother volume                   
                     isBoolean,           // boolean operations (true, false) 
                     0,                   // copy number                          
                     fCheckOverlaps);     // check overlaps                         

   // beam right 
   G4ThreeVector P_ur = G4ThreeVector(-xbm,ybm,zbm);

   new G4PVPlacement(rm,                  // rotation [relative to mother]             
                     P_ur,                // position [relative to mother]          
                     logicCoilBMNT,       // logical volume                            
                     "physPUCoilMNT_ubr", // name                                       
                     logicMother,         // logical mother volume                   
                     isBoolean,           // boolean operations                   
                     1,                   // copy number                          
                     fCheckOverlaps);     // check overlaps                         

   //---- downstream 
   G4ThreeVector P_dl = G4ThreeVector(xbm,ybm,-zbm);

   // beam left
   new G4PVPlacement(rm180,               // rotation [relative to mother]             
                     P_dl,                // position [relative to mother]          
                     logicCoilBMNT,       // logical volume                            
                     "physPUCoilMNT_dbl", // name                                       
                     logicMother,         // logical mother volume                   
                     isBoolean,           // boolean operations (true, false)
                     2,                   // copy number                          
                     fCheckOverlaps);     // check overlaps                         

   // beam right 
   G4ThreeVector P_dr = G4ThreeVector(-xbm,ybm,-zbm);

   new G4PVPlacement(rm,                  // rotation [relative to mother]             
                     P_dr,                // position [relative to mother]          
                     logicCoilBMNT,       // logical volume                            
                     "physPUCoilMNT_dbr", // name                                       
                     logicMother,         // logical mother volume                   
                     isBoolean,           // boolean operations (true, false)
                     3,                   // copy number                          
                     fCheckOverlaps);     // check overlaps                         

   //---- pickup coil
   G4VisAttributes *visCoil = new G4VisAttributes(); 
   visCoil->SetColour( G4Colour(255,140,0) );  // dark orange
 
   G4LogicalVolume *logicCoil = new G4LogicalVolume(coil,GetMaterial("Copper"),"logicCoil");
   logicCoil->SetVisAttributes(visCoil);

   // place it to be flush against the mount assembly 
   G4double xbm_c = xbm - pu.x_len/2. - cmul.x_len; 
   G4double ybm_c = y0;   
   G4ThreeVector P_pul = G4ThreeVector(xbm_c,ybm_c,zbm); 

   // beam left
   new G4PVPlacement(rm180,               // rotation [relative to mother]             
                     P_pul,               // position [relative to mother]          
                     logicCoil,           // logical volume                            
                     "physPUCoil_ubl",    // name                                       
                     logicMother,         // logical mother volume                   
                     isBoolean,           // boolean operations (true, false) 
                     0,                   // copy number                          
                     fCheckOverlaps);     // check overlaps                         

   // beam right 
   G4ThreeVector P_pur = G4ThreeVector(-xbm_c,ybm_c,zbm);

   new G4PVPlacement(rm,                  // rotation [relative to mother]             
                     P_pur,               // position [relative to mother]          
                     logicCoil,           // logical volume                            
                     "physPUCoil_ubr",    // name                                       
                     logicMother,         // logical mother volume                   
                     isBoolean,           // boolean operations                   
                     1,                   // copy number                          
                     fCheckOverlaps);     // check overlaps                         

   //---- downstream 
   G4ThreeVector P_pdl = G4ThreeVector(xbm_c,ybm_c,-zbm);

   // beam left
   new G4PVPlacement(rm180,               // rotation [relative to mother]             
                     P_pdl,               // position [relative to mother]          
                     logicCoil,           // logical volume                            
                     "physPUCoil_dbl",    // name                                       
                     logicMother,         // logical mother volume                   
                     isBoolean,           // boolean operations (true, false)
                     2,                   // copy number                          
                     fCheckOverlaps);     // check overlaps                         

   // beam right 
   G4ThreeVector P_pdr = G4ThreeVector(-xbm_c,ybm_c,-zbm);

   new G4PVPlacement(rm,                  // rotation [relative to mother]             
                     P_pdr,               // position [relative to mother]          
                     logicCoil,           // logical volume                            
                     "physPUCoil_dbr",    // name                                       
                     logicMother,         // logical mother volume                   
                     isBoolean,           // boolean operations (true, false)
                     3,                   // copy number                          
                     fCheckOverlaps);     // check overlaps                         
 
}
//______________________________________________________________________________
void He3TargetDetectorConstruction::BuildHelmholtzCoils(int config,const std::string type,
                                                        G4LogicalVolume *logicMother){
   // Helmholtz coils for B fields 
   // - config: kSBS_GEN_677, kSBS_GEN_1018, kSBS_GEN_368, kSBS_GEN_146
   //           Different rotation angle based on index number  
   // - types: maj = large radius coil pair
   //          min = small radius coil pair 
   //          rfy = RF coil pair, aligned along the vertical (y) axis  
   // - distance between coils D = 0.5(rmin+rmax), roughly the major radius of the tube   
   //   - coil 1 (placed at -D/2) 
   //   - coil 2 (placed at +D/2)
   // - drawing number: A09016-03-08-0000

   char partName[14];
   char coilName_n[200],coilName_p[200];   // shape name 
   char coilName_nl[200],coilName_pl[200]; // logic name 
   char coilName_np[200],coilName_pp[200]; // phys name

   // shell names  
   char shellName_n[200],shellName_p[200];    
   char shellName_nl[200],shellName_pl[200];  
   char shellName_np[200],shellName_pp[200];  

   sprintf(partName   ,"helmholtz_%s",type.c_str());  
   sprintf(coilName_n ,"%s_n"        ,partName);  
   sprintf(coilName_p ,"%s_p"        ,partName);  
   sprintf(coilName_nl,"logic_%s_n"  ,partName);  
   sprintf(coilName_pl,"logic_%s_p"  ,partName);  
   sprintf(coilName_np,"phys_%s_n"   ,partName);  
   sprintf(coilName_pp,"phys_%s_p"   ,partName);  

   sprintf(shellName_n ,"shell_%s_n"        ,partName);  
   sprintf(shellName_p ,"shell_%s_p"        ,partName);  
   sprintf(shellName_nl,"shell_logic_%s_n"  ,partName);  
   sprintf(shellName_pl,"shell_logic_%s_p"  ,partName);  
   sprintf(shellName_np,"shell_phys_%s_n"   ,partName);  
   sprintf(shellName_pp,"shell_phys_%s_p"   ,partName);  

   // load parameter data 
   partParameters_t cn,cp;
   GetPart(partName,cn); 
   cp = cn;
   cn.name = coilName_n;  
   cp.name = coilName_p;  

   // coil parameters   
   G4double D      = 0.5*(cn.r_min + cn.r_max);          // helmholtz separation D = R = 0.5(rmin + rmax) 
   // G4double D      = cn.r_tor;          // helmholtz separation D = R = 0.5(rmin + rmax) 
   G4double shWall = 0;  

   if( type.compare("maj")==0 ) shWall = 5.0*mm;         // FIXME: Estimates for now! 
   if( type.compare("min")==0 ) shWall = 5.0*mm;         // FIXME: Estimates for now! 
   if( type.compare("rfy")==0 ) shWall = 0.030*25.4*mm; 

   // cylindrical geometry 
   // create the shell first 
    partParameters_t cns;
    cns.name     = shellName_n; 
    cns.r_min    = cn.r_min - shWall;  
    cns.r_max    = cn.r_max + shWall; 
    cns.length   = cn.length + 2.*shWall; // this is so we have the wall on both sides
    cns.startPhi = 0.*deg;
    cns.dPhi     = 360.*deg;   
  
    G4Tubs *cnsTube = new G4Tubs(cns.name,
                                 cns.r_min    ,cns.r_max,
                                 cns.length/2.,
                                 cns.startPhi ,cns.dPhi);
 
    partParameters_t cps;
    cps.name     = shellName_p; 
    cps.r_min    = cp.r_min - shWall;  
    cps.r_max    = cp.r_max + shWall; 
    cps.length   = cp.length + 2.*shWall; // this is so we have the wall on both sides
    cps.startPhi = 0.*deg;
    cps.dPhi     = 360.*deg;   
 
    G4Tubs *cpsTube = new G4Tubs(cps.name,
 	                         cps.r_min    ,cps.r_max,
 	                         cps.length/2.,
 	                         cps.startPhi ,cps.dPhi);

    // now create the core we subtract from the above 
    partParameters_t cns_core = cn;
    cns_core.name = cn.name + "_core"; 
    G4Tubs *cnsTube_core = new G4Tubs(cns_core.name,
                                      cns_core.r_min    ,cns_core.r_max,
                                      cns_core.length/2.,
                                      cns_core.startPhi ,cns_core.dPhi);

    partParameters_t cps_core = cp;
    cps_core.name = cp.name + "_core"; 
    G4Tubs *cpsTube_core = new G4Tubs(cps_core.name,
                                      cps_core.r_min    ,cps_core.r_max,
                                      cps_core.length/2.,
                                      cps_core.startPhi ,cps_core.dPhi);

 
    // subtract the core 
    G4SubtractionSolid *coilShell_n = new G4SubtractionSolid("cns_sub",cnsTube,cnsTube_core,0,G4ThreeVector(0,0,0));
    G4SubtractionSolid *coilShell_p = new G4SubtractionSolid("cps_sub",cpsTube,cpsTube_core,0,G4ThreeVector(0,0,0));

//    // torus geometry
//    // - some adjustments here: we should make the *thickness* of the core 
//    //   equal to the thickness specified in the drawings. 
//    //   this is better than the estimate of the derived r_max of the torus centered on 
//    //   the radius of curvature  
//    cn.r_max = cn.length/2.;    
//    cp.r_max = cp.length/2.;    
// 
//    // create the shell first 
//    partParameters_t cns;
//    cns.name     = shellName_n; 
//    cns.r_min    = cn.r_max;  
//    cns.r_max    = cn.r_max + shWall;
//    cns.r_tor    = cn.r_tor;  
//    cns.length   = cn.length;
//    cns.startPhi = 0.*deg;
//    cns.dPhi     = 360.*deg;   
//  
//    G4Torus *cnsTube = new G4Torus(cns.name,
//                                   cns.r_min   ,cns.r_max,cns.r_tor,
//                                   cns.startPhi,cns.dPhi);
// 
//    partParameters_t cps;
//    cps.name     = shellName_p; 
//    cps.r_min    = cp.r_max;  
//    cps.r_max    = cp.r_max + shWall;
//    cps.r_tor    = cp.r_tor;  
//    cps.length   = cp.length;
//    cps.startPhi = 0.*deg;
//    cps.dPhi     = 360.*deg;   
// 
//    G4Torus *cpsTube = new G4Torus(cps.name,
// 	                          cps.r_min   ,cps.r_max,cps.r_tor,
// 	                          cps.startPhi,cps.dPhi);

   // create a union of the coils
   G4ThreeVector Ps = G4ThreeVector(0.*cm,0.*cm,D);  // separated by z = D  
   G4UnionSolid *coilShells; 
   // coilShells = new G4UnionSolid("coilShells",cnsTube,cpsTube,0,Ps); // no rotation
   coilShells = new G4UnionSolid("coilShells",coilShell_n,coilShell_p,0,Ps); // no rotation

   // now for the logical volume 
   G4VisAttributes *visCoilShell = new G4VisAttributes();
   visCoilShell->SetForceWireframe(); 
   if(type.compare("maj")==0) visCoilShell->SetColour( G4Colour::Red()   ); 
   if(type.compare("rfy")==0) visCoilShell->SetColour( G4Colour::Green() ); 
   if(type.compare("min")==0) visCoilShell->SetColour( G4Colour::Blue()  ); 

   G4LogicalVolume *logicCoilShell = new G4LogicalVolume(coilShells,GetMaterial("NEMAG10"),"logicCoilShell");
   logicCoilShell->SetVisAttributes(visCoilShell);  

   // place the volume according to the input parameters  
   G4double x = cn.x;  
   G4double y = cn.y;  
   G4double z = cn.z;  

   // additional rotation to match engineering drawings (number A09016-03-08-0000) 
   G4double drx=0,dry=0,drz=0;
   if( type.compare("maj")==0 || type.compare("min")==0 ){
      if(config==kSBS_GEN_677)  dry = 10.0*deg; 
      if(config==kSBS_GEN_1018) dry = 10.0*deg; 
      if(config==kSBS_GEN_146)  dry = 43.5*deg; 
      if(config==kSBS_GEN_368)  dry = 43.5*deg; 
   } 

   G4double ph_x = cn.rx + drx; 
   G4double ph_y = cn.ry + dry; 
   G4double ph_z = cn.rz + drz; 
 
   // total rotation
   G4double COS_TOT = cos(ph_y); 
   G4double SIN_TOT = sin(ph_y);

   // adjust for y rotation.  
   // FIXME: does this exactly follow rotated coordinates?
   // rotation about y 
   // x' =  xcos + zsin 
   // z' = -xsin + zcos
 
   if(type.compare("maj")==0){ 
      x += (D/2)*fabs(SIN_TOT); 
      z += (D/2)*fabs(COS_TOT);
   }else if(type.compare("rfy")==0){
      y -= D/2.;
   }else if(type.compare("min")==0){
      x += (D/2)*fabs(SIN_TOT); 
      z -= (D/2)*fabs(COS_TOT);
   } 
 
   G4ThreeVector P_c    = G4ThreeVector(x,y,z);   
   G4RotationMatrix *rm = new G4RotationMatrix();
   rm->rotateX(ph_x); 
   rm->rotateY(ph_y); 
   rm->rotateZ(ph_z);

   char coilShellName[200];
   sprintf(coilShellName,"%s_shell",partName);  

   bool isBoolean = true; 

   new G4PVPlacement(rm,               // rotation                                        
                     P_c,              // position                                             
                     logicCoilShell,   // logical volume                                      
                     coilShellName,    // name                                                  
                     logicMother,      // logical mother volume is the target chamber          
                     isBoolean,        // boolean operations                        
                     0,                // copy number                                   
                     fCheckOverlaps);  // check overlaps                              

   // aluminum core -- goes *inside* the shell  
   G4VisAttributes *visCoil = new G4VisAttributes();
   // visCoil->SetColour( G4Colour(255,140,0) );  // dark orange 
   visCoil->SetColour( G4Colour::Grey() );  

   // cylindrical geometry 
   G4Tubs *cnTube = new G4Tubs(cn.name,
	                       cn.r_min,cn.r_max,
	                       cn.length/2.,
	                       cn.startPhi,cn.dPhi);

   G4Tubs *cpTube = new G4Tubs(cp.name,
	                       cp.r_min,cp.r_max,
	                       cp.length/2.,
	                       cp.startPhi,cp.dPhi);

   // torus geometry  
//    G4Torus *cnTube = new G4Torus(cn.name,
//                                  cn.r_min   ,cn.r_max,cn.r_tor,
//                                  cn.startPhi,cn.dPhi);
// 
//    G4Torus *cpTube = new G4Torus(cp.name,
// 	                         cp.r_min   ,cp.r_max,cp.r_tor,
// 	                         cp.startPhi,cp.dPhi);

   // create a union of the coils
   G4ThreeVector P      = G4ThreeVector(0.*cm,0.*cm,D);  // separated by z = D  
 
   G4UnionSolid *coils; 
   coils = new G4UnionSolid("coils",cnTube,cpTube,0,P); // no rotation

   // now for the logical volume 
   G4LogicalVolume *logicCoils = new G4LogicalVolume(coils,GetMaterial("Aluminum"),"logicCoils");
   logicCoils->SetVisAttributes(visCoil);  

   // NOTE: the position is the "origin" because we already rotated 
   //       the shell relative to the mother logical volume which is the physical world 
   G4ThreeVector Pcc = G4ThreeVector(0,0,0); 
   new G4PVPlacement(0,                // rotation [relative to mother]                                       
                     Pcc,              // position [relative to mother]                                            
                     logicCoils,       // logical volume                                      
                     coilName_pp,      // name                                                  
                     logicCoilShell,   // logical mother volume          
                     isBoolean,        // no boolean operations                        
                     0,                // copy number                                   
                     fCheckOverlaps);  // check overlaps                              

}
//______________________________________________________________________________
G4LogicalVolume *He3TargetDetectorConstruction::BuildGlassCell(){
   // construct the glass cell of the target
   // FIXME: Check transfer tube coordinates and dimensions (see JT file) 
   //        - Vertical tubes may be separated by 1.8 in = 2.286 cm  
   //        - Horizontal tubes may be higher in y than originally determined.  Posts will need to be taller  

   //---- pumping chamber ----
   partParameters_t pumpCh;
   GetPart("pumpingChamber",pumpCh); 

   G4Sphere *pumpChamberShape = new G4Sphere(pumpCh.name,
	                                     pumpCh.r_min     ,pumpCh.r_max,
	                                     pumpCh.startPhi  ,pumpCh.dPhi, 
	                                     pumpCh.startTheta,pumpCh.dTheta); 

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
	                                   tgtCh.startPhi ,tgtCh.dPhi); 

   G4ThreeVector P_tc = G4ThreeVector(tgtCh.x,tgtCh.y,tgtCh.z); 
   G4RotationMatrix *rm_tc = new G4RotationMatrix();
   rm_tc->rotateX(tgtCh.rx); 
   rm_tc->rotateY(tgtCh.ry); 
   rm_tc->rotateZ(tgtCh.rz); 

   // //---- end window on target chamber, downstream ----  

   // partParameters_t ewDn;
   // GetPart("endWindow_dn",ewDn); 

   // G4Sphere *endWindowShapeDn = new G4Sphere(ewDn.name,
   //       ewDn.r_min     ,ewDn.r_max,
   //       ewDn.startPhi  ,ewDn.dPhi, 
   //       ewDn.startTheta,ewDn.dTheta); 

   // G4ThreeVector P_ewd = G4ThreeVector(ewDn.x,ewDn.y,ewDn.z); 
   // G4RotationMatrix *rm_ewd = new G4RotationMatrix();
   // rm_ewd->rotateX(ewDn.rx); 
   // rm_ewd->rotateY(ewDn.ry); 
   // rm_ewd->rotateZ(ewDn.rz); 

   // //---- end window on target chamber, upstream ----  

   // partParameters_t ewUp;
   // GetPart("endWindow_up",ewUp); 

   // G4Sphere *endWindowShapeUp = new G4Sphere(ewUp.name,
   //       ewUp.r_min     ,ewUp.r_max,
   //       ewUp.startPhi  ,ewUp.dPhi, 
   //       ewUp.startTheta,ewUp.dTheta); 

   // G4ThreeVector P_ewu = G4ThreeVector(ewUp.x,ewUp.y,ewUp.z); 
   // G4RotationMatrix *rm_ewu = new G4RotationMatrix();
   // rm_ewu->rotateX(ewUp.rx); 
   // rm_ewu->rotateY(ewUp.ry); 
   // rm_ewu->rotateZ(ewUp.rz); 

   //---- transfer tube elbow, downstream 

   partParameters_t tted; 
   GetPart("transTubeEl_dn",tted); 

   G4Torus *transTubeElDnShape = new G4Torus(tted.name,
	                                     tted.r_min   ,tted.r_max,tted.r_tor,
	                                     tted.startPhi,tted.dPhi); 

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
	 tteu.startPhi  ,tteu.dPhi); 

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
	 ttedl.startPhi  ,ttedl.dPhi); 

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
	 tteul.startPhi  ,tteul.dPhi); 

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
	 tts.startPhi  ,tts.dPhi,
	 tts.startTheta,tts.dTheta); 

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
	 ttuz.startPhi ,ttuz.dPhi); 

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
	 ttuy.startPhi ,ttuy.dPhi); 

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
	 ttdz.startPhi ,ttdz.dPhi); 

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
	 ttdby.startPhi ,ttdby.dPhi); 

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
	 ttday.startPhi ,ttday.dPhi); 

   G4ThreeVector P_ttday = G4ThreeVector(ttday.x,ttday.y,ttday.z); 
   G4RotationMatrix *rm_ttday = new G4RotationMatrix();
   rm_ttday->rotateX(ttday.rx); 
   rm_ttday->rotateY(ttday.ry); 
   rm_ttday->rotateZ(ttday.rz); 

   //---- transfer tube post: downstream, along y ---- 
   partParameters_t ttpdy; 
   GetPart("transTubePost_dn",ttpdy);

   G4Tubs *transTubePostDnShape = new G4Tubs(ttpdy.name,
	 ttpdy.r_min    ,ttpdy.r_max,
	 ttpdy.length/2.,
	 ttpdy.startPhi ,ttpdy.dPhi); 

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
	 ttpuy.startPhi ,ttpuy.dPhi); 

   G4ThreeVector P_ttpuy = G4ThreeVector(ttpuy.x,ttpuy.y,ttpuy.z); 
   G4RotationMatrix *rm_ttpuy = new G4RotationMatrix();
   rm_ttpuy->rotateX(ttpuy.rx); 
   rm_ttpuy->rotateY(ttpuy.ry); 
   rm_ttpuy->rotateZ(ttpuy.rz); 

   // Create unions to make this a single continuous piece of material  
   G4UnionSolid *glassCell; // this is the top level object everything becomes 

   // build the target chamber + endcaps
   // // target chamber + upstream window 
   // glassCell = new G4UnionSolid("gc_tc_ewu",targetChamberShape,endWindowShapeUp,rm_ewu,P_ewu);
   // // downstream window  
   // glassCell = new G4UnionSolid("gc_tc_ewud",glassCell,endWindowShapeDn,rm_ewd,P_ewd);  

   // transfer tube posts
   // upstream  
   glassCell = new G4UnionSolid("gc_tc_ewud_pu" ,targetChamberShape,transTubePostUpShape,rm_ttpuy,P_ttpuy);  
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

   G4LogicalVolume *logicGlassCell = new G4LogicalVolume(glassCell,GetMaterial("GE180"),"logicGlassCell");

   return logicGlassCell; 
}
//______________________________________________________________________________
G4Material *He3TargetDetectorConstruction::GetMaterial(G4String name){

  std::map<G4String, G4Material*>::iterator it = fMaterialsMap.find( name );

  if( it != fMaterialsMap.end() ){
    return fMaterialsMap[name];
  } else {
    std::cout << "[He3TargetDetectorConstruction::GetMaterial]: ERROR Material " << name.data() 
              << " not found! " << std::endl;
    exit(1);
    return NULL;
  }

}
//______________________________________________________________________________
int He3TargetDetectorConstruction::ConstructMaterials(){

   // Get nist material manager
   G4NistManager* nist = G4NistManager::Instance();

   // world material 
   G4Material *G4Air = nist->FindOrBuildMaterial("G4_AIR");
   fMaterialsMap["G4air"] = G4Air; 

   // Define elements, compounds here that we'll need 
   // We don't need to supply the molar mass, because G4NistManager does it for us!
   G4int Z,N,ncomponents,nel; 
   G4double abundance;

   G4Isotope *iso_3He = new G4Isotope( "He3", Z=2, N=3 );

   G4Element *el3He = new G4Element("Helium3","3He",ncomponents=1); //Define isotopically pure Helium-3 
   el3He->AddIsotope( iso_3He, abundance=100.0*perCent );

   G4Element *elH  = nist->FindOrBuildElement("H"); 
   G4Element *elO  = nist->FindOrBuildElement("O");
   G4Element *elC  = nist->FindOrBuildElement("C");
   G4Element *elAl = nist->FindOrBuildElement("Al");
   G4Element *elSi = nist->FindOrBuildElement("Si");
   G4Element *elCa = nist->FindOrBuildElement("Ca");
   G4Element *elSr = nist->FindOrBuildElement("Sr");
   G4Element *elBa = nist->FindOrBuildElement("Ba");
   G4Element *elCu = nist->FindOrBuildElement("Cu");
   G4Element *elMn = nist->FindOrBuildElement("Mn"); 
   G4Element *elFe = nist->FindOrBuildElement("Fe"); 
   G4Element *elS  = nist->FindOrBuildElement("S"); 
   G4Element *elP  = nist->FindOrBuildElement("P"); 

   // Cu
   G4Material *Cu = new G4Material("Copper",8.96*g/cm3,1.); 
   Cu->AddElement(elCu,1);  
   fMaterialsMap["Copper"] = Cu;

   // Al 
   G4Material *Al = new G4Material("Aluminum",2.7*g/cm3,1.); 
   Al->AddElement(elAl,1); 
   fMaterialsMap["Aluminum"] = Al; 

   G4Material* NEMAG10 = new G4Material("NEMAG10",1.70*g/cm3,nel=4);
   NEMAG10->AddElement(elSi, 1);
   NEMAG10->AddElement(elO , 2);
   NEMAG10->AddElement(elC , 3);
   NEMAG10->AddElement(elH , 3);
   fMaterialsMap["NEMAG10"] = NEMAG10; 

   // GE180   
   G4double bigden = 1e9*g/cm3; // why so big? To make these materials not weigh as much in the physics?  
   // gather necessary molecules and compounds  
   // SiO2 60.3%
   G4Material* SiO2 = new G4Material("GE180_SiO2", 2.2*g/cm3, 2 );
   SiO2->AddElement(elSi, 1);
   SiO2->AddElement(elO, 2);
   fMaterialsMap["GE180_SiO2"] = SiO2;
   // BaO  18.2%
   G4Material* BaO = new G4Material("GE180_BaO", bigden, 2 );
   BaO->AddElement(elBa, 1);
   BaO->AddElement(elO, 1);
   fMaterialsMap["GE180_BaO"] = BaO;
   // Al2O3 14.3%
   G4Material* Al2O3 = new G4Material("GE180_Al2O3", bigden, 2 );
   Al2O3->AddElement(elAl, 2);
   Al2O3->AddElement(elO, 3);
   fMaterialsMap["GE180_Al2O3"] = Al2O3;
   // CaO   6.5%
   G4Material* CaO = new G4Material("GE180_CaO", bigden, 2 );
   CaO->AddElement(elCa, 1);
   CaO->AddElement(elO, 1);
   fMaterialsMap["GE180_CaO"] = CaO;
   // SrO   0.25%
   G4Material* SrO = new G4Material("GE180_SrO", bigden, 2 );
   SrO->AddElement(elSr, 1);
   SrO->AddElement(elO, 1);
   fMaterialsMap["GE180_SrO"] = SrO;

   // Density 2.76 g/cm^3
   // Index of Refraction 1.536
   G4Material* GE180 = new G4Material("GE180", 2.76*g/cm3, 5);
   GE180->AddMaterial(SiO2 , 0.6039);
   GE180->AddMaterial(BaO  , 0.1829);
   GE180->AddMaterial(Al2O3, 0.1439);
   GE180->AddMaterial(CaO  , 0.0659);
   GE180->AddMaterial(SrO  , 0.0034);
   fMaterialsMap["GE180"] = GE180;

   //---- polarized 3He ----
   G4double gasden = 10.77*atmosphere*(3.016*g/Avogadro)/(300*kelvin*k_Boltzmann);
   G4Material *pol3He = new G4Material("pol3He", gasden, 1 );
   pol3He->AddElement(el3He, 1);
   fMaterialsMap["He3"] = pol3He;  

   // for the support materials 
   G4Material *Steel = nist->FindOrBuildMaterial("G4_STAINLESS-STEEL");
   fMaterialsMap["Stainless_Steel"] = Steel; 

   // AISI 1008 carbon steel
   // details from https://www.azom.com/article.aspx?ArticleID=6538
   // assuming the composition with 0.06% C, 0.38% Mn, 0.01% Si for a density of 7.782*g/cm3   
   G4Material *Carbon_Steel_1008 = new G4Material("Carbon_Steel_1008",7.872*g/cm3,6); 
   Carbon_Steel_1008->AddElement(elFe,0.9905); // avg of 0.9931 and 0.9970  
   Carbon_Steel_1008->AddElement(elMn,0.0038);  
   Carbon_Steel_1008->AddElement(elC ,0.0006); 
   Carbon_Steel_1008->AddElement(elS ,0.0005); 
   Carbon_Steel_1008->AddElement(elP ,0.0004); 
   Carbon_Steel_1008->AddElement(elSi,0.0001); 
   fMaterialsMap["Carbon_Steel_1008"] = Carbon_Steel_1008; 

   // print for confirmation
   if(fDebug){
      std::cout << "Constructed materials: " << std::endl; 
      for(auto it=fMaterialsMap.cbegin(); it!=fMaterialsMap.cend(); ++it){
	 std::cout << (*it).first << ": " << (*it).second << std::endl;
      }
   }

   return 0;
}
//______________________________________________________________________________
int He3TargetDetectorConstruction::GetPart(const std::string partName,partParameters_t &data){
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
      std::cout << " THETA start theta: " << data.startTheta << " d theta: " << data.dTheta << std::endl; 
      std::cout << " PHI start phi: " << data.startPhi   << " d phi: "   << data.dPhi << std::endl; 
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
      exit(1);
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
   double x_len=0,y_len=0,z_len=0;
   double length=0,startTheta=0,dTheta=0,startPhi=0,dPhi=0;
   double x=0,y=0,z=0,rx=0,ry=0,rz=0;
   double LEN_UNIT=1.;
   double DEG_UNIT=1.;

   // now parse the data
   int rc=0;
   for(int i=0;i<NROW;i++){
      // split the row into a vector which represents the columns 
      rc = SplitString(',',row[i],col);
      if(rc!=0){
         G4cout << "[He3TargetDetectorConstruction::ReadData]: Cannot parse string " << row[i] << G4endl;
	 exit(1);
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
      x_len           = std::atof( col[8].c_str()  );  
      y_len           = std::atof( col[9].c_str()  );  
      z_len           = std::atof( col[10].c_str() );  
      startTheta      = std::atof( col[11].c_str() );  
      dTheta          = std::atof( col[12].c_str() );  
      startPhi        = std::atof( col[13].c_str() );  
      dPhi            = std::atof( col[14].c_str() );  
      x               = std::atof( col[15].c_str() );  
      y               = std::atof( col[16].c_str() );  
      z               = std::atof( col[17].c_str() );  
      rx              = std::atof( col[18].c_str() );  
      ry              = std::atof( col[19].c_str() );  
      rz              = std::atof( col[20].c_str() );
      // convert to units
      if( dataPt.len_unit.compare("mm")==0 )  LEN_UNIT = mm;
      if( dataPt.len_unit.compare("cm")==0)   LEN_UNIT = cm; 
      if( dataPt.len_unit.compare("m")==0)    LEN_UNIT = m; 
      if( dataPt.len_unit.compare("in")==0)   LEN_UNIT = 25.4*mm; 
      if( dataPt.ang_unit.compare("deg")==0 ) DEG_UNIT = deg;
      if( dataPt.ang_unit.compare("rad")==0)  DEG_UNIT = rad; 
      // store data 
      dataPt.r_tor      = LEN_UNIT*r_tor;  
      dataPt.r_min      = LEN_UNIT*r_min;  
      dataPt.r_max      = LEN_UNIT*r_max;  
      dataPt.length     = LEN_UNIT*length; 
      dataPt.x_len      = LEN_UNIT*x_len;  
      dataPt.y_len      = LEN_UNIT*y_len;  
      dataPt.z_len      = LEN_UNIT*z_len;  
      dataPt.startTheta = DEG_UNIT*startTheta;  
      dataPt.dTheta     = DEG_UNIT*dTheta;  
      dataPt.startPhi   = DEG_UNIT*startPhi;  
      dataPt.dPhi       = DEG_UNIT*dPhi;  
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
	     << " start theta: "          << dataPt.startTheta << " dtheta: " << dataPt.dTheta 
	     << " start phi: "            << dataPt.startPhi   << " dphi: "   << dataPt.dPhi 
	     << " x: "                    << dataPt.x          << " y: "      << dataPt.y 
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
int He3TargetDetectorConstruction::SplitString(const char delim,const std::string inStr,
                                               std::vector<std::string> &out){
   // split a string by a delimiter
   std::stringstream ss(inStr);
   while( ss.good() ){
      std::string substr;
      std::getline(ss,substr,delim);
      out.push_back(substr);
   }
   return 0;
}
