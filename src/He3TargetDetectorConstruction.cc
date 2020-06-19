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
  fScoringVolume(0),fDebug(false),fCheckOverlaps(true)
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
  BuildHelmholtzCoils("maj",logicWorld);  
  BuildHelmholtzCoils("rfy",logicWorld);  
  BuildHelmholtzCoils("min",logicWorld); 

  // shield 
  BuildShield(logicWorld); 

  // support ladder 
  BuildLadderPlate(logicWorld); 

  // always return the physical World
  return physWorld;
}
//______________________________________________________________________________
void He3TargetDetectorConstruction::BuildPolarizedHe3(G4LogicalVolume *logicMother){
  // build the polarized 3He logical volume and place it in the logicMother 
  // TODO: Check volume overlap with glass cell  

  //---- target chamber component  
  partParameters_t tc; 
  GetPart("targetChamber",tc); 
  tc.r_max = tc.r_min; 
  tc.r_min = 0.*cm;
 
  G4Tubs *tcShape = new G4Tubs("He3_tc",
                               tc.r_min    ,tc.r_max,
                               tc.length/2.,
                               tc.startPhi ,tc.dPhi);

  // end window (upstream)
  // combine info from lip, rlip, end cap to one object here
  // try to make this a simple union of a cylinder and a hemisphere.  
  // TODO: may need something more involved later... 

  // build the "base" -- combines end window main shaft + lip   
  partParameters_t mshu;
  GetPart("ew_mainShaft_up",mshu);

  partParameters_t lipu;
  GetPart("ew_lip_up",lipu);

  partParameters_t baseu; 
  baseu.name       = "He3_base_u"; 
  baseu.r_max      = tc.r_max; // mshu.r_min; 
  baseu.r_min      = 0.*cm;
  baseu.length     = mshu.length + lipu.length;  
  baseu.x          = 0.*cm; 
  baseu.y          = 0.*cm; 
  baseu.z          = (-1.)*( 0.5*tc.length + 0.5*baseu.length ); 
  baseu.startTheta = 0.*deg; 
  baseu.dTheta     = 0.*deg; 
  baseu.startPhi   = 0.*deg; 
  baseu.dPhi       = 360.*deg; 

  G4Tubs *base_u = new G4Tubs(baseu.name,
	                      baseu.r_min    ,baseu.r_max,
	                      baseu.length/2.,
	                      baseu.startPhi ,baseu.dPhi);

  G4ThreeVector P_baseu      = G4ThreeVector(baseu.x,baseu.y,baseu.z);  
  G4RotationMatrix *rm_baseu = new G4RotationMatrix();
  rm_baseu->rotateX(baseu.rx);   
  rm_baseu->rotateY(baseu.ry);   
  rm_baseu->rotateZ(baseu.rz);   
  
  // hemisphere
  partParameters_t hspu; 
  hspu.name       = "He3_hspu"; 
  hspu.r_max      = tc.r_max;  
  hspu.r_min      = 0.*cm;
  hspu.length     = 0.*cm;  
  hspu.x          = 0.*cm; 
  hspu.y          = 0.*cm; 
  hspu.z          = baseu.z - 0.5*hspu.r_max; 
  hspu.rx         = 0.*deg; 
  hspu.ry         = 180.*deg; 
  hspu.rz         = 0.*deg; 
  hspu.startTheta = 0.*deg; 
  hspu.dTheta     = 90.*deg; 
  hspu.startPhi   = 0.*deg; 
  hspu.dPhi       = 360.*deg; 

  G4Sphere *hSphere_u = new G4Sphere(hspu.name,
                                     hspu.r_min      ,hspu.r_max,
                                     hspu.startPhi   ,hspu.dPhi,
                                     hspu.startTheta ,hspu.dTheta);

  G4ThreeVector P_hspu      = G4ThreeVector(hspu.x,hspu.y,hspu.z);  
  G4RotationMatrix *rm_hspu = new G4RotationMatrix();
  rm_hspu->rotateX(hspu.rx);   
  rm_hspu->rotateY(hspu.ry);   
  rm_hspu->rotateZ(hspu.rz);   

  // downstream 
  partParameters_t based = baseu; 
  based.name = "He3_base_d"; 
  based.z *= (-1.); 

  G4Tubs *base_d = new G4Tubs(based.name,
 	                      based.r_min    ,based.r_max,
	                      based.length/2.,
	                      based.startPhi ,based.dPhi);

  G4ThreeVector P_based      = G4ThreeVector(based.x,based.y,based.z);  
  G4RotationMatrix *rm_based = new G4RotationMatrix();
  rm_based->rotateX(based.rx);   
  rm_based->rotateY(based.ry);   
  rm_based->rotateZ(based.rz);   

  partParameters_t hspd = hspu;
  hspd.name  = "He3_hspd";  
  hspd.z    *= (-1.); 
  hspd.ry    = 0.*deg;

  G4Sphere *hSphere_d = new G4Sphere(hspd.name,
                                     hspd.r_min      ,hspd.r_max,
                                     hspd.startPhi   ,hspd.dPhi,
                                     hspd.startTheta ,hspd.dTheta);

  G4ThreeVector P_hspd      = G4ThreeVector(hspd.x,hspd.y,hspd.z);  
  G4RotationMatrix *rm_hspd = new G4RotationMatrix();
  rm_hspd->rotateX(hspd.rx);   
  rm_hspd->rotateY(hspd.ry);   
  rm_hspd->rotateZ(hspd.rz);   

  // Union solid 
  // use same rotation and positional vectors as glass shell!  
  G4UnionSolid *he3Tube;
  // target chamber + base upstream  
  he3Tube = new G4UnionSolid("tc_bu"     ,tcShape,base_u   ,rm_baseu,P_baseu);
  // assembly + hsphere upstream 
  he3Tube = new G4UnionSolid("tc_bhu"    ,he3Tube,hSphere_u,rm_hspu,P_hspu);  
  // assembly + base downstream 
  he3Tube = new G4UnionSolid("tcbhu_hspd",he3Tube,base_d   ,rm_based,P_based);  
  // assembly + hsphere downstream 
  he3Tube = new G4UnionSolid("heTube"    ,he3Tube,hSphere_d,rm_hspd,P_hspd);  

  // logical volume of He3
  G4LogicalVolume *logicHe3 = new G4LogicalVolume(he3Tube,GetMaterial("He3"),"logicHe3");  

  // set the color of He3 
  G4VisAttributes *visHe3 = new G4VisAttributes(); 
  visHe3->SetColour( G4Colour::Yellow() );
  visHe3->SetForceWireframe(true);  
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
void He3TargetDetectorConstruction::BuildShield(G4LogicalVolume *logicMother){
   // shield box for the target magnetic field 
   // - material: 1008 steel
   // - the shield is actually two layers
   //   - each layer is 0.25" thick
   //   - outer surfaces of layers separated by 1.29" 
   //     => 1.29 - 0.25 = 1.04" center-to-center distance 
   // - Layer 1 = outside layer, Layer 2 = inside layer 
   // - drawing number: A09016-03-05-0000 
   
   G4double gap = 1.04*2.54*cm; // 0.52" center-to-center distance 

   // offset is somewhat arbitrary, choosing to center the entire enclosure based on its full extended length 
   G4double delta_x = 51.21*cm; // 24.19*25.4*mm - 9.5*cm;    // offset of whole enclosure in x direction  
   G4double delta_y = 0.*mm;                     // offset of whole enclosure in y direction  
   G4double delta_z = 0.*mm;                     // offset of whole enclosure in z direction 

   //---- shield door along x, upstream 
   partParameters_t sdu1;  
   GetPart("shield_door_up",sdu1);
   partParameters_t sdu2 = sdu1; 

   // set up for centering device on the midpoint between the panels 
   sdu1.z -= gap/2.; 
   sdu2.z += gap/2.;

   // // layer 1
   // G4Box *shieldDoor_up1     = new G4Box("sdu1",sdu1.x_len/2.,sdu1.y_len/2.,sdu1.z_len/2.);
   // G4ThreeVector P_sdu1      = G4ThreeVector(sdu1.x,sdu1.y,sdu1.z);  
   // G4RotationMatrix *rm_sdu1 = new G4RotationMatrix(); 
   // rm_sdu1->rotateX(sdu1.rx); rm_sdu1->rotateY(sdu1.ry); rm_sdu1->rotateZ(sdu1.rz);   

   // // layer 2 
   // G4Box *shieldDoor_up2     = new G4Box("sdu2",sdu2.x_len/2.,sdu2.y_len/2.,sdu2.z_len/2.);
   // G4ThreeVector P_sdu2      = G4ThreeVector(sdu2.x,sdu2.y,sdu2.z);  
   // G4RotationMatrix *rm_sdu2 = new G4RotationMatrix(); 
   // rm_sdu2->rotateX(sdu2.rx); rm_sdu2->rotateY(sdu2.ry); rm_sdu2->rotateZ(sdu2.rz);   
   
   //---- shield pane along x, upstream  
   partParameters_t spd1;  
   GetPart("shield_pane_dn",spd1);
   partParameters_t spd2 = spd1; 

   // finish definition of structure offset 
   delta_z = spd1.z;  

   // note the sign inversion here! 
   spd1.z += gap/2.; 
   spd2.z -= gap/2.; 

   // layer 1
   G4Box *shieldPane_dn1     = new G4Box("spd1",spd1.x_len/2.,spd1.y_len/2.,spd1.z_len/2.);
   G4ThreeVector P_spd1      = G4ThreeVector(spd1.x,spd1.y,spd1.z);  
   G4RotationMatrix *rm_spd1 = new G4RotationMatrix(); 
   rm_spd1->rotateX(spd1.rx); rm_spd1->rotateY(spd1.ry); rm_spd1->rotateZ(spd1.rz);   

   // layer 2 
   G4Box *shieldPane_dn2     = new G4Box("spd2",spd2.x_len/2.,spd2.y_len/2.,spd2.z_len/2.);
   G4ThreeVector P_spd2      = G4ThreeVector(spd2.x,spd2.y,spd2.z);  
   G4RotationMatrix *rm_spd2 = new G4RotationMatrix(); 
   rm_spd2->rotateX(spd2.rx); rm_spd2->rotateY(spd2.ry); rm_spd2->rotateZ(spd2.rz);

   //---- shield pane, beam left [along z]  
   partParameters_t spl1;  
   GetPart("shield_pane_zl",spl1);
   partParameters_t spl2 = spl1; 

   // adjust in x direction because we are using the beam right pane as our origin 
   spl1.x -= delta_x; 
   spl2.x -= delta_x; 
   // adjust for gap between layers 
   spl1.x += gap/2.; 
   spl2.x -= gap/2.; 

   // layer 1
   // beam-left parts, panes along z 
   G4Box *shieldPane_l1      = new G4Box("spl1",spl1.x_len/2.,spl1.y_len/2.,spl1.z_len/2.);

   G4ThreeVector P_spl1      = G4ThreeVector(spl1.x,spl1.y,spl1.z);  
   G4RotationMatrix *rm_spl1 = new G4RotationMatrix(); 
   rm_spl1->rotateX(spl1.rx); rm_spl1->rotateY(spl1.ry); rm_spl1->rotateZ(spl1.rz); 

   // layer 2 
   G4Box *shieldPane_l2      = new G4Box("spl2",spl2.x_len/2.,spl2.y_len/2.,spl2.z_len/2.);

   G4ThreeVector P_spl2      = G4ThreeVector(spl2.x,spl2.y,spl2.z);  
   G4RotationMatrix *rm_spl2 = new G4RotationMatrix(); 
   rm_spl2->rotateX(spl2.rx); rm_spl2->rotateY(spl2.ry); rm_spl2->rotateZ(spl2.rz);

   // beam-right parts

   //---- shield pane, beam left [along z]  
   partParameters_t spr1;  
   GetPart("shield_pane_zr",spr1);
   partParameters_t spr2 = spr1; 

   // adjust in x direction because we are using the beam right pane as our origin 
   spr1.x -= delta_x; 
   spr2.x -= delta_x; 

   // adjust for gap between layers 
   spr1.x += gap/2.; 
   spr2.x -= gap/2.; 

   // std::cout << "spr1 beam right (x,y,z) = " << spr1.x/cm << "," << spr1.y/cm << "," << spr1.z/cm << std::endl;
   // std::cout << "spr2 beam right (x,y,z) = " << spr2.x/cm << "," << spr2.y/cm << "," << spr2.z/cm << std::endl;

   G4Box *shieldPane_r1      = new G4Box("spr1",spr1.x_len/2.,spr1.y_len/2.,spr1.z_len/2.); 
   G4ThreeVector P_spr1      = G4ThreeVector(spr1.x,spr1.y,spr1.z); 
   G4RotationMatrix *rm_spr1 = new G4RotationMatrix(); 
   rm_spr1->rotateX(spr1.rx); rm_spr1->rotateY(spr1.ry); rm_spr1->rotateZ(spr1.rz);
 
   G4Box *shieldPane_r2      = new G4Box("spr2",spr2.x_len/2.,spr2.y_len/2.,spr2.z_len/2.);  
   G4ThreeVector P_spr2      = G4ThreeVector(spr2.x,spr2.y,spr2.z); 
   G4RotationMatrix *rm_spr2 = new G4RotationMatrix(); 
   rm_spr2->rotateX(spr2.rx); rm_spr2->rotateY(spr2.ry); rm_spr2->rotateZ(spr2.rz);

   //---- shield pane, upstream [along x]
   partParameters_t spu1; 
   GetPart("shield_pane_up",spu1); 
   partParameters_t spu2 = spu1; 

   // adjust for gap  
   spu1.z += gap/2.; 
   spu2.z -= gap/2.;
   // adjust for z positioning relative to downstream pane 
   spu1.z -= delta_z; 
   spu2.z -= delta_z; 

   // layer 1
   G4Box *shieldPane_up1     = new G4Box("spu1",spu1.x_len/2.,spu1.y_len/2.,spu1.z_len/2.);
   G4ThreeVector P_spu1      = G4ThreeVector(spu1.x,spu1.y,spu1.z);  
   G4RotationMatrix *rm_spu1 = new G4RotationMatrix(); 
   rm_spu1->rotateX(spu1.rx); rm_spu1->rotateY(spu1.ry); rm_spd1->rotateZ(spu1.rz);   

   // layer 2 
   G4Box *shieldPane_up2     = new G4Box("spu2",spu2.x_len/2.,spu2.y_len/2.,spu2.z_len/2.);
   G4ThreeVector P_spu2      = G4ThreeVector(spu2.x,spu2.y,spu2.z);  
   G4RotationMatrix *rm_spu2 = new G4RotationMatrix(); 
   rm_spu2->rotateX(spu2.rx); rm_spu2->rotateY(spu2.ry); rm_spu2->rotateZ(spu2.rz);

   // std::cout << "spu1 beam right (x,y,z) = " << spu1.x/cm << "," << spu1.y/cm << "," << spu1.z/cm << std::endl;
   // std::cout << "spu2 beam right (x,y,z) = " << spu2.x/cm << "," << spu2.y/cm << "," << spu2.z/cm << std::endl;

   // collect into single object
   G4RotationMatrix *rm0 = new G4RotationMatrix(); 
   rm0->rotateX(0); rm0->rotateY(0); rm0->rotateZ(0); 

   // build the object along the x axis.  Then center and rotate it 
   G4UnionSolid *enclosure; 
   // [downstream] pane inner and outer  
   G4ThreeVector P_pd2     = G4ThreeVector(0.,0,gap);                              // relative to pane 1  
   enclosure = new G4UnionSolid("sd_p12",shieldPane_dn1,shieldPane_dn2,0,P_pd2);  
   // // [downstream] door, inner panel 
   // G4ThreeVector P_dd1 = G4ThreeVector(-(sdu1.x_len/2.+spd1.x_len/2.),0,0);    
   // enclosure = new G4UnionSolid("sd_p12_d1",enclosure,shieldDoor_up1,0,P_dd1);   
   // // [downstream] door, outer panel   
   // G4ThreeVector P_dd2 = G4ThreeVector(-(sdu1.x_len/2.+spd1.x_len/2.),0,gap);   
   // enclosure = new G4UnionSolid("sd_p12_d12",enclosure,shieldDoor_up2,0,P_dd2);   
   // [beam left] pane, along z  
   enclosure = new G4UnionSolid("sd_p12_d12_l1" ,enclosure,shieldPane_l1,rm_spl1,P_spl1);   
   enclosure = new G4UnionSolid("sd_p12_d12_l12",enclosure,shieldPane_l2,rm_spl2,P_spl2);   
   // [beam right] pane, along z  
   enclosure = new G4UnionSolid("sd_p12_d12_l12_r1" ,enclosure,shieldPane_r1,rm_spr1,P_spr1);   
   enclosure = new G4UnionSolid("sd_p12_d12_l12_r12",enclosure,shieldPane_r2,rm_spr2,P_spr2);  
   // [upstream] pane, inner  
   enclosure = new G4UnionSolid("sd_p12_d12_l12_r12_u1",enclosure,shieldPane_up1,rm_spu1,P_spu1);  
   // [upstream] pane, outer  
   enclosure = new G4UnionSolid("shield"               ,enclosure,shieldPane_up2,rm_spu2,P_spu2);  

   // visualization 
   G4VisAttributes *vis = new G4VisAttributes();
   // vis->SetForceWireframe();
   vis->SetColour( G4Colour::Magenta() ); 

   // logical volume 
   G4LogicalVolume *logicEnclosure = new G4LogicalVolume(enclosure,GetMaterial("Stainless_Steel"),"logicEnclosure");
   logicEnclosure->SetVisAttributes(vis);  

   // adjust for y rotation.  
   // FIXME: I don't like that this doesn't follow the rotated coordinates...
   // FIXME: Don't know the actual rotation angle!
   G4double dry     = -32.0*deg;          // 43.5? rotation angle of whole unit about y axis 
   G4double RY      = dry;
   G4double COS     = cos(RY); 
   G4double SIN     = sin(RY);

   G4double len    = spl1.x_len - delta_x; 
   G4double dx     =  (len/2.)*fabs(COS);
   G4double dz     = -(len/2.)*fabs(SIN);
 
   // place it.  these coordinates center the shield on the target
   // FIXME: Can we remove the fudge factor? 
   G4double shield_x = delta_x + dx - 15.*cm;    
   G4double shield_y = delta_y;     
   G4double shield_z = delta_z + dz;  

   G4ThreeVector P_s      = G4ThreeVector(shield_x,shield_y,shield_z); 
   G4RotationMatrix *rm_s = new G4RotationMatrix();
   rm_s->rotateX(0.*deg); rm_s->rotateY(RY); rm_s->rotateZ(0.*deg); 
 
   new G4PVPlacement(rm_s,                // rotation [relative to mother]    
                     P_s,                 // position [relative to mother] 
                     logicEnclosure,      // logical volume    
                     "physEnclosure",     // name                          
                     logicMother,         // logical mother volume            
                     false,               // no boolean operations          
                     0,                   // copy number                   
                     fCheckOverlaps);     // check overlaps               
 
}
//______________________________________________________________________________
void He3TargetDetectorConstruction::BuildLadderPlate(G4LogicalVolume *logicMother){
   // ladder plate  
   // - this comes close to the beam path
   // - drawing number: A09016-03-04-0601

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

   G4double lx = lvu.x; // FIXME: 5 cm beam right
   G4double ly = lvu.y; 
   G4double lz = lvu.z; 
    
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
void He3TargetDetectorConstruction::BuildHelmholtzCoils(const std::string type,G4LogicalVolume *logicMother){
   // Helmholtz coils for B fields 
   // - types: maj = large radius coil pair
   //          min = small radius coil pair 
   //          rfy = RF coil pair, aligned along the vertical (y) axis  
   // - distance between coils D = 0.5(rmin+rmax), roughly the major radius of the tube   
   //   - coil 1 (placed at -D/2) 
   //   - coil 2 (placed at +D/2)
   // drawing number: A09016-03-08-0000

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

   partParameters_t cn,cp;
   GetPart(partName,cn); 
   cp = cn;
   cn.name = coilName_n;  
   cp.name = coilName_p;  

   // coil parameters   
   // G4double D      = 0.5*(cn.r_min + cn.r_max);          // helmholtz separation D = R = 0.5(rmin + rmax) 
   G4double D      = cn.r_tor;          // helmholtz separation D = R = 0.5(rmin + rmax) 
   G4double shWall = 0;  

   if( type.compare("maj")==0 ) shWall = 5.0*mm;         // FIXME: Estimates for now! 
   if( type.compare("min")==0 ) shWall = 5.0*mm;         // FIXME: Estimates for now! 
   if( type.compare("rfy")==0 ) shWall = 0.030*25.4*mm; 

   // cylindrical geometry 
   // create the shell first 
//    partParameters_t cns;
//    cns.name     = shellName_n; 
//    cns.r_min    = cn.r_min - shWall;  
//    cns.r_max    = cn.r_max + shWall; 
//    cns.length   = cn.length;
//    cns.startPhi = 0.*deg;
//    cns.dPhi     = 360.*deg;   
//  
//    G4Tubs *cnsTube = new G4Tubs(cns.name,
//                                 cns.r_min    ,cns.r_max,
//                                 cns.length/2.,
//                                 cns.startPhi ,cns.dPhi);
// 
//    partParameters_t cps;
//    cps.name     = shellName_p; 
//    cps.r_min    = cp.r_min - shWall;  
//    cps.r_max    = cp.r_max + shWall; 
//    cps.length   = cp.length;
//    cps.startPhi = 0.*deg;
//    cps.dPhi     = 360.*deg;   
// 
//    G4Tubs *cpsTube = new G4Tubs(cps.name,
// 	                        cps.r_min    ,cps.r_max,
// 	                        cps.length/2.,
// 	                        cps.startPhi ,cps.dPhi);

   // torus geometry
   // - some adjustments here: we should make the *thickness* of the core 
   //   equal to the thickness specified in the drawings. 
   //   this is better than the estimate of the derived r_max of the torus centered on 
   //   the radius of curvature  
   cn.r_max = cn.length/2.;    
   cp.r_max = cp.length/2.;    

   // create the shell first 
   partParameters_t cns;
   cns.name     = shellName_n; 
   cns.r_min    = cn.r_max;  
   cns.r_max    = cn.r_max + shWall;
   cns.r_tor    = cn.r_tor;  
   cns.length   = cn.length;
   cns.startPhi = 0.*deg;
   cns.dPhi     = 360.*deg;   
 
   G4Torus *cnsTube = new G4Torus(cns.name,
                                  cns.r_min   ,cns.r_max,cns.r_tor,
                                  cns.startPhi,cns.dPhi);

   partParameters_t cps;
   cps.name     = shellName_p; 
   cps.r_min    = cp.r_max;  
   cps.r_max    = cp.r_max + shWall;
   cps.r_tor    = cp.r_tor;  
   cps.length   = cp.length;
   cps.startPhi = 0.*deg;
   cps.dPhi     = 360.*deg;   

   G4Torus *cpsTube = new G4Torus(cps.name,
	                          cps.r_min   ,cps.r_max,cps.r_tor,
	                          cps.startPhi,cps.dPhi);

   // create a union of the coils
   G4ThreeVector Ps = G4ThreeVector(0.*cm,0.*cm,D);  // separated by z = D  
   G4UnionSolid *coilShells; 
   coilShells = new G4UnionSolid("coilShells",cnsTube,cpsTube,0,Ps); // no rotation

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
   if( type.compare("maj")==0 || type.compare("min")==0 ) dry = 43.5*deg;  
 
   // rotation about y 
   // x' =  xcos + zsin 
   // z' = -xsin + zcos
 
   G4double ph = cn.ry + dry;  // total angular rotation!  
   G4double COS = cos(ph); 
   G4double SIN = sin(ph);

   // adjust for y rotation.  
   // FIXME: I don't like that this doesn't follow the rotated coordinates...  
   if(type.compare("maj")==0){ 
      x += (D/2)*fabs(COS); 
      z += (D/2)*fabs(SIN);
   }else if(type.compare("rfy")==0){
      y -= D/2.;
   }else if(type.compare("min")==0){
      x += (D/2)*fabs(COS); 
      z -= (D/2)*fabs(SIN);
   } 
 
   G4ThreeVector P_c    = G4ThreeVector(x,y,z);   
   G4RotationMatrix *rm = new G4RotationMatrix();
   rm->rotateX(cn.rx+drx); 
   rm->rotateY(cn.ry+dry); 
   rm->rotateZ(cn.rz+drz);

   char coilShellName[200];
   sprintf(coilShellName,"%s_shell",partName);  

   new G4PVPlacement(rm,               // rotation                                        
                     P_c,              // position                                             
                     logicCoilShell,   // logical volume                                      
                     coilShellName,    // name                                                  
                     logicMother,      // logical mother volume is the target chamber          
                     false,            // no boolean operations                        
                     0,                // copy number                                   
                     fCheckOverlaps);  // check overlaps                              

   // copper coil -- goes *inside* the shell  
   G4VisAttributes *visCoil = new G4VisAttributes();
   visCoil->SetColour( G4Colour(255,140,0) );  // dark orange 

   // cylindrical geometry 
//    G4Tubs *cnTube = new G4Tubs(cn.name,
//                                cn.r_min,cn.r_max,
//                                cn.length/2.,
//                                cn.startPhi,cn.dPhi);
// 
//    G4Tubs *cpTube = new G4Tubs(cp.name,
// 	                       cp.r_min,cp.r_max,
// 	                       cp.length/2.,
// 	                       cp.startPhi,cp.dPhi);
 
   // torus geometry  
   G4Torus *cnTube = new G4Torus(cn.name,
                                 cn.r_min   ,cn.r_max,cn.r_tor,
                                 cn.startPhi,cn.dPhi);

   G4Torus *cpTube = new G4Torus(cp.name,
	                         cp.r_min   ,cp.r_max,cp.r_tor,
	                         cp.startPhi,cp.dPhi);

   // create a union of the coils
   G4ThreeVector P      = G4ThreeVector(0.*cm,0.*cm,D);  // separated by z = D  
 
   G4UnionSolid *coils; 
   coils = new G4UnionSolid("coils",cnTube,cpTube,0,P); // no rotation

   // now for the logical volume 
   G4LogicalVolume *logicCoils = new G4LogicalVolume(coils,GetMaterial("Copper"),"logicCoils");
   logicCoils->SetVisAttributes(visCoil);  

   // NOTE: the position is the "origin" because we already rotated 
   //       the shell relative to the mother logical volume which is the physical world 
   G4ThreeVector Pcc = G4ThreeVector(0,0,0); 
   new G4PVPlacement(0,                // rotation [relative to mother]                                       
                     Pcc,              // position [relative to mother]                                            
                     logicCoils,       // logical volume                                      
                     coilName_pp,      // name                                                  
                     logicCoilShell,   // logical mother volume          
                     false,            // no boolean operations                        
                     0,                // copy number                                   
                     fCheckOverlaps);  // check overlaps                              

}
//______________________________________________________________________________
G4LogicalVolume *He3TargetDetectorConstruction::BuildGlassCell(){
   // construct the glass cell of the target 

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
