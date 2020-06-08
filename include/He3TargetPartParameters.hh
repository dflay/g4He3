#ifndef HE3_PART_PARAMETERS_HH
#define HE3_PART_PARAMETERS_HH

// data struct to read in geometry dimensions, positional and rotational information

typedef struct partParameters { 
   std::string name;     // part name     
   std::string shape;    // shape (sphere, torus, tube, etc)    
   std::string len_unit; // length units  
   std::string ang_unit; // angular units  
   double r_tor;         // major radius (torus only)    
   double r_min;         // minimum radius (tube, sphere, torus)   
   double r_max;         // maximum radius (tube, sphere, torus)     
   double length;        // full length of object (tube only)  
   double startTheta;    // extruding polar angle (about y)  
   double dTheta;        // polar angle step size 
   double startPhi;      // extruding azimuthal angle (about z)  
   double dPhi;          // azimuthal angle step size
   double x;             // x position relative to mother volume 
   double y;             // y position relative to mother volume  
   double z;             // z position relative to mother volume  
   double rx;            // rotation about x axis relative to mother volume    
   double ry;            // rotation about y axis relative to mother volume   
   double rz;            // rotation about z axis relative to mother volume  

   // constructor 
   partParameters(): 
      name("none"),shape("none"),len_unit("none"),ang_unit("none"),
      r_tor(0),r_min(0),r_max(0),length(0),startTheta(0),dTheta(0),startPhi(0),dPhi(0),
      x(0),y(0),z(0),rx(0),ry(0),rz(0)
   { }
 
} partParameters_t;   

#endif 
