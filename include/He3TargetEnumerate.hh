#ifndef HE3_ENUMERATE_HH
#define HE3_ENUMERATE_HH

// switches for Helmholtz coils or shielding 
// - 146  => Q2 = 1.46  (GeV/c)^2  
// - 368  => Q2 = 3.68  (GeV/c)^2  
// - 677  => Q2 = 6.77  (GeV/c)^2  
// - 1018 => Q2 = 10.18 (GeV/c)^2
// for shielding:  
// - full => full window cut (remove panel 1, 2, 3, and door) 
// - new  => new design from Bert Metzger (6/2020)  

enum expConf { 
   kSBS_GEN_146  = 1,
   kSBS_GEN_368  = 2,
   kSBS_GEN_677  = 3,
   kSBS_GEN_1018 = 4,
   kSBS_GEN_full = 5,
   kSBS_GEN_new  = 6
}; 

#endif 
