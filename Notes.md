# Notes 

# Helmholtz Coils

## Materials

Choose the main material to be aluminum.  Assume copper wiring is a small effect
 
## Cylindrical Geometry Dimensions 

Dimensions of coil geometry, implemented as cylinders (G4Tubs).  G10 wall is actually a G4Tubs 
with a hollow core, and a wall thickness as specified below    

| Label | rmin (cm) | rmax (cm) | length (cm) | G10 Wall (cm) |  
|-------|-----------|-----------|-------------|---------------|
| major | 72.23     | 79.38     | 8.09        | 5.0           |
|-------|-----------|-----------|-------------|---------------|
| minor | 63.18     | 67.79     | 6.48        | 5.0           |
|-------|-----------|-----------|-------------|---------------|
| rfy   | 48.89     | 52.07     | 1.41        | 0.076         |

## Torus Geometry Dimensions 

Dimensions of coil geometry, implemented as toruses (G4Torus).  G10 wall is actually a G4Torus 
with a hollow core, and a wall thickness as specified below    

| Label | rtor (cm) | rmin (cm) | rmax (cm) | length (cm) | G10 Wall (cm) |  
|-------|-----------|-----------|-----------|-------------|---------------|
| major | 75.81     | 0.0       | 7.15      | 8.09        | 5.0           |
|-------|-----------|-----------|-----------|-------------|---------------|
| minor | 65.48     | 0.0       | 4.61      | 6.48        | 5.0           |
|-------|-----------|-----------|-----------|-------------|---------------|
| rfy   | 50.48     | 0.0       | 3.18      | 1.41        | 0.076         |

