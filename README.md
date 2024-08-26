# SeaMotions

# Introduction

SeaMotions is a software designed to simulate the interaction in between the waves and floating bodies. It is based on the Boundary Element Method and employs special Green functions to accelerate the computation of the sources on the frequency domain.

## Features

### First order module
- Calculation of AddedMass/Damping coefficients (Hydromechanic Forces)
- Calculation of wave diffraction and Froude-Krylov forces (Wave Exciting Forces)
- Multi-Body formulation which allows to simulate OWCs
- Calculation of wave height along the domain
- Calculation of relative wave height
- Calculation of velocities and pressures along the domain

### Second order module
- Calculation of Mean Drift (Second Order force estimated through first order contributions)
- QTF coefficients
- Second order potential
    - Pinkster
    - Indirect method (X.Chen)

The calculation of the second order quantities is performed in serial.

## Contact
For any software related issues contact with the software developer: [Sergio Fernandez Ruano](https://ihcantabria.com/directorio-personal/sergio-fernandez-ruano/)