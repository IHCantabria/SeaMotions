
# Change Log
All notable changes to this project will be documented in this file.
 

 
## [Unreleased] - yyyy-mm-dd
 
Up to date
 
### Added

 
### Changed
 
### Fixed

## [1.1.0] - 2024-08-26

Reviewed code to work with OWC
 
### Added

- New output channel to storage the mesh in the output file.
- Possibility to fix the body to earth coordinate system so
  there is no contribution of the hydromechanic forces for 
  this body
 
### Changed

- Serialized Froude-Krylov calculations inside QTF forces estimations
- Reformulated Froude-Krylov potential and its derivatives to not 
  use cosh formulation and to use exp formulation. It allows to have 
  a more robust behaviour when big function arguments are used
 
### Fixed

- Problems with domain decomposition across processors
- Problems related to the output file system

## [1.0.11] - 2024-07-26

First release version in Linux environment
 
### Added
 
### Changed
 
### Fixed

- Internal issues related to software compiling in linux environment
- Code clean up of unused variables
- Kochin functions module dependency function cos_alpha

## [1.0.0] - 2024-07-23

This is the first production version of the software. The most relevant features implemented are:  
- Calculation of hydromechanic coefficients (Freq Domain)  
- Calculation of first order wave exciting forces: diffraction and Froude-Krylov (Freq Domain)  
- Calculation of RAO operator (Freq Domain)  
- Calculation of second order wave exciting forces using two methods (Freq Domain):  
    - Pinkster  
    - Indirect  
- Calculation of wave relative height (Freq Domain)  
- Parallel execution using MPI (Freq Domain)  
- Output system based on HDF5 technology (Freq Domain).  
 
### Added
 
### Changed
 
### Fixed

## Format
The format is based on [Keep a Changelog](http://keepachangelog.com/)