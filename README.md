# PET-codesharing

This directory contains MATLAB and C/MEX code for PET/CT simulation, image processing, and positron range modeling. The codebase supports research and development in PET imaging, including simulation of detector effects, positron range blurring, and attenuation correction.

## Directory Structure and Key Files

- **MATLAB Scripts**
  - `mainRichCodeshare.m`: Main script for running PET/CT simulation and image processing workflows.
  - `ct2mu.m`: Converts CT images to attenuation coefficients for PET attenuation correction.
  - `CTIm2RhoEstIm.m`: Converts CT images (HU) to estimated tissue density images.
  - `imPosRangeFilt_tripleImage2Dfaster.m`: Applies local positron range blurring to PET images using tissue-specific kernels.
  - `preComputeAllPosRangeKernelParams.m`: Precomputes kernel parameters and masks for positron range filtering.
  - `getPosRangeKernelSizeGa68.m`: Calculates kernel size for positron range blurring based on aPSF threshold.
  - `aPSFEstimate2D.m`: Computes a 2D anisotropic point spread function using exponential models.
  - `detResLoss2D.m`, `nonColl2D.m`: Simulate detector resolution loss and non-collinearity effects.
  - `checkDetectorRadii.m`: Checks detector positions relative to the PET ring.
  - `ParseChemicalFormula.m`, `PhysProps.m`, `PhotonAttenuationQ.m`: Utilities for chemical formula parsing and physical property lookup.

- **C/MEX Files**
  - `backProjLInv_BoxChunk.c`, `forwardProjLInv_BoxChunk.c`: High-performance MEX functions for forward and back projection in PET image reconstruction.
  - `getAttenuationMEX.c`: MEX function to compute attenuation coefficients for lines of response (LORs) in PET.

- **Data Files**
  - `PETCT_2D_Huber_25.mat`: Example PET/CT data for simulation.

## Getting Started

1. **Requirements:**
   - MATLAB (recommended: R2020a or newer)
   - C compiler for building MEX files (e.g., GCC, MSVC)

2. **Compiling MEX Files:**
   - In MATLAB, run:
     ```matlab
     mex backProjLInv_BoxChunk.c
     mex forwardProjLInv_BoxChunk.c
     mex getAttenuationMEX.c
     ```

3. **Running Simulations:**
   - Edit and run `mainRichCodeshare.m` to execute the PET/CT simulation pipeline.

## Notes
- Some scripts require specific data files (e.g., `PETCT_2D_Huber_25.mat`).
- Utility functions for physical properties and attenuation are adapted from NIST and other sources.

## License
This code is intended for academic and research use. Please cite appropriately if used in publications.