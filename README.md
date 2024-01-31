# CARMA-1Dexoplanet
This repository is for running the baseline case of 1D cloud microphysical experiment with CARMA.

## File Description
-- tests/carma_nuc_grow.F90: This is the "main code" which defines the cloud formation process. It specifies temperature profile, moisture and aerosol transport and aerosol distributions, numerical details like time step and running length...   
-- run/carma: This is where output data is stored. Before running the model, an input file for T-P profile should be added into this folder. The file "tests/clima_mod.F90" works to read the input file. It should be edited if T-P profile changes.  
-- source/base: This folder is not originally created by the author of this repository. It is a standard version for CARMA, which is well documented in https://wiki.ucar.edu/display/CARMA/CARMA+Home.  
-- tests/clima_mod.F90: As described, this works to read the atmosphere profile input file. It should be updated if the input file is in a different form.   
-- source/base/carma_planet_mod.F90: Most of the planetary properties are specified in this file.  

## Running the model
1. Building:
```bash
./make-carma.csh clean
./make-carma.csh
```
2. Running
```bash
./run-all
```
To tune the cloud formation process, most of the edits can be done in source/base/carma_planet_mod.F90 and tests/carma_nuc_grow.F90.

## System requirements:
intel-mpi & ifort in gcc
