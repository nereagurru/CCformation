This readme file was generated on 2025-07-06 by Nerea Gurrutxaga

## GENERAL INFORMATION

Title of code: mcdust for CC formation
Description: this code calculates the dust evolution of particles 

Author Information
Name: Nerea Gurrutxaga
ORCID: 0009-0008-3256-9564
Institution: Max Planck Institute for Solar System Research
Email: gurrutxaga@mps.mpg.de


## SHARING/ACCESS Code

The code will be made public ones the paper is public

## FILE OVERVIEW

Main files:
- src: all .F90 programs
- disk: data from diskevol (download data inside disk file to run this code)
- setup: compilation and parameter information
- outputs: for data storaged
- scripts: generate all figures from the main manuscript
- obj: storage of temporary ".mod" and ".o" files

## Prerequisites

`gfortran` `hdf5-serial` `python`

To install the required software in Ubuntu (this requires root permissions):
`sudo apt-get install gfortran`
`sudo apt-get install libhdf5-serial-dev`

Python is not required to run the code. But if you want to use the routines to read/write data from the simulation you will need a python installation.

## To compile the code: 

in the /setups/global/ file: `make`

This uses the setup files from the default run in `/setups/global/` 
 An executable in the name of the setup file will be created in the root directory which is `global` in this case.

parameter file: setups/global/setup.par


## To run the code: 

in the /setups/global/ file: `export OMP_NUM_THREADS=number_of_CPUs`

in the /setups/global/ file: `./global setup.par`

To clear the run

in the /setups/global/ file: `make clean`

Once the global output is generated, run local simulation:

in the /setups/local/ file: `make`
in the /setups/local/ file: `./global setup.par`
in the /setups/local/ file: `make clean`



