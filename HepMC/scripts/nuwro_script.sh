#!/bin/bash

# Useage: 
# Copy this script somewhere onto /pnfs/uboone/resilient/ 
# Add to xml worklflow in the <initsource> field       

# Use PROCESS to select a unique hepmc file from list
echo Process num is ${PROCESS}

let FILENO="${PROCESS}+1"

export HEPMCFILE=$(sed -n ${FILENO}'p' < hepmc_files.list)
echo Hepmc file is ${HEPMCFILE}

# Copy the hepmc file onto the grid node
ifdh cp ${HEPMCFILE} nuwro_input.hepmc

#cat nuwro_input.hepmc 
#cat hepmc_POT.meta
