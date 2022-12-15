#!/bin/bash

export DETECTOR=hadron_endcap
export DETECTOR_PATH=${EIC_SHELL_PREFIX}/share/${DETECTOR}
export DETECTOR_VERSION=master
export BEAMLINE_CONFIG=ip6
export BEAMLINE_CONFIG_VERSION=master

export JUGGLER_DETECTOR=$DETECTOR
export JUGGLER_DETECTOR_VERSION=$DETECTOR_VERSION
export JUGGLER_DETECTOR_PATH=$DETECTOR_PATH
export JUGGLER_BEAMLINE_CONFIG=$BEAMLINE_CONFIG
export JUGGLER_BEAMLINE_CONFIG_VERSION=$BEAMLINE_CONFIG_VERSION
export JUGGLER_INSTALL_PREFIX=${EIC_SHELL_PREFIX}

export LD_LIBRARY_PATH=${EIC_SHELL_PREFIX}/lib:$LD_LIBRARY_PATH
export PATH=${EIC_SHELL_PREFIX}/bin:$PATH

#HDF5                                                                                                                                                                    
export HDF5_DIR=${EIC_SHELL_PREFIX}/../../to_hdf5/hdf5-1.8.22/hdf5                                                                        
export PATH=$PATH:${HDF5_DIR}/bin                                                                                             
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${HDF5_DIR}/lib