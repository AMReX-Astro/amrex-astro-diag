PRECISION = DOUBLE
PROFILE = FALSE

DEBUG = FALSE

DIM = 2

COMP = g++

BL_NO_FORT = TRUE

USE_MPI = FALSE
USE_OMP = FALSE

USE_REACT = TRUE
USE_CXX_EOS = TRUE

MAX_ZONES := 16384

DEFINES += -DNPTS_MODEL=$(MAX_ZONES)

# programs to be compiled
EBASE := sprad

# EOS and network
EOS_DIR := helmholtz

# this should match what we used for the simulation
NETWORK_DIR := nova

# this should be RKC so we have access to the circle theorem function
INTEGRATOR_DIR := RKC

Bpack := ./Make.package
Blocs := . ..

EXTERN_SEARCH += . ..

include $(MICROPHYSICS_HOME)/Make.Microphysics
