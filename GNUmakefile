# $Id: GNUmakefile,v 1.2 2000-10-19 12:22:10 stanaka Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := exampleB4a
G4TARGET := $(name)
G4EXLIB := true


ifndef G4INSTALL
  G4INSTALL =/opt/geant4.10.1-install/share/Geant4-10.1.1/geant4make
endif

.PHONY: all
all: lib bin

# ROOT support
#CPPFLAGS += $(shell $(ROOTSYS)/bin/root-config --cflags)
#CPPFLAGS += -I$(ROOTSYS)/include
#EXTRALIBS += $(shell $(ROOTSYS)/bin/root-config --libs) 
#CPPFLAGS += 'root-config --cflags'
#EXTRALIBS += 'root-config --nonew --libs'

include $(G4INSTALL)/config/binmake.gmk

visclean:
	rm -f g4*.prim g4*.eps g4*.wrl
	rm -f .DAWN_*

