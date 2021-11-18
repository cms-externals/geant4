# ----------------------------------------------------------------
# Makes test program in environment variable G4TARGET.
# ----------------------------------------------------------------

ifndef G4TARGET
  G4TARGET := $(TESTTARGET)
endif

ifndef G4INSTALL
  G4INSTALL = ../../../..
endif

G4EXEC_BUILD = true

include $(G4INSTALL)/config/architecture.gmk

# Override some variables for binmake.gmk.
#
INCFLAGS := -I$(G4BASE)/geometry/management/include \
            -I$(G4BASE)/geometry/volumes/include \
            -I$(G4BASE)/geometry/navigation/include \
            -I$(G4BASE)/geometry/magneticfield/include \
            -I$(G4BASE)/geometry/solids/CSG/include \
            -I$(G4BASE)/geometry/solids/Boolean/include \
            -I$(G4BASE)/materials/include \
            -I$(G4BASE)/global/management/include \
            -I$(G4BASE)/global/HEPRandom/include \
            -I$(G4BASE)/global/HEPGeometry/include \
            -I$(G4BASE)/graphics_reps/include \
            -I$(G4BASE)/intercoms/include

GLOBLIBS := $(shell test -f $(G4LIB)/$(G4SYSTEM)/libG4global.a && echo yes)

ifdef GLOBLIBS
  LDLIBS   := -lG4geometry -lG4graphics_reps -lG4materials \
              -lG4intercoms -lG4global
else
  LDLIBS   := -lG4csg -lG4geomBoolean -lG4navigation -lG4volumes \
              -lG4magneticfield -lG4geometrymng \
              -lG4materials -lG4graphics_reps \
              -lG4intercoms -lG4UIbasic -lG4UIcommon \
              -lG4hepnumerics -lG4heprandom -lG4globman 
endif

#  For external CLHEP must not include G4clhep
#ifndef CLHEP_BASE_DIR
LDLIBS   += -lG4clhep
# endif            

include $(G4INSTALL)/config/binmake.gmk

.PHONY: clean_all

clean_all:
	@for i in *.cc; do \
	(target=`basename $$i .cc`; \
	echo  "Removing $$target ... "; \
	$(MAKE) clean G4SKIP_DEPEND=1 G4TARGET=$$target); done
	@echo  "Wiping output file ... "
	@$(RM) -rf test.out
