#******************************************************************************
# May 2004, Emilio Gallicchio
# Copyright 2004 Schrodinger, Inc
# Source level make file for AGBNP library.
#******************************************************************************

BUILD_ROOT = $(SCHRODINGER)
SOURCE_ROOT = $(SCHRODINGER_SRC)

DIRNAME = src/libagbnp
SRC_DIR = $(SCHRODINGER_SRC)/impact-src/$(DIRNAME)
BUILD_DIR = $(SCHRODINGER)/$(IMPACT)/$(DIRNAME)

# build tools
BUILD_TOOLS = $(SOURCE_ROOT)/mmshare/build_tools
PREPMAKE = $(BUILD_TOOLS)/prepmake2 $(BUILD_TOOLS)

# To create the Makefile:
#	cp Makefile.pro Makefile # (only if there's no Makefile present)
#	make pmakefile VERSION=SGI-g # (SGI-g comes from the following list)
#	make version		 # (sanity check; should say "SGI-g")
#	make clean		 # (in case .o files were compiled otherwise)
#	make install		 # (makes and installs libraries)

# sets of files:

OBJ =  libagbnp.$(O) libagbnp_pbc.$(O)
SRC =  $(SRC_DIR)/libagbnp.c $(SRC_DIR)/agbnp.h $(SRC_DIR)/agbnp_private.h 

MISC = Makefile.pro Makefile

#### Local macros
FC = $(F77)
LOCAL_LIBS =  -lm -lmmnblist -lnblist
IMPACT_INC_DIR = $(SCHRODINGER)/$(IMPACT)/include
IMPACT_LIB_DIR = $(SCHRODINGER)/$(IMPACT)/lib/$(OS_CPU)
#### Global macros

#Linux.omp: CC = icc
#Linux.omp: LOCAL_CFLAGS = -openmp

#OpenMP compilation with gomp
#Linux-x86_64.ompgcc: LOCAL_CFLAGS = -fopenmp

###########################################################################
LIBEXT = $(STATIC_LIBEXT)

all: build_tree install

pmakefile: env_vars
	$(PREPMAKE) $(SRC_DIR) $(VERSION)
	mv -f Makefile.tmp Makefile

version:
	@$(MAKE_VERSION) VERSION=$(VERSION)

env_vars:
	@if test ! "$(SCHRODINGER)" ; then \
        echo `pwd`/Makefile: 'Error: SCHRODINGER is undefined!'; \
        exit 1; \
    fi
	@if test ! "$(SCHRODINGER_SRC)" ; then \
        echo `pwd`/Makefile: 'Error: SCHRODINGER_SRC is undefined!'; \
        exit 1; \
    fi
	@if test ! "$(VERSION)" ; then \
        echo `pwd`/Makefile: 'Error: VERSION is undefined!'; \
        exit 1; \
    fi

globals: env_vars
	@if test ! "$(MY_CP)" ; then \
        echo `pwd`/Makefile: 'Error: MY_CP is undefined!'; \
        echo `pwd`/Makefile: 'Run "make pmakefile" to define global macros'; \
        exit 1; \
    fi

build_tree: globals
	$(DO_NOTHING)

install: $(IMPACT_LIB_DIR)/libagbnp.$(LIBEXT) $(IMPACT_INC_DIR)/agbnp.h

$(IMPACT_LIB_DIR)/libagbnp.$(LIBEXT): libagbnp.$(LIBEXT)
	$(MY_CP) . $(IMPACT_LIB_DIR) libagbnp.$(LIBEXT)

$(IMPACT_INC_DIR)/agbnp.h: $(SRC_DIR)/agbnp.h
	$(MY_CP) $(SRC_DIR) $(IMPACT_INC_DIR) agbnp.h

agbnp.tgz: $(IMPACT_INC_DIR)/agbnp.h  $(IMPACT_LIB_DIR)/libagbnp.$(LIBEXT) $(DEV_INC_DIR)/libnblist.h  $(DEV_LIB_DIR)/libnblist.$(LIBEXT) $(SRC_DIR)/agbnp.param $(SRC_DIR)/README
	rm -rf agbnp ; mkdir agbnp ; \
        cp  $(IMPACT_INC_DIR)/agbnp.h  $(IMPACT_LIB_DIR)/libagbnp.$(LIBEXT) \
        $(DEV_INC_DIR)/libnblist.h  $(DEV_LIB_DIR)/libnblist.$(LIBEXT) \
        $(SRC_DIR)/agbnp.param $(SRC_DIR)/README agbnp/ ; \
        tar zcvf agbnp.tgz agbnp ; rm -rf agbnp/


clean: 
	rm -rf *.$(O)

clean_all: clean
	rm -rf *.$(LIBEXT) *.h

clean_installed: clean_all
	rm -rf $(IMPACT_INC_DIR)/agbnp.h $(IMPACT_LIB_DIR)/libagbnp.$(LIBEXT)

libagbnp.$(STATIC_LIBEXT):  $(OBJ)
	$(AR) $(ARFLAGS) $(STATIC_OUT_FLAG)libagbnp.$(STATIC_LIBEXT) $(OBJ)

#libagbnp.$(SHARED_LIBEXT):  libagbnp.$(O)
#	$(SHARED_LINKER) -o libagbnp.$(SHARED_LIBEXT) $(SHARED_LFLAGS) $(OBJ)

libagbnp.$(O): $(SRC_DIR)/libagbnp.c $(IMPACT_INC_DIR)/agbnp.h agbnp_private.h
	$(CC) $(CFLAGS) -DUSE_MMNBLIST -I$(IMPACT_INC_DIR) -I$(DEV_INC_DIR) $(OBJ_OUT_FLAG)libagbnp.$(O) $(SRC_DIR)/libagbnp.c


libagbnp_pbc.$(O): $(SRC_DIR)/libagbnp.c $(IMPACT_INC_DIR)/agbnp.h agbnp_private.h
	$(CC) $(CFLAGS) -DUSE_PBC -DUSE_ALT -DUSE_MMNBLIST -I$(IMPACT_INC_DIR) -I$(DEV_INC_DIR) $(OBJ_OUT_FLAG)libagbnp_pbc.$(O) $(SRC_DIR)/libagbnp.c

#Settings for pgcc with OpenMP:
#
#Compilation:
# pgcc -mp -fast -Minline=name:ogauss_incremental
# pgcc -mp -fast -Minline=name:ogauss_incremental_pbc
#
#add to VERS_SYSLIBS:
# $(PGI)/linux86/lib/libpgmp.a $(PGI)/linux86/lib/libpgthread.a
#
#Linux.ifort.openmp: F77 = ifort -openmp
#Linux.ifort.openmp: F90 = ifort -openmp
agbnp_private.h: $(SRC_DIR)/agbnp_private.h
	$(MY_CP) $(SRC_DIR) . agbnp_private.h


agbnp_test: $(SRC_DIR)/agbnp_test.c $(IMPACT_INC_DIR)/agbnp.h $(IMPACT_LIB_DIR)/libagbnp.$(LIBEXT) $(DEV_INC_DIR)/libnblist.h $(DEV_LIB_DIR)/libnblist.$(LIBEXT)
	$(CC) -I$(IMPACT_INC_DIR) -I$(DEV_INC_DIR) $(SRC_DIR)/agbnp_test.c -L$(IMPACT_LIB_DIR) -L$(DEV_LIB_DIR) -lm -lagbnp -lnblist -o agbnp_test

agbnp.dat: $(SRC_DIR)/agbnp.dat
	$(MY_CP) $(SRC_DIR) $(BUILD_DIR) agbnp.dat

test: agbnp_test agbnp.dat
	./agbnp_test < agbnp.dat
