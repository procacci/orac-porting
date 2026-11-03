#################################################################
# *** $RCSfile: Makefile.pro,v $
# *** $Revision: 1.13 $
# *** $Date: 2009/08/24 21:19:23 $
#################################################################

####### USING THIS MAKEFILE: ########

# Makefile.pro contains code to make Makefiles for a variety of machines.

BUILD_ROOT = $(SCHRODINGER)
SOURCE_ROOT = $(SCHRODINGER_SRC)

DIRNAME = mmlibs/mmnblist
SRC_DIR = $(SOURCE_ROOT)/mmshare/$(DIRNAME)
BUILD_DIR = $(BUILD_ROOT)/$(MMSHARE)/$(DIRNAME)
CODE_DOC_DIR=$(BUILD_ROOT)/$(MMSHARE)/code_docs/$(DIRNAME)

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

OBJ =  mmnblist.$O mmnblistg.$O mmnblist_ilist.$O
SRC =  $(SRC_DIR)/mmnblist.c $(SRC_DIR)/mmnblist.h \
       $(SRC_DIR)/mmnblist_private.h \
       $(SRC_DIR)/mmnblistg.c $(SRC_DIR)/mmnblistg.h \
       $(SRC_DIR)/mmnblistg_private.h \
       $(SRC_DIR)/mmnblist_ilist.c $(SRC_DIR)/mmnblist_ilist.h \
       $(SRC_DIR)/libnblist.c $(SRC_DIR)/libnblist.h

MISC = Makefile.pro Makefile

LD = $(CC)
#Linux-ia64: LD = $(CXX)
#Linux-x86_64: LD = $(CXX)

#### Global macros

#### Local macros
FC = $(F77)

LOCAL_LIBS =  -lm -lmmerr -lmmlist -lmmbs -lm2io -lmmct -lmmhash -lmmat -lmmtoken -lmmfile -lmmio -lmmnblist -lnblist


###########################################################################

LIBEXT = $(STATIC_LIBEXT)

all: libnblist.$(LIBEXT) libmmnblist.$(LIBEXT)

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

install: $(DEV_LIB_DIR)/libnblist.$(LIBEXT) $(DEV_INC_DIR)/libnblist.h \
         $(DEV_LIB_DIR)/libmmnblist.$(LIBEXT) $(DEV_INC_DIR)/mmnblist.h \
         $(DEV_INC_DIR)/mmnblistg.h

includes: $(DEV_INC_DIR)/libnblist.h $(DEV_INC_DIR)/mmnblist.h $(DEV_INC_DIR)/mmnblistg.h

code_docs:
	@if test ! -d $(CODE_DOC_DIR) ; then \
		mkdir -p $(CODE_DOC_DIR) ; \
	fi
	perl $(BUILD_TOOLS)/createDoxygenConfig.pl \
	--docdir=$(CODE_DOC_DIR) \
	--name=MMnblist \
	--version=$(MMSHARE_VERSION)
	cd $(SRC_DIR) ; \
	doxygen $(BUILD_DIR)/doxygen.conf ; \
	cd $(BUILD_DIR)

$(DEV_LIB_DIR)/libnblist.$(LIBEXT): libnblist.$(LIBEXT)
	$(MY_CP) . $(DEV_LIB_DIR) libnblist.$(LIBEXT)

$(DEV_INC_DIR)/libnblist.h: $(SRC_DIR)/libnblist.h
	$(MY_CP) $(SRC_DIR) $(DEV_INC_DIR) libnblist.h

$(DEV_LIB_DIR)/libmmnblist.$(LIBEXT): libmmnblist.$(LIBEXT)
	$(MY_CP) . $(DEV_LIB_DIR) libmmnblist.$(LIBEXT)

$(DEV_INC_DIR)/mmnblist.h: $(SRC_DIR)/mmnblist.h
	$(MY_CP) $(SRC_DIR) $(DEV_INC_DIR) mmnblist.h

$(DEV_INC_DIR)/mmnblistg.h: $(SRC_DIR)/mmnblistg.h
	$(MY_CP) $(SRC_DIR) $(DEV_INC_DIR) mmnblistg.h

# test - no diff means a good test:
test: mmnblist_test$(EXE) prot.dat
	$(SCHRODINGER)/run ./mmnblist_test > mmnblist_test.out
	$(DIFF) mmnblist_test.out $(SRC_DIR)/mmnblist_test.out.std $(COMPARE)

prot.dat: $(SRC_DIR)/prot.dat
	$(MY_CP) $(SRC_DIR) . prot.dat

mmnblist_test$(EXE):  mmnblist_test.$(O) $(DEV_LIB_DIR)/libmmnblist.$(LIBEXT) $(DEV_LIB_DIR)/libnblist.$(LIBEXT)
	$(LD) $(LFLAGS) mmnblist_test.$(O) $(LIBS) $(MACH_SYSLIBS) \
	$(EXE_OUT_FLAG) mmnblist_test


mmnblist_test.$(O): $(SRC_DIR)/mmnblist_test.c $(DEV_INC_DIR)/mmnblist.h $(DEV_INC_DIR)/mmnblistg.h
	$(CC) $(CFLAGS) $(SRC_DIR)/mmnblist_test.c

clean: 
	rm -rf *.$O mmnblist_test.out mmnblist_test

clean_all: clean
	rm -rf *.$(LIBEXT) *.h prot.dat

clean_installed: clean_all
	rm -rf $(DEV_LIB_DIR)/libmmnblist.$(LIBEXT) $(DEV_INC_DIR)/mmnblist.h $(DEV_INC_DIR)/mmnblistg.h  $(DEV_LIB_DIR)/libnblist.$(LIBEXT) $(DEV_INC_DIR)/libnblist.h

libnblist.$(STATIC_LIBEXT):  libnblist.$O
	$(AR) $(ARFLAGS) $(STATIC_OUT_FLAG)libnblist.$(STATIC_LIBEXT) libnblist.$O

libnblist.$(SHARED_LIBEXT):  libnblist.$O
	$(SHARED_LINKER) -o libnblist.$(LIBEXT) $(SHARED_LFLAGS) libnblist.$O

libmmnblist.$(STATIC_LIBEXT): 	$(OBJ) $(DEV_LIB_DIR)/libnblist.$(LIBEXT) $(DEV_INC_DIR)/libnblist.h
	$(AR) $(ARFLAGS) $(STATIC_OUT_FLAG)libmmnblist.$(STATIC_LIBEXT) $(OBJ)

libmmnblist.$(SHARED_LIBEXT): 	$(OBJ) $(DEV_LIB_DIR)/libnblist.$(LIBEXT) $(DEV_INC_DIR)/libnblist.h
	$(SHARED_LINKER) -o libmmnblist.$(LIBEXT) $(SHARED_LFLAGS) $(OBJ)

libnblist.$O: $(SRC_DIR)/libnblist.c $(DEV_INC_DIR)/libnblist.h
	$(CC) $(CFLAGS) $(SRC_DIR)/libnblist.c

mmnblist.$O: $(SRC_DIR)/mmnblist.c $(DEV_INC_DIR)/mmnblist.h $(DEV_INC_DIR)/libnblist.h mmnblist_private.h
	$(CC) $(CFLAGS) $(SRC_DIR)/mmnblist.c

mmnblistg.$O: $(SRC_DIR)/mmnblistg.c $(DEV_INC_DIR)/mmnblistg.h $(DEV_INC_DIR)/mmnblist.h  $(DEV_INC_DIR)/libnblist.h mmnblistg_private.h mmnblist_private.h mmnblist_ilist.h 
	$(CC) $(CFLAGS) $(SRC_DIR)/mmnblistg.c

mmnblist_ilist.$O: $(SRC_DIR)/mmnblist_ilist.c  mmnblist_ilist.h
	$(CC) $(CFLAGS) $(SRC_DIR)/mmnblist_ilist.c

mmnblist_private.h: $(SRC_DIR)/mmnblist_private.h
	$(MY_CP) $(SRC_DIR) . mmnblist_private.h

mmnblistg_private.h: $(SRC_DIR)/mmnblistg_private.h
	$(MY_CP) $(SRC_DIR) . mmnblistg_private.h

mmnblist_ilist.h: $(SRC_DIR)/mmnblist_ilist.h
	$(MY_CP) $(SRC_DIR) . mmnblist_ilist.h


