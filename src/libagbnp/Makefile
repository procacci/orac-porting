#******************************************************************************
# May 2004, Emilio Gallicchio
# Copyright 2004 Schrodinger, Inc
# Source level make file for AGBNP library.
#******************************************************************************
SRC_DIR = .
BUILD_DIR = .
SRC = $(SRC_DIR)/libagbnp.c $(SRC_DIR)/agbnp.h $(SRC_DIR)/agbnp_private.h 
OBJ = $(SRC_DIR)/libagbnp.o $(SRC_DIR)/libagbnp_pbc.o

CC = gcc 
LD = $(CC)
FC = $(F77)

###########################################################################

STATIC_LIBEXT = a
LIBEXT = $(STATIC_LIBEXT)
ARFLAGS = rs
CFLAGS = -c -O
LIBNBLIST_DIR = ../libnblist

all: libagbnp.$(LIBEXT)

clean: 
	rm -rf *.o agbnp_test

clean_all: clean
	rm -rf *.$(LIBEXT) 

#OpenMP compilation with gomp
#Linux-x86_64.ompgcc: LOCAL_CFLAGS = -fopenmp

libagbnp.$(STATIC_LIBEXT):  $(OBJ)
	$(AR) $(ARFLAGS) libagbnp.$(STATIC_LIBEXT) $(OBJ)

libagbnp.o: $(SRC_DIR)/libagbnp.c $(SRC_DIR)/agbnp.h agbnp_private.h
	$(CC) $(CFLAGS) -I$(LIBNBLIST_DIR) -DUSE_MMNBLIST $(SRC_DIR)/libagbnp.c

libagbnp_pbc.o: $(SRC_DIR)/libagbnp.c $(SRC_DIR)/agbnp.h agbnp_private.h
	$(CC) $(CFLAGS) -DUSE_PBC -DUSE_ALT -I$(LIBNBLIST_DIR) -DUSE_MMNBLIST -o libagbnp_pbc.o $(SRC_DIR)/libagbnp.c

agbnp_test: $(SRC_DIR)/agbnp_test.c $(SRC_DIR)/agbnp.h $(SRC_DIR)/libagbnp.$(LIBEXT) $(LIBNBLIST_DIR)/libnblist.h $(LIBNBLIST_DIR)/libnblist.$(LIBEXT)
	$(CC) -I$(LIBNBLIST_DIR) $(SRC_DIR)/agbnp_test.c -L$(SRC_DIR) -L$(LIBNBLIST_DIR) -lm -lagbnp -lnblist -o agbnp_test

test: agbnp_test agbnp.dat
	./agbnp_test < agbnp.dat
