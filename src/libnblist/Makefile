SRC_DIR = ./
SRC = $(SRC_DIR)/libnblist.c $(SRC_DIR)/libnblist.h

CC = gcc 
LD = $(CC)
FC = $(F77)

###########################################################################

STATIC_LIBEXT = a
LIBEXT =  $(STATIC_LIBEXT)
ARFLAGS = rs
CFLAGS = -c -O

all: libnblist.$(LIBEXT)

clean: 
	rm -rf *.o

clean_all: clean
	rm -rf *.$(LIBEXT) 


libnblist.$(STATIC_LIBEXT):  libnblist.$O
	$(AR) $(ARFLAGS) libnblist.$(STATIC_LIBEXT) libnblist.o

libnblist.$O: $(SRC_DIR)/libnblist.c $(SRC_DIR)/libnblist.h
	$(CC) $(CFLAGS) $(SRC_DIR)/libnblist.c
