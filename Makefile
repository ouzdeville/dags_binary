# Project: shrivastava
# Makefile created by ouzdeville
CC       = gcc
OBJ      = gf.o matrix.o rng.o key_gen.o util.o kem.o PQCgenKAT_kem.o  encapsulation.o poly.o decoding.o decapsulation.o 
LINKOBJ  = gf.o matrix.o rng.o key_gen.o util.o kem.o PQCgenKAT_kem.o  encapsulation.o poly.o decoding.o decapsulation.o 
LIBS     = -L/usr/lib -L. -lssl -lcrypto -lkeccak
INCS     = 
CXXINCS  = 
BIN      = PQCgenKAT_kem
LFLAGS=
CFLAGS= -c -Wall -I. 
RM       = rm -f

.PHONY: all all-before all-after clean clean-custom

all: all-before $(BIN) all-after

clean: clean-custom
	${RM} $(OBJ) $(BIN)

$(BIN): $(OBJ)
	$(CC) -O3 $(LINKOBJ) -o $(BIN) $(LIBS)

%.o: %.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -c -o $@ $<
	
gen: clean-custom
	$(CC) main_genparams.c -o main_genparams