# Makefile for chrom-gc

CPP  = g++
CC   = g++

CXXINCLUDEDIR = 
CXXFLAGS = $(CXXINCLUDEDIR) -D_FILE_OFFSET_BITS=64 -Wall -ggdb -g3 -fno-inline-small-functions -O0 -fno-inline -fno-eliminate-unused-debug-types
RM = rm -f

OBJ  = chrom-gc.o \
	   Chromosome_dsbreak.o \
	   Chromosome_mutate.o \
	   Chromosome_repair0.o \
	   Chromosome_repair1.o

HEADER = Chromosome.h \
         GC.h \
         Histogram.h \
         RandBinomial.h \
         RandGeometric.h \
         RandUniform.h \
         RandUniform_GSL.h \
         SequenceRuns.h \
         VectorUtility.h

BIN  = chrom-gc

.PHONY: all all-before all-after clean clean-custom

all: all-before $(BIN) all-after


$(BIN): $(OBJ)
	$(CPP) $(CXXFLAGS) $(OBJ) -o $@ $(LIBS)

$(OBJ): $(HEADER)

clean: 
	$(RM) $(OBJ) $(BIN)

