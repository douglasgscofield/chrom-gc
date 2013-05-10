# Makefile for chrom-gc

CPP  = g++
CC   = g++

CXXINCLUDEDIR = 
CXXLDFLAGS = -g3
CXXFLAGS = $(CXXINCLUDEDIR) -Wall -g3 -fno-inline 
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
	$(CPP) $(CXXLDFLAGS) $(OBJ) -o $@ $(LIBS)

$(OBJ): $(HEADER)

clean: 
	$(RM) $(OBJ) $(BIN)

