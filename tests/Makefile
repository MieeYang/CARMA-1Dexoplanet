# Sub makefile for test files

# This is intended to be included by another makefile to actually
# build the system. It has all the dependency information for the
# files in the test tree.

# NOTE: In the future (or in parnet models) this could be replaced
# by automatic dependency generation.

TEST_OBJ = clima_mod.o carma_nuc_grow.o

TEST_DOC = clima_mod.html carma_nuc_grow.html

clima_mod.o : clima_mod.F90
	$(FORTRAN) $(FFLAGS) -c $<

carma_nuc_grow.o : carma_nuc_grow.F90 clima_mod.mod carma_mod.mod carma_constants_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(FFLAGS) -c $<

