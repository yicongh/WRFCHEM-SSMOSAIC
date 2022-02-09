# ==========================================================================
#           THIS IS THE MAKEFILE FOR THE SIMPLESOM-MOSAIC MODULE
# ==========================================================================

# OBJECT FILES:
OBJs = main.ss_mosaic.o param.ss_mosaic.o sub.step_soacond.o \
       sub.map_inout.o sub.step_olig.o sub.step_hetchem.o

# DEFAULT TARGET:
# ==========================================================================
default: $(OBJs)
	@cd ./objs/; ar -ru ../../../main/libwrflib.a $(OBJs)
	@cp ./objs/*.mod ../../chem/

# RECIPE PATTERN:
# ==========================================================================
%.o : %.f90 makefile
	@ifort -c -I ../../frame/ -I ./objs/ -o ./objs/$@ -module ./objs/ $<

# DEPENDENCIES:
# ==========================================================================
main.ss_mosaic.o: param.ss_mosaic.o sub.map_inout.o sub.step_soacond.o \
	          sub.step_olig.o sub.step_hetchem.o

param.ss_mosaic.o:

sub.step_soacond.o: param.ss_mosaic.o

sub.step_olig.o: param.ss_mosaic.o

sub.step_hetchem.o: param.ss_mosaic.o

sub.map_inout.o: param.ss_mosaic.o

# UTILITIES:
# ==========================================================================
clean:
	@rm ./objs/*
