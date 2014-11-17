.SUFFIXES: 
.SUFFIXES: .f90 .o

FC = ifort
FFLAGS = -CB
INCLUDE = 
LIBS = 

EXE = WHAM.x

MODULES = bin.mod constant.mod precision_m.mod react_coord_bin.mod simulation.mod snapshot.mod wham.mod

OBJS = precision_m.o lib.o constant.o snapshot.o bin.o react_coord_bin.o simulation.o WHAM.o WHAM_caller.o

all:	${EXE}


$(EXE):$(OBJS) ${MODULES}
	$(FC) -o $@ $(FFLAGS) $(OBJS) $(LIBS)

%.o %.mod:%.f90
	$(FC) -c $(FFLAGS) $(INCLUDE) $<

include .depend

depend .depend:
	makedepf90 *.f90 > .depend

clean:
	/bin/rm -f $(EXE) $(OBJS) ${MODULES}

