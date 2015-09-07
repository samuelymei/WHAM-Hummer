.SUFFIXES: 
.SUFFIXES: .f .f90 .o

FC = ifort
FFLAGS = -CB
INCLUDE = 
LIBS = 

EXE = WHAM.x

MODULES = control.mod bin.mod constant.mod precision_m.mod react_coord_bin.mod simulation.mod snapshot.mod wham.mod

OBJS = control.o precision_m.o lib.o constant.o snapshot.o bin.o react_coord_bin.o simulation.o WHAM.o WHAM_caller.o lnsrch.o dfpmin.o lbfgs.o

all:	${EXE}


$(EXE):$(OBJS) ${MODULES}
	$(FC) -o $@ $(FFLAGS) $(OBJS) $(LIBS)

lbfgs.o:lbfgs.f
	$(FC) -c $(FFLAGS) $(INCLUDE) $<

lnsrch.o:lnsrch.f
	$(FC) -c $(FFLAGS) $(INCLUDE) $<

dfpmin.o:dfpmin.f
	$(FC) -c $(FFLAGS) $(INCLUDE) $<

%.o %.mod:%.f90
	$(FC) -c $(FFLAGS) $(INCLUDE) $<

include .depend

depend .depend:
	makedepf90 *.f90 *.f > .depend

clean:
	/bin/rm -f $(EXE) $(OBJS) ${MODULES} lbfgs.o

