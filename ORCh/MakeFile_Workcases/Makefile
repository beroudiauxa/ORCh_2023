
ORCh=../../ORCh




all: conditions mainProgram

conditions : 
	gcc -c conditions.cpp

mainProgram: ../../ORCh/main/mainProgram.cpp
	mpic++ -w -std=c++11 -L$(MPI_LIB) -I$(BOOSTPATH) -L$(BOOSTPATH) -I$(GTCOMB_CT_HOME)/include -L$(GTCOMB_CT_HOME)/build/lib -L$(GTCOMB_CT_HOSTTYPE) -o mainProgram $(ORCh)/main/mainProgram.cpp  conditions.o  $(ORCh)/tools/tools.o $(ORCh)/tools/gnuplot.o $(ORCh)/tools/fitness_criteria.o $(ORCh)/lumping/lumping.o $(ORCh)/lumping/lumping_species.o $(ORCh)/cantera/computeMultipleInlet.o $(ORCh)/cantera/computePremixedFlames.o $(ORCh)/cantera/computeAutoIgnition.o $(ORCh)/cantera/Analytic_function.o $(ORCh)/cantera/particle.o $(ORCh)/cantera/flamemodel.o $(ORCh)/drgep/drgep.o $(ORCh)/read_write/species.o $(ORCh)/read_write/reaction.o $(ORCh)/read_write/read.o $(ORCh)/read_write/write.o $(ORCh)/read_write/write_QSS.o $(ORCh)/read_write/write_QSS_FORTRAN.o $(ORCh)/read_write/QSSscenario.o $(ORCh)/optimisation/OptimScenario.o $(ORCh)/optimisation/optimisation.o $(ORCh)/optimisation/CreateAnalyticDirectory.o $(ORCh)/main/mainDRGEPSpecies.o $(ORCh)/main/mainDRGEPReactions.o $(ORCh)/main/mainComputeTrajectories.o $(ORCh)/main/mainLumping.o $(ORCh)/tools/outputs.o $(ORCh)/emst/emst_subs.o $(ORCh)/emst/emst.o -O3 -Wall -lm -lcantera -luser -lifcore



clean:
	rm -f mainProgram *.o
