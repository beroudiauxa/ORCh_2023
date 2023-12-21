ORCh=../.

all: conditions mainProgram

conditions : 
	gcc -c conditions.cpp

mainProgram: ../main/mainProgram.cpp
	mpic++ -pthread -L$(SUNDIALS_LIB) -I$(SUNDIALS_INC) -L$(MPI_LIB) -I$(BOOSTPATH) -I$(ORCH_BASE)/Cantera/include -L$(CT_LIB) -L$(GTCOMB_CT_HOSTTYPE) -L$(TENSORFLOW_CAPI_PATH)/lib -I$(TENSORFLOW_CAPI_PATH)/include -L$(OPENCV_PATH)/lib -I$(OPENCV_PATH)/include -o mainProgram $(ORCh)/main/mainProgram.cpp  conditions.o  $(ORCh)/tools/tools.o $(ORCh)/tools/gnuplot.o $(ORCh)/tools/fitness_criteria.o $(ORCh)/lumping/lumping.o $(ORCh)/lumping/lumping_species.o $(ORCh)/cantera/computeMultipleInlet.o $(ORCh)/cantera/computePremixedFlames.o $(ORCh)/cantera/computeAutoIgnition.o $(ORCh)/cantera/Analytic_function.o $(ORCh)/cantera/particle.o $(ORCh)/cantera/flamemodel.o $(ORCh)/drgep/drgep.o $(ORCh)/read_write/species.o $(ORCh)/read_write/reaction.o $(ORCh)/read_write/read.o $(ORCh)/read_write/write.o $(ORCh)/read_write/write_QSS.o $(ORCh)/read_write/write_QSS_FORTRAN.o $(ORCh)/read_write/QSSscenario.o $(ORCh)/optimisation/OptimScenario.o $(ORCh)/optimisation/optimisation.o $(ORCh)/optimisation/CreateAnalyticDirectory.o $(ORCh)/main/mainDRGEPSpecies.o $(ORCh)/main/mainDRGEPReactions.o $(ORCh)/main/mainComputeTrajectories.o $(ORCh)/main/mainLumping.o $(ORCh)/tools/outputs.o $(ORCh)/emst/emst_subs.o $(ORCh)/emst/emst.o $(ORCH)/cppflow/src/Model.o $(ORCH)/cppflow/src/Tensor.o -O3 -lm -lhdf5 -lhdf5_hl -lcantera -luser -lifcore -ltensorflow -ltensorflow_framework -lopencv_core -lsundials_cvodes -lsundials_nvecserial -lgfortran

clean:
	rm -f mainProgram *.o
