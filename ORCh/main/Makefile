

all: Main

Main: 
	mpic++ -I$(ORCH_BASE)/Cantera/include -I$(TENSORFLOW_CAPI_PATH)/include mainDRGEPSpecies.cpp mainDRGEPReactions.cpp mainComputeTrajectories.cpp mainLumping.cpp -c -O3

clean:
	rm *.o
