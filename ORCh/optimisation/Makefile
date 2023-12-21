

all: Main

Main: 
	mpic++ -I$(ORCH_BASE)/Cantera/include OptimScenario.cpp optimisation.cpp CreateAnalyticDirectory.cpp -c -O3

clean:
	rm *.o
