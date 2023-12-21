

all: Main

Main: 
	mpic++ -I$(ORCH_BASE)/Cantera/include tools.cpp gnuplot.cpp  outputs.cpp fitness_criteria.cpp -c -O3

clean:
	rm *.o
