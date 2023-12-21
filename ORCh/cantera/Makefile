

all: Main

Main: 
	mpic++ -I$(ORCH_BASE)/Cantera/include -I$(TENSORFLOW_CAPI_PATH)/include -I$(OPENCV_PATH)/include computeMultipleInlet.cpp computePremixedFlames.cpp computeAutoIgnition.cpp Analytic_function.cpp particle.cpp flamemodel.cpp -c -O3
clean:
	rm *.o
