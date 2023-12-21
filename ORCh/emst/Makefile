

all: Main

Main: 
	mpifort -I$(CT_INC) emst_subs.f emst.f -c -O3

clean:
	rm *.o *.mod
