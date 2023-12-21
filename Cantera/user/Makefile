#/bin/sh
###############################################################
#  $Author: hkmoffa $
#  $Date: 2008/01/01 19:19:50 $
#  $Revision: 1.7 $
#
#  Copyright 2002 California Institute of Technology
#
###############################################################

###############################################################
# list your object files here
# This is probably the only thing you need to change in this file.

OBJS    = user.o

###############################################################

INCDIR = $(GTCOMB_CT_HOME)/include

.SUFFIXES :
.SUFFIXES : .cpp .d .o

OBJDIR = .
PIC_FLAG=-fPIC
CXX_FLAGS = -O3 -Wall   $(CXX_OPT) $(PIC_FLAG)


CXX_INCLUDES = -I$(GTCOMB_CT_HOSTTYPE)/

LIB = $(GTCOMB_CT_HOSTTYPE)/libuser.a

DEPENDS = $(OBJS:.o=.d)

%.d: Makefile %.o
	g++ -MM $(CXX_INCLUDES) $*.cpp > $*.d

.cpp.o: user.h
	g++ -std=c++11 -w -c $< $(CXX_FLAGS) $(CXX_INCLUDES) -I$(BOOSTPATH)  -I../include/  -I../src/ 

.f.o:
	g77 -c $< $(F77_FLAGS) 

#all lib: $(LIB)
#	cp user.h ../../build/include/cantera/

$(LIB): $(OBJS)
	ar ruv $(LIB) $(OBJS) > /dev/null

clean:
	$(RM) *.o *~ $(LIB)

depends:
	make .depends

.depends: $(DEPENDS)
	cat *.d > .depends

$(OBJS): Makefile

install:

TAGS: 
	etags *.h *.cpp

ifeq ($(wildcard .depends), .depends)
include .depends
endif

