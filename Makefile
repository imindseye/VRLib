####################################################################

SHELL = /bin/sh

#######################################

LIBNAME = volren
RM = rm -f 
AR = ar cq

C++ = g++
CCFLAGS = -g  -DREAL=float 

TOP = ..

INCLUDE = -I. 

OBJS = Map.o Trans_Stack.o render.o image.o  render_aux.o image_composite.o
  
SRCS = Map.C Trans_Stack.C render.C image.C  render_aux.C image_composite.C

.SUFFIXES: .C
.C.o:
	$(C++) $(CCFLAGS) $(INCLUDE) -c $<

default: all

all: lib$(LIBNAME).a  testmain 

lib$(LIBNAME).a : $(OBJS) render.h
	$(RM) $@
	$(AR) $@ $(OBJS)

## a simple test program 
testmain: testmain.o lib$(LIBNAME).a 
	$(C++) -o testmain testmain.o -L. -l$(LIBNAME) -lm 

###########################################################

clean:
	rm -f *.o *.a

depend: Makefile $(SRCS)
	makedepend -fMakefile $(CCFLAGS) $(INCLUDE) $(SRCS)

# DO NOT DELETE THIS LINE -- make depend depends on it.
