# Makefile to build lbstime.a library and driver routines with PAPI 
#
# Andrew J. Pounds, Ph.D.
# Departments of Chemistry and Computer Science
# Mercer University
# Spring 2018
#

# Use gfortran to avoid g77/gcc incompatibilities on Fedora
# g77 works fine on systems where gcc/g77 were linked against the same library
CC = gcc 
CFLAGS = -mtune=nocona -march=nocona
PROF =  

LIBS =  liblbstime.a -lpapi

OBJS = cputime.o walltime.o  

all: lib cdriver

cputime.o : cputime.cc   
	$(CC) $(CFLAGS) -c cputime.cc  

walltime.o : walltime.cc   
	$(CC) $(CFLAGS) -c walltime.cc  

cdriver.o : cdriver.c   
	$(CC) $(CFLAGS) -c cdriver.c   

# Don't forget the -lstdc++
cdriver : cdriver.o  $(OBJS) 
	$(CC) $(CFLAGS) -o cdriver cdriver.o  $(LIBS) -lstdc++   

# Default Targets for Cleaning up the Environment
clean :
	rm *.o
	rm *.a

pristine :
	rm *.o
	rm *.a
	touch *.cc *.c *.f 
	echo ${SHELL}
	`if [ -e "cdriver" ]  ; then rm cdriver; fi`

ctags :
	ctags *.cc *.c *.f

# Target for making the library

lib: $(OBJS) 
	ar -rc liblbstime.a $(OBJS) 
	ranlib liblbstime.a
