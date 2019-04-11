CFLAGS=-ansi -Wall -pedantic -g
g=g++

MDP= genericMDP discountedMDP  averageMDP totalRewardMDP
SOLUTION = solutionMDP feedbackSolutionMDP nonStationarySolutionMDP
UTILITES = alglin
INCLUDEDIR = /header
MARMOTECORE = binarySequence binarySimplex discountedMDP_LPsolve eventMixture marmoteBox marmoteInterval marmoteSet multiDimHomTransition sparseMatrix transitionStructure
WRAPPER = totalRewardMDP_wrap sparseMatrix_wrap marmoteInterval_wrap  

MARMOTECOREO =\
	$(addsuffix .o, $(MARMOTECORE) ) 

OBJECTS=\
	$(addsuffix .o, $(MDP) ) \
	$(addsuffix .o, $(SOLUTION) )  \
	$(addsuffix .o, $(UTILITES) ) \



POINTH=\
	$(addprefix header/,$(addsuffix .h, $(MDP) )) \
	$(addprefix header/,$(addsuffix .h, $(SOLUTION) ))  \
	$(addprefix header/,$(addsuffix .h, $(UTILITES) )) \
	$(addprefix header/,$(addsuffix .h, $(MARMOTECORE) ))


all: $(MDP) $(SOLUTION) $(UTILITES) $(MARMOTECORE) librairies1 librairies $(WRAPPER)

genericMDP: genericMDP.cpp
	$g $(CFLAGS) -c $<  -I$(INCLUDEDIR)

discountedMDP: discountedMDP.cpp 
	$g $(CFLAGS) -c $<  -I$(INCLUDEDIR)

averageMDP: averageMDP.cpp
	$g $(CFLAGS) -c $<  -I$(INCLUDEDIR)

totalRewardMDP: totalRewardMDP.cpp
	$g -fpic $(CFLAGS) -c $<  -I$(INCLUDEDIR)

solutionMDP: solutionMDP.cpp
	$g $(CFLAGS) -c $<  -I$(INCLUDEDIR)

feedbackSolutionMDP: feedbackSolutionMDP.cpp
	$g $(CFLAGS) -c $< -I$(INCLUDEDIR)

nonStationarySolutionMDP: nonStationarySolutionMDP.cpp
	$g $(CFLAGS) -c $<  -I$(INCLUDEDIR)

alglin : alglin.c
	$g $(CFLAGS) -c $< -I$(INCLUDEDIR) 

binarySequence: binarySequence.cpp
	$g $(CFLAGS) -c $<  -I$(INCLUDEDIR)

binarySimplex: binarySimplex.cpp
	$g $(CFLAGS) -c $<  -I$(INCLUDEDIR)

discountedMDP_LPsolve: discountedMDP_LPsolve.cpp
	$g $(CFLAGS) -c $<  -I$(INCLUDEDIR)	

eventMixture: eventMixture.cpp 
	$g $(CFLAGS) -c $<  -I$(INCLUDEDIR)

marmoteBox: marmoteBox.cpp
	$g $(CFLAGS) -c $<  -I$(INCLUDEDIR)

marmoteInterval: marmoteInterval.cpp
	$g -fpic $(CFLAGS) -c $<  -I$(INCLUDEDIR)

marmoteSet: marmoteSet.cpp
	$g $(CFLAGS) -c $<  -I$(INCLUDEDIR)

multiDimHomTransition: multiDimHomTransition.cpp
	$g $(CFLAGS) -c $<  -I$(INCLUDEDIR)

sparseMatrix: sparseMatrix.cpp
	$g -fpic $(CFLAGS) -c $<  -I$(INCLUDEDIR)

transitionStructure: transitionStructure.cpp
	$g $(CFLAGS) -c $<  -I$(INCLUDEDIR)

clean:
	/bin/rm -v *.o *.so *.py 

librairies1: $(MARMOTECOREO)
	ar crv libMarmoteCORE.a $(MARMOTECOREO)

librairies: $(OBJECTS) 
	ar crv libMarmoteMDP.a $(OBJECTS)

totalRewardMDP_wrap : totalRewardMDP_wrap1
	$g -shared  totalRewardMDP_wrap.o totalRewardMDP.o -L/libMarmoteMDP.a  -L/libMarmoteCore.a   -o _totalReward.so -lstdc++

totalRewardMDP_wrap1 :totalRewardMDP_wrap2 
	$g -fpic -c  totalRewardMDP_wrap.cxx header/totalRewardMDP.h -L/libMarmoteMDP.a -L/libMarmoteCore.a -I/usr/include/python3.6

sparseMatrix_wrap : sparseMatrix_wrap1
	$g -shared  sparseMatrix_wrap.o sparseMatrix.o -L/libMarmoteMDP.a  -L/libMarmoteCore.a   -o _sparseMatrix.so -lstdc++

sparseMatrix_wrap1 : sparseMatrix_wrap2 
	$g -fpic -c  sparseMatrix_wrap.cxx header/sparseMatrix.h -L/libMarmoteMDP.a -L/libMarmoteCore.a -I/usr/include/python3.6
 
marmoteInterval_wrap : marmoteInterval_wrap1
	$g -shared  marmoteInterval_wrap.o marmoteInterval.o -L/libMarmoteMDP.a  -L/libMarmoteCore.a   -o _marmoteInterval.so -lstdc++

marmoteInterval_wrap1 : marmoteInterval_wrap2 
	$g -fpic -c  marmoteInterval_wrap.cxx header/marmoteInterval.h -L/libMarmoteMDP.a -L/libMarmoteCore.a -I/usr/include/python3.6

totalRewardMDP_wrap2 : 
	swig -c++ -python totalRewardMDP.i

sparseMatrix_wrap2 :
	swig -c++ -python sparseMatrix.i

marmoteInterval_wrap2 :
	swig -c++ -python marmoteInterval.i

sspLib: sspLib.i
	swig -c++ -python sspLib.i


