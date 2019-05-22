CFLAGS=-ansi -Wall -pedantic -g
g=g++

MDP= genericMDP discountedMDP  averageMDP totalRewardMDP
SOLUTION = solutionMDP feedbackSolutionMDP nonStationarySolutionMDP
UTILITES = alglin
INCLUDEDIR = /header
DISTRIBUTION = Distribution discreteDistribution diracDistribution
MARMOTECORE = binarySequence binarySimplex discountedMDP_LPsolve eventMixture marmoteBox marmoteInterval marmoteSet multiDimHomTransition sparseMatrix transitionStructure
WRAPPER = totalRewardMDP_wrap sparseMatrix_wrap marmoteInterval_wrap  

MARMOTECOREO =\
	$(addsuffix .o, $(MARMOTECORE) ) 

DISTRIBUTION0 = \
	$(addsuffix .o, $(DISTRIBUTION) )

OBJECTS=\
	$(addsuffix .o, $(MDP) ) \
	$(addsuffix .o, $(SOLUTION) )  \
	$(addsuffix .o, $(UTILITES) ) \



POINTH=\
	$(addprefix header/,$(addsuffix .h, $(MDP) )) \
	$(addprefix header/,$(addsuffix .h, $(SOLUTION) ))  \
	$(addprefix header/,$(addsuffix .h, $(UTILITES) )) \
	$(addprefix header/,$(addsuffix .h, $(MARMOTECORE) ))


all: $(DISTRIBUTION) $(MDP)  $(SOLUTION) $(UTILITES) $(MARMOTECORE) librairies1 librairies pymarmote

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

diracDistribution : diracDistribution.cpp
	$g $(CFLAGS) -c $<  -I$(INCLUDEDIR)

discreteDistribution : discreteDistribution.cpp
	$g $(CFLAGS) -c $<  -I$(INCLUDEDIR)

Distribution : Distribution.cpp
	$g $(CFLAGS) -c $<  -I$(INCLUDEDIR)


clean:
	/bin/rm -v *.o *.so *.py 

librairies1: $(MARMOTECOREO) $(DISTRIBUTION)
	ar crv libMarmoteCORE.a $(MARMOTECOREO) $(DISTRIBUTION0)

librairies: $(OBJECTS)  $(MARMOTECOREO) $(DISTRIBUTION)
	ar crv libMarmoteMDP.a $(OBJECTS) $(MARMOTECOREO) $(DISTRIBUTION0)


pymarmote : pymarmote_wrap.o
	$g -shared pymarmote_wrap.o  -o _pymarmote.so

pymarmote_wrap.o : pymarmote_wrap.cxx
	$g $(CFLAGS) -c $<  -I/usr/include/python3.6
 

pymarmote_wrap.cxx : 
	swig -c++ -python pymarmote.i
