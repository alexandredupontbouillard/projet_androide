CFLAGS=-ansi -Wall -fPIC -O3
GPPCOMPILER=g++

INCLUDEDIR=/header
MarmoteCoreDIR=./
MarmoteMDPDIR=./

#Si on a besoin de qu'une seule librairie
LIBRARY=$(addprefix -l, MarmoteCore)

DISTRIBS=Distribution diracDistribution  discreteDistribution geometricDistribution
TRANSITIONS=sparseMatrix transitionStructure
SETS=marmoteSet marmoteBox marmoteInterval 
MDP=solutionMDP feedbackSolutionMDP genericMDP totalRewardMDP 
UTILS = alglin
OBJECTS=\
	$(addsuffix .o, $(DISTRIBS) )  \
	$(addsuffix .o, $(TRANSITIONS) )  \
	$(addsuffix .o, $(SETS) ) \
	$(addsuffix .o, $(MDP) ) \
	$(addsuffix .o, $(UTILS))

 


all: $(DISTRIBS) $(TRANSITIONS) $(SETS) $(MDP) $(UTILS) libraries generate_wrapper compile_wrapper generate_module

## le terme $< represente le premier fichier de dependance celui juste apres les deux points
Distribution: Distribution.cpp
	$(GPPCOMPILER) $(CFLAGS) -c $<  -I$(INCLUDEDIR)

diracDistribution: diracDistribution.cpp
	$(GPPCOMPILER) $(CFLAGS) -c $< -I$(INCLUDEDIR)
	
discreteDistribution: discreteDistribution.cpp
	$(GPPCOMPILER) $(CFLAGS) -c $< -I$(INCLUDEDIR)

geometricDistribution: geometricDistribution.cpp
	$(GPPCOMPILER) $(CFLAGS) -c $< -I$(INCLUDEDIR)

	
sparseMatrix : sparseMatrix.cpp
	$(GPPCOMPILER) $(CFLAGS) -c $< -I$(INCLUDEDIR)

transitionStructure: transitionStructure.cpp
	$(GPPCOMPILER) $(CFLAGS) -c $< -I$(INCLUDEDIR)

marmoteSet: marmoteSet.cpp
	$(GPPCOMPILER) $(CFLAGS) -c $< -I$(INCLUDEDIR)
	
marmoteBox: marmoteBox.cpp
	$(GPPCOMPILER) $(CFLAGS) -c $< -I$(INCLUDEDIR)
	
marmoteInterval: marmoteInterval.cpp
	$(GPPCOMPILER) $(CFLAGS) -c $< -I$(INCLUDEDIR)

solutionMDP : solutionMDP.cpp
	$(GPPCOMPILER) $(CFLAGS) -c $< -I$(INCLUDEDIR)

genericMDP:genericMDP.cpp
	$(GPPCOMPILER) $(CFLAGS) -c $< -I$(INCLUDEDIR)

feedbackSolutionMDP: feedbackSolutionMDP.cpp
	$(GPPCOMPILER) $(CFLAGS) -c $< -I$(INCLUDEDIR)

totalRewardMDP : totalRewardMDP.cpp
	$(GPPCOMPILER) $(CFLAGS) -c $< -I$(INCLUDEDIR)

alglin : alglin.c
	$(GPPCOMPILER) $(CFLAGS) -c $< -I$(INCLUDEDIR)
#$@ represente le nom de la cible 
testSP: testSP.cpp
	$(GPPCOMPILER) $(CFLAGS) $< -o $@ -I$(INCLUDEDIR) -L$(MarmoteCoreDIR) $(LIBRARY) 


objects: $(OBJECTS)

libraries: $(OBJECTS)
	ar crv libMarmoteCore.a $(OBJECTS)	


generate_wrapper: 
	swig -Wall -python -c++ pyMarmoteMDP.i

#inutile 
compile_file: sparseMatrix.cpp
	$(GPPCOMPILER) $(CFLAGS) -c $< -I/usr/include/python3.6 -I$(INCLUDEDIR)

compile_wrapper: pyMarmoteMDP_wrap.cxx
	$(GPPCOMPILER) $(CFLAGS) -c $< -I/usr/include/python3.6 -I$(INCLUDEDIR)

generate_module: pyMarmoteMDP_wrap.o
	g++ -shared pyMarmoteMDP_wrap.o -o _pyMarmoteMDP.so -L$(MarmoteCoreDIR) $(LIBRARY)

clean:
	@/bin/rm -rvf $(OBJECTS)  testSP libMarmoteCore.a
	@/bin/rm -rvf *.o *.py *.cxx
