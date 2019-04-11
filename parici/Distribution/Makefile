CFLAGS=-ansi -Wall -pedantic -g

DISTRIBS=bernoulliDistribution Distribution \
	uniformDistribution diracDistribution \
	exponentialDistribution discreteDistribution  \
	uniformDiscreteDistribution geometricDistribution \
	poissonDistribution

OBJECTS=$(addprefix ../lib/, $(addsuffix .o, $(DISTRIBS) ) )

all: objects libraries

objects: $(OBJECTS)

../lib/%.o: %.cpp
	g++ $(CFLAGS) -c $< -o $@

libraries:
	ar rucv ../lib/distribution.a $(OBJECTS)

clean:
	/bin/rm $(OBJECTS) ../lib/distribution.a
