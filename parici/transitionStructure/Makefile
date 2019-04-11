CFLAGS=-ansi -Wall -pedantic -g

TRANSITIONS=transitionStructure sparseMatrix multiDimHomTransition

OBJECTS=$(addprefix ../lib/, $(addsuffix .o, $(TRANSITIONS) ) )

all: objects

objects: $(OBJECTS)

../lib/%.o: %.cpp
	g++ $(CFLAGS) -c $< -o $@

clean:
	/bin/rm $(OBJECTS)
