CC=g++
CCOPT = -m64 -O -fPIC -fexceptions -DNDEBUG -DIL_STD
CFLAGS = -c -Wall -g $(CCOPT)
LDFLAGS = -lz -lm -g 
SOURCES = main.cc common.cc fasta.cc geneannotation.cc
OBJECTS = $(SOURCES:.cpp=.o)
EXECUTABLE = circ

all: $(SOURCES) $(EXECUTABLE)
	rm -rf *.o

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) -o $@ $(LDFLAGS)

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@ 
clean:
	rm -f *.o *~ 
