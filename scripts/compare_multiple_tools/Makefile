all: OPTIMIZE_FLAGS build
build: cleanexe mydiff 
clean: cleanobj

CXX         ?= g++

INCS        := 
CXXFLAGS    := -Wall $(INCS) -std=c++11
LDFLAGS     := 

MYOBJ        = diff.o 


mydiff: $(MYOBJ)
	$(CXX) -w $(MYOBJ) -o $@ ${LDFLAGS} ${LIBS}

cleanobj:
	@rm -f $(MYOBJ)

cleanexe:
	@rm -f mydiff

OPTIMIZE_FLAGS:
	$(eval CXXFLAGS = $(CXXFLAGS) -O3)

