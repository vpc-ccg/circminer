all: OPTIMIZE_FLAGS build
debug: DEBUG_FLAGS OPTIMIZE_FLAGS build
profile: PROFILE_FLAGS DEBUG_FLAGS OPTIMIZE_FLAGS build
valgrind: OPTIMIZE_FLAGS DEBUG_FLAGS build
build: eval cleanobj

CC          ?= gcc
CXX         ?= g++

SRCDIR      := src
MRSDIR      := $(SRCDIR)/mrsfast

INCS        := 
#LIBS        := -lz
LIBS        := -lz -lm -lpthread
CFLAGS      := -w 
CXXFLAGS    := -w $(INCS) -std=c++14
LDFLAGS     := 

MYOBJ        = eval.o \
               cigar.o 

eval: $(MYOBJ)
	$(CXX) -w $(MYOBJ) $(BWAOBJ) $(MRSOBJ) -o $@ ${LDFLAGS} ${LIBS}

cleanobj:
	@rm -f $(MYOBJ)

clean:
	@rm -f $(MYOBJ) $(MRSOBJ)

OPTIMIZE_FLAGS:
	$(eval CFLAGS = $(CFLAGS) -O3)
	$(eval CXXFLAGS = $(CXXFLAGS) -O3)

DEBUG_FLAGS:
	$(eval CFLAGS = $(CFLAGS) -g -DDEBUG=1)
	$(eval CXXFLAGS = $(CXXFLAGS) -g -DDEBUG=1)

PROFILE_FLAGS:
	$(eval LIBS = $(LIBS) -lprofiler)
