all: OPTIMIZE_FLAGS build
debug: DEBUG_FLAGS build
profile: PROFILE_FLAGS DEBUG_FLAGS OPTIMIZE_FLAGS build
valgrind: OPTIMIZE_FLAGS DEBUG_FLAGS build
build: mrsfast circRNA cleanobj

CC          ?= gcc
CXX         ?= g++

SRCDIR      := src
MRSDIR      := $(SRCDIR)/mrsfast

INCS        := 
#LIBS        := -lz
LIBS        := -lz -lm -lpthread
CFLAGS      := -w 
CXXFLAGS    := -w $(INCS)
LDFLAGS     := 

MYOBJ        = $(SRCDIR)/circRNA.o \
               $(SRCDIR)/filter.o \
               $(SRCDIR)/match_read.o \
               $(SRCDIR)/fastq_parser.o \
               $(SRCDIR)/commandline_parser.o \
			   $(SRCDIR)/chain.o \
			   $(SRCDIR)/gene_annotation.o \
			   $(SRCDIR)/align.o \
			   $(SRCDIR)/common.o \
			   $(SRCDIR)/hash_table.o \
			   $(SRCDIR)/process_circ.o \
			   $(SRCDIR)/extend.o

#MRSOBJ       = $(MRSDIR)/*.o
MRSOBJ       = $(MRSDIR)/[!base]*.o

circRNA: $(MYOBJ)
	$(CXX) -w $(MYOBJ) $(BWAOBJ) $(MRSOBJ) -o $@ ${LDFLAGS} ${LIBS}

mrsfast:
	@$(MAKE) -C $(MRSDIR)

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
