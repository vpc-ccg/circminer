all: OPTIMIZE_FLAGS build
debug: DEBUG_FLAGS build
profile: PROFILE_FLAGS build
build: bwalib mrsfast circRNA cleanobj

CC          ?= gcc
CXX         ?= g++

SRCDIR      := src
BWADIR      := $(SRCDIR)/bwa
MRSDIR      := $(SRCDIR)/mrsfast

INCS        := -I$(BWADIR)
#LIBS        := -lz
LIBS        := -lz -lm -lpthread
CFLAGS      := -w 
CXXFLAGS    := -w $(INCS)
LDFLAGS     := 

MYOBJ        = $(SRCDIR)/circRNA.o \
               $(SRCDIR)/bwt.o \
               $(SRCDIR)/filter.o \
               $(SRCDIR)/match_read.o \
               $(SRCDIR)/fastq_parser.o \
               $(SRCDIR)/commandline_parser.o \
			   $(SRCDIR)/chain.o \
			   $(SRCDIR)/gene_annotation.o \
			   $(SRCDIR)/align.o \
			   $(SRCDIR)/fragment_list.o \
			   $(SRCDIR)/common.o

BWAOBJ       = $(BWADIR)/*.o
#MRSOBJ       = $(MRSDIR)/*.o
MRSOBJ       = $(MRSDIR)/[!base]*.o

circRNA: $(MYOBJ)
	$(CXX) -w $(MYOBJ) $(BWAOBJ) $(MRSOBJ) -o $@ ${LDFLAGS} ${LIBS}

bwalib:
	@$(MAKE) -C $(BWADIR)

mrsfast:
	@$(MAKE) -C $(MRSDIR)

cleanobj:
	@rm -f $(MYOBJ)

clean:
	@rm -f $(MYOBJ) $(BWAOBJ) $(MRSOBJ)

OPTIMIZE_FLAGS:
	$(eval CFLAGS = $(CFLAGS) -O3)
	$(eval CXXFLAGS = $(CXXFLAGS) -O3)

DEBUG_FLAGS:
	$(eval CFLAGS = $(CFLAGS) -g)
	$(eval CXXFLAGS = $(CXXFLAGS) -g)

PROFILE_FLAGS:
	$(eval LIBS = $(LIBS) -lprofiler)
