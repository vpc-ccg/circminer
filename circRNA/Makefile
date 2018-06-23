all: build
build: bwalib circRNA clean

CC          ?= gcc
CXX         ?= g++

SRCDIR      := src
BWADIR      := $(SRCDIR)/bwa

INCS        := -I$(BWADIR)
LIBS        := -lz
#LIBS        := -lz -lm -lpthread
CFLAGS      := -w -O3 
CXXFLAGS    := -w -O3 $(INCS) -g
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

circRNA: $(MYOBJ)
	$(CXX) -w $(MYOBJ) $(BWAOBJ) -o $@ ${LDFLAGS} ${LIBS}

bwalib:
	@$(MAKE) -C $(BWADIR)

clean:
	@rm -f $(MYOBJ)

