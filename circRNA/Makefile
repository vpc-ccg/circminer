all: build
build: bwalib circRNA clean

CC          ?= gcc
CXX         ?= g++

SRCDIR      := src
BWADIR      := $(SRCDIR)/bwa

INCS        := -I$(BWADIR)
LIBS        := -lz
#LIBS        := -lz -lm -lpthread
CFLAGS      := -w -O4
CXXFLAGS    := -w -O4 $(INCS)
LDFLAGS     := 

MYOBJ        = $(SRCDIR)/circRNA.o \
               $(SRCDIR)/bwt.o \
               $(SRCDIR)/fastq_parser.o \

BWAOBJ       = $(BWADIR)/*.o

circRNA: $(MYOBJ)
	$(CXX) -w $(MYOBJ) $(BWAOBJ) -o $@ ${LDFLAGS} ${LIBS}

bwalib:
	@$(MAKE) -C $(BWADIR)

clean:
	@rm -f $(MYOBJ)

