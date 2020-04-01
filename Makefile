all: OPTIMIZE_FLAGS build
debug: DEBUG_FLAGS build
verbose: VERBOSE_FLAGS OPTIMIZE_FLAGS build
profile: PROFILE_FLAGS_GP DEBUG_FLAGS OPTIMIZE_FLAGS build
valgrind: OPTIMIZE_FLAGS DEBUG_FLAGS build
#build: cleanexe $(EDLIB_SRC_PATH) mrsfast circminer cleanobj
build: cleanexe $(EDLIB_SRC_PATH) mrsfast circminer 
clean: cleanobj

CC          ?= gcc
CXX         ?= g++

SRCDIR      := src
LIBDIR      := lib
MRSDIR      := $(SRCDIR)/mrsfast

INCS        := 
#LIBS        := -lz
LIBS        := -lz -lm -lpthread
CFLAGS      := -w 
CXXFLAGS    := -Wall $(INCS) -std=c++11
LDFLAGS     := 

MYOBJ        = $(SRCDIR)/circminer.o \
               $(SRCDIR)/utils.o \
               $(SRCDIR)/output.o \
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
			   $(SRCDIR)/extend.o \
			   $(SRCDIR)/genome.o \
			   $(SRCDIR)/edlib.o \

MRSOBJ       = $(MRSDIR)/[!base]*.o


EDLIB_LIB_HED_PATH = $(LIBDIR)/edlib/edlib/include/edlib.h
EDLIB_LIB_SRC_PATH = $(LIBDIR)/edlib/edlib/src/edlib.cpp

EDLIB_SRC_PATH = $(SRCDIR)/edlib.cpp
EDLIB_HED_PATH = $(SRCDIR)/edlib.h

$(EDLIB_LIB_SRC_PATH):
		@echo Please clone the repository with --recursive option!; exit 1;

$(EDLIB_LIB_HED_PATH):
		@echo Please clone the repository with --recursive option!; exit 1;

$(EDLIB_SRC_PATH): $(EDLIB_LIB_SRC_PATH) $(EDLIB_HED_PATH)
		@cp -v $(EDLIB_LIB_SRC_PATH) $(EDLIB_SRC_PATH)

$(EDLIB_HED_PATH): $(EDLIB_LIB_HED_PATH)
		@cp -v $(EDLIB_LIB_HED_PATH) $(EDLIB_HED_PATH)

circminer: $(MYOBJ)
	$(CXX) -w $(MYOBJ) $(BWAOBJ) $(MRSOBJ) -o $@ ${LDFLAGS} ${LIBS}

mrsfast:
	@$(MAKE) -C $(MRSDIR)

cleanobj:
	@rm -f $(MYOBJ) $(MRSOBJ)

cleanexe:
	@rm -f circminer

OPTIMIZE_FLAGS:
	$(eval CFLAGS = $(CFLAGS) -O3)
	$(eval CXXFLAGS = $(CXXFLAGS) -O3)

DEBUG_FLAGS:
	$(eval CFLAGS = $(CFLAGS) -g)
	$(eval CXXFLAGS = $(CXXFLAGS) -g)

VERBOSE_FLAGS:
	$(eval CFLAGS = $(CFLAGS) -DDEBUG=1)
	$(eval CXXFLAGS = $(CXXFLAGS) -DDEBUG=1)

PROFILE_FLAGS:
	$(eval LIBS = $(LIBS) -lprofiler)

PROFILE_FLAGS_GP:
	$(eval CFLAGS = $(CFLAGS) -pg)
	$(eval CXXFLAGS = $(CXXFLAGS) -pg)
	$(eval LIBS = $(LIBS) -pg)
