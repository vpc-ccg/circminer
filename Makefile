CC          ?= gcc
CXX         ?= g++

INCS        := 
LIBS        := -lz -lm -lpthread
CFLAGS      := -w 
#CXXFLAGS    := -Wall $(INCS) -std=c++14
CXXFLAGS    := -w $(INCS) -std=c++14
LDFLAGS     := 

SRCDIR      := src
LIBDIR      := lib
OBJDIR      := obj
MRSDIR      := $(SRCDIR)/mrsfast

EXE			:= circminer
SRCEXT      := cpp

MYSRC        = circminer.cpp \
               utils.cpp \
               output.cpp \
               filter.cpp \
               match_read.cpp \
               fastq_parser.cpp \
               commandline_parser.cpp \
			   chain.cpp \
			   gene_annotation.cpp \
			   align.cpp \
			   common.cpp \
			   hash_table.cpp \
			   process_circ.cpp \
			   extend.cpp \
			   genome.cpp \
			   edlib.cpp \

MRSOBJ       = $(MRSDIR)/[!base]*.o

SOURCES = $(patsubst %, $(SRCDIR)/%, $(MYSRC))
MYOBJ = $(SOURCES:$(SRCDIR)/%.$(SRCEXT)=$(OBJDIR)/%.o)

MRSLIBDIR = $(LIBDIR)/mrsfast

MRS_CP_FILES := Common.c \
				Common.h \
				RefGenome.c \
				RefGenome.h \

MRS_SRC_FILES := $(addprefix $(MRSDIR)/, $(MRS_CP_FILES))

EDLIB_LIB_HED_PATH = $(LIBDIR)/edlib/edlib/include/edlib.h
EDLIB_LIB_SRC_PATH = $(LIBDIR)/edlib/edlib/src/edlib.cpp

EDLIB_SRC_PATH = $(SRCDIR)/edlib.cpp
EDLIB_HED_PATH = $(SRCDIR)/edlib.h

LOGGER_LIB_HED_PATH = $(LIBDIR)/util-logger/include/logger.h

LOGGER_HED_PATH = $(SRCDIR)/logger.h

.PHONY: all debug verbose profile valgrind build clean cleanobj cleanlib cleanexe dirs 
.PHONY: OPTIMIZE_FLAGS DEBUG_FLAGS VERBOSE_FLAGS PROFILE_FLAGS PROFILE_FLAGS_GP

all: OPTIMIZE_FLAGS build
debug: DEBUG_FLAGS build
verbose: VERBOSE_FLAGS OPTIMIZE_FLAGS build
profile: PROFILE_FLAGS_GP DEBUG_FLAGS OPTIMIZE_FLAGS build
valgrind: OPTIMIZE_FLAGS DEBUG_FLAGS build
build: cleanexe dirs mrsfast $(EXE)
clean: cleanobj cleanlib

$(MRSLIBDIR):
	@echo Please clone the repository with --recursive option!; exit 1;

$(EDLIB_LIB_SRC_PATH):
	@echo Please clone the repository with --recursive option!; exit 1;

$(EDLIB_LIB_HED_PATH):
	@echo Please clone the repository with --recursive option!; exit 1;

$(LOGGER_LIB_HED_PATH):
	@echo Please clone the repository with --recursive option!; exit 1;

$(EDLIB_SRC_PATH): $(EDLIB_LIB_SRC_PATH) $(EDLIB_HED_PATH)
	@cp -v $(EDLIB_LIB_SRC_PATH) $(EDLIB_SRC_PATH)

$(EDLIB_HED_PATH): $(EDLIB_LIB_HED_PATH)
	@cp -v $(EDLIB_LIB_HED_PATH) $(EDLIB_HED_PATH)

$(LOGGER_HED_PATH): $(LOGGER_LIB_HED_PATH)
	@cp -v $(LOGGER_LIB_HED_PATH) $(LOGGER_HED_PATH)

dirs:
	@mkdir -p $(OBJDIR)

$(EXE): $(LOGGER_HED_PATH) $(EDLIB_SRC_PATH) $(MYOBJ)
	$(CXX) $(MYOBJ) $(MRSOBJ) -o $@ ${LDFLAGS} ${LIBS}

$(OBJDIR)/%.o: $(SRCDIR)/%.$(SRCEXT)
	$(CXX) $(CXXFLAGS) $(LIBS) -c $< -o $@

mrsfast: $(MRS_SRC_FILES)
	@$(MAKE) -C $(MRSDIR)

$(MRSDIR)/%: $(MRSLIBDIR)/%
	@cp -v $< $@

cleanobj:
	@rm -fv $(MYOBJ) $(MRSOBJ)

cleanlib:
	@rm -fv $(LOGGER_HED_PATH) $(EDLIB_SRC_PATH) $(EDLIB_HED_PATH) $(MRS_SRC_FILES)

cleanexe:
	@rm -fv $(EXE)

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
