MRSFAST_VERSION := "3.4.0"

CC ?= gcc

LIBS=-lz -lm -pthread -lpthread
CFLAGS=-fno-pic -DMRSFAST_VERSION=\"$(MRSFAST_VERSION)\"

SRCEXT	:= c
OBJDIR	:= obj
SOURCES	= Sort.c Common.c RefGenome.c HashTable.c
OBJECTS	= $(SOURCES:%.$(SRCEXT)=$(OBJDIR)/%.o) 

.PHONY: all debug profile valgrind build clean dirs
.PHONY: OPTIMIZE_FLAGS DEBUG_FLAGS PROFILE_FLAGS

all: OPTIMIZE_FLAGS build
debug: DEBUG_FLAGS build
profile: PROFILE_FLAGS build
valgrind: OPTIMIZE_FLAGS DEBUG_FLAGS build
build: dirs $(OBJECTS)

dirs:
	@mkdir -p $(OBJDIR)

$(OBJDIR)/%.o: %.$(SRCEXT)
	$(CC) $(CFLAGS) $(LIBS) -c $< -o $@ 

clean:
	@rm -fv $(OBJECTS)
	
DEBUG_FLAGS:
	$(eval CFLAGS = $(CFLAGS) -ggdb)
	$(eval LIBS = $(LIBS) -ggdb)

OPTIMIZE_FLAGS:
	$(eval CFLAGS = $(CFLAGS) -O2)

PROFILE_FLAGS:
	$(eval CFLAGS = $(CFLAGS) -pg -g)
	$(eval LIBS = $(LIBS) -pg -g)

