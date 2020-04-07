MRSFAST_VERSION := "3.4.0"
BUILD_DATE := "$(shell date)"

CC ?= gcc

LIBS=-lz -lm -pthread -lpthread
CFLAGS=-fno-pic -DMRSFAST_VERSION=\"$(MRSFAST_VERSION)\" -DBUILD_DATE=\"$(BUILD_DATE)\"

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
build: SSE_FLAGS dirs $(OBJECTS)

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

SSE_FLAGS:
ifeq ($(shell uname -s),Linux)
ifeq ($(with-sse4),no)
		$(shell echo "-DSSE4=0")
else
        	$(eval CFLAGS = $(CFLAGS) \
        	$(shell gv=`gcc -dumpversion`; \
            	    sc=`grep -c "sse4" /proc/cpuinfo`; \
                	echo $$sc.$$gv | awk -F. '{if($$1>0 && $$2>=4 && $$3>=4) print "-DSSE4=1 -msse4.2"; else print "-DSSE4=0"}'))
endif
else
ifeq ($(with-sse4),no)
		$(shell echo "-DSSE4=0")
else
        $(eval CFLAGS = $(CFLAGS) \
        $(shell gv=`gcc -dumpversion`; \
                sc=`sysctl -n machdep.cpu.features | grep -c "SSE4"` ;\
                echo $$sc.$$gv | awk -F. '{if($$1>0 && $$2>=4 && $$3>=4) print "-DSSE4=1 -msse4.2"; else print "-DSSE4=0"}'))
endif
endif
