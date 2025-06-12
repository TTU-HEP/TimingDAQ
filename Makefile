CXX = $(shell root-config --cxx)
LD = $(shell root-config --ld)
INC = $(shell pwd)

LDFLAGS := $(shell root-config --glibs)
CPPFLAGS := $(shell root-config --cflags) -I$(INC)/include
#CPPFLAGS += -g -std=c++14
CPPFLAGS += -g
#CPPFLAGS += -g -fsanitize=address -Wall -Wextra -Wno-sign-compare

ifeq ($(shell uname), Darwin)
	CPPFLAGS += -rpath $(shell root-config --prefix)/lib
endif

TARGETS = NetScopeStandaloneDat2Root
SRC = src/Configuration.cc src/DatAnalyzer.cc 

all : $(TARGETS)

$(TARGETS) : %Dat2Root : $(SRC:.cc=.o) src/%Analyzer.o ./%Dat2Root.cc
	@echo Building $@
	$(LD) $(CPPFLAGS) -o $@ $^ $(LDFLAGS)

%.o : %.cc
	@echo $@
	$(CXX) $(CPPFLAGS) -o $@ -c $<
clean :
	rm -rf *.o src/*.o $(TARGETS) *~ *.dSYM
