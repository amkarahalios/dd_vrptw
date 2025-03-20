SYSTEM     = x86-64_linux
LIBFORMAT  = static_pic

#------------------------------------------------------------
#
# When you adapt this makefile to compile your CPLEX programs
# please copy this makefile and set CPLEXDIR and CONCERTDIR to
# the directories where CPLEX and CONCERT are installed.
#
#------------------------------------------------------------

CPLEXDIR      = /opt/ibm/ILOG/CPLEX_Studio221/cplex
CONCERTDIR    = /opt/ibm/ILOG/CPLEX_Studio221/concert
VRPTWSEPDIR    = $(HOME)/dd_vrptw/cvrpsep
CLIQUERDIR    = $(HOME)/dd_vrptw/Cliquer/src

# ---------------------------------------------------------------------
# Compiler selection
# ---------------------------------------------------------------------

CCC = g++ -O3
#CCC = g++ -pg -g

# ---------------------------------------------------------------------
# Compiler options
# ---------------------------------------------------------------------

CCOPT = -m64 -fPIC -fno-strict-aliasing -fexceptions -DNDEBUG -DIL_STD

# ---------------------------------------------------------------------
# Link options and libraries
# ---------------------------------------------------------------------

CPLEXBINDIR   = $(CPLEXDIR)/bin/$(BINDIST)
CPLEXLIBDIR   = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CPLEXVRPTWSEPDIR = $(VRPTWSEPDIR)/obj/
CLIQUERLIBDIR = $(CLIQUERDIR)/lib/
CONCERTLIBDIR = $(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)

# For dynamic linking
CPLEXBINDIR   = $(CPLEXDIR)/bin/$(SYSTEM)
CPLEXLIB      = cplex$(dynamic:yes=1290)

CCLNDIRS  = -L$(CPLEXLIBDIR) -L$(CONCERTLIBDIR) $(dynamic:yes=-L$(CPLEXBINDIR)) -L$(CPLEXVRPTWSEPDIR) -L$(CLIQUERLIBDIR)
CLNDIRS   = -L$(CPLEXLIBDIR) $(dynamic:yes=-L$(CPLEXBINDIR))
CCLNFLAGS = -lconcert -lilocplex -l$(CPLEXLIB) -lm -lpthread -ldl -lcvrpsep -lcliquer
CLNFLAGS  = -l$(CPLEXLIB) -lm -lpthread -ldl

CONCERTINCDIR = $(CONCERTDIR)/include
CPLEXINCDIR   = $(CPLEXDIR)/include

CCFLAGS = $(CCOPT) -I$(CPLEXINCDIR) -I$(CONCERTINCDIR) -I$(VRPTWSEPDIR) -I$(CLIQUERDIR)

all: solver

solver: solvevrptw.o vrptwdecisiondiagram.o vrptwcolgen.o vrptwddsolver.o vrptw.o
	$(CCC) $(CCFLAGS) $(CCLNDIRS) -o solver solvevrptw.o vrptwdecisiondiagram.o vrptwcolgen.o vrptwddsolver.o vrptw.o $(CCLNFLAGS)

solvevrptw.o: solvevrptw.cpp
	$(CCC) -c $(CCFLAGS) solvevrptw.cpp -o solvevrptw.o

vrptwdecisiondiagram.o: vrptwdecisiondiagram.cpp
	$(CCC) -c $(CCFLAGS) vrptwdecisiondiagram.cpp -o vrptwdecisiondiagram.o

vrptwcolgen.o: vrptwcolgen.cpp
	$(CCC) -c $(CCFLAGS) vrptwcolgen.cpp -o vrptwcolgen.o

vrptwddsolver.o: vrptwddsolver.cpp
	$(CCC) -c $(CCFLAGS) vrptwddsolver.cpp -o vrptwddsolver.o

vrptw.o: vrptw.cpp
	$(CCC) -c $(CCFLAGS) vrptw.cpp -o vrptw.o

clean:
	rm -f *.o solver
