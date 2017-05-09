INCLUDE   = -I. -I/home/saito/local/include -I/usr/include -I/usr/local/include

LIBPATH   = -L. -L/home/saito/local/lib -L/usr/lib -L/usr/local/lib

LDLIBS    = -lm -lboost_program_options -lboost_thread

DEFINE    = 
#DEFINE   += -DNDEBUG

CPPFLAGS += $(INCLUDE) $(DEFINE) -Wall -Wno-sign-compare -O3
#CPPFLAGS += -static 
#CPPFLAGS += -pg

LDFLAGS  += $(LIBPATH)

CXX ?= g++

PROG1 = 3DColocNormalize
OBJ1  = 3DColocNormalize.o Framework3DColocNormalize.o Utility.o \
	Heatmap.o
HEAD1 = Framework3DColocNormalize.hh Utility.hh \
	SparseMatrix.hh Heatmap.hh

PROG2 = 3DColocTest
OBJ2  = 3DColocTest.o Framework3DColocTest.o Utility.o \
	HypothesisTest.o Heatmap.o Bed.o StatCalculator.o
HEAD2 = Framework3DColocTest.hh Utility.hh \
	HypothesisTest.hh SparseMatrix.hh Heatmap.hh Bed.hh StatCalculator.hh

PROG3 = 3DColocSearch
OBJ3  = 3DColocSearch.o Framework3DColocSearch.o Utility.o \
	HypothesisTest.o Heatmap.o Bed.o StatCalculator.o
HEAD3 = Framework3DColocSearch.hh Utility.hh \
	HypothesisTest.hh SparseMatrix.hh Heatmap.hh Bed.hh StatCalculator.hh

PROG4 = 3DColocLook
OBJ4  = 3DColocLook.o Framework3DColocLook.o Utility.o \
	HypothesisTest.o Heatmap.o Bed.o StatCalculator.o
HEAD4 = Framework3DColocLook.hh Utility.hh \
	HypothesisTest.hh SparseMatrix.hh Heatmap.hh Bed.hh StatCalculator.hh

PROG5 = 3DColocDiff
OBJ5  = 3DColocDiff.o Framework3DColocDiff.o Utility.o \
	Heatmap.o
HEAD5 = Framework3DColocDiff.hh Utility.hh \
	SparseMatrix.hh Heatmap.hh

#------------- Make rules -----------------

.cc.o:
	$(CXX) $(CPPFLAGS) -c $< -o $@

.cpp.o:
	$(CXX) $(CPPFLAGS) -c $< -o $@

#------------- Commands -------------------

.PHONY: all clean

all: $(PROG1) $(PROG2) $(PROG3) $(PROG4) $(PROG5) Makefile

clean:
	\rm -f *.o *~ $(PROG1) $(PROG2) $(PROG3) $(PROG4) $(PROG5)
	chmod 644 *.cc *.hh

$(PROG1): $(OBJ1)
	$(CXX) $(CPPFLAGS) -o $(PROG1) $(OBJ1) $(LDFLAGS) $(LDLIBS)
$(OBJ1): $(HEAD1)

$(PROG2): $(OBJ2)
	$(CXX) $(CPPFLAGS) -o $(PROG2) $(OBJ2) $(LDFLAGS) $(LDLIBS)
$(OBJ2): $(HEAD2)

$(PROG3): $(OBJ3)
	$(CXX) $(CPPFLAGS) -o $(PROG3) $(OBJ3) $(LDFLAGS) $(LDLIBS)
$(OBJ3): $(HEAD3)

$(PROG4): $(OBJ4)
	$(CXX) $(CPPFLAGS) -o $(PROG4) $(OBJ4) $(LDFLAGS) $(LDLIBS)
$(OBJ4): $(HEAD4)

$(PROG5): $(OBJ5)
	$(CXX) $(CPPFLAGS) -o $(PROG5) $(OBJ5) $(LDFLAGS) $(LDLIBS)
$(OBJ5): $(HEAD5)
