CFLAGS = -std=c++11 -O3 -DNDEBUG -fopenmp

all: marvin marvin_prep

HTSDIR=htslib-1.2.1
include $(HTSDIR)/htslib.mk
HTSLIB = $(HTSDIR)/libhts.a
IFLAGS = -I$(HTSDIR)  -I./
LFLAGS = -lz -lm
CXX=g++

default: CFLAGS = -std=c++11 -O3 -DNDEBUG -fopenmp
default: all

debug: CFLAGS = -std=c++11 -g -O1 -fopenmp
debug: all

profile: CFLAGS = -std=c++11 -pg -O3 -fopenmp
profile: all

marvin_prep: marvin_prep.cpp $(HTSLIB)
	$(CXX) $(CFLAGS) -o marvin_prep marvin_prep.cpp $(IFLAGS) $(HTSLIB) $(LFLAGS) 
marvin: marvin.cpp $(HTSLIB)
	$(CXX) $(CFLAGS) -o marvin marvin.cpp  $(IFLAGS) $(HTSLIB) $(LFLAGS) 	
clean:
	rm marvin marvin_prep
