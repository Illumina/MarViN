HTSLIB_PATH=htslib-1.2.1
HTSLIB = $(HTSLIB_PATH)/libhts.a
IFLAGS = -I$(HTSLIB_PATH) 
LFLAGS = -lz -lm
CXX=g++

CFLAGS = -std=c++11 -O3 -DNDEBUG -fopenmp

default: CFLAGS = -std=c++11 -O3 -DNDEBUG -fopenmp
default: all

debug: CFLAGS = -std=c++11 -g -O1 -fopenmp
debug: all

profile: CFLAGS = -std=c++11 -pg -O3 -fopenmp
profile: all

marvin_prep: marvin_prep.cpp 
	$(CXX) $(CFLAGS) -o marvin_prep marvin_prep.cpp $(IFLAGS) $(HTSLIB) $(LFLAGS) 
marvin: marvin.cpp
	$(CXX) $(CFLAGS) -o marvin marvin.cpp  $(IFLAGS) $(HTSLIB) $(LFLAGS) 	
		
all: marvin marvin_prep

clean:
	rm marvin marvin_prep
