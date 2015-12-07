EXEC   = blend_peaks

OPTIMIZE =  -O2 -Wno-write-strings 


OBJS   = main.o box_collision.o write_shock_catalogues.o timer.o rng.o read_athena_header.o blend_peaks.o kdtree2.o


CXX    = mpicxx
CC     = mpicxx

INCL   = box_collision.hpp shock_data_types.hpp write_shock_catalogues.hpp timer.h rng.h read_athena_header.hpp read_athena_tracers.h spanning_tree_data_types.hpp blend_peaks.hpp kdtree2.hpp

LIBS   = -lm -lgsl -lgslcblas -lmpi  -stdlib=libstdc++

CFLAGS = $(OPTIMIZE) -stdlib=libstdc++
CXXFLAGS = $(OPTIMIZE) -stdlib=libstdc++

$(EXEC): $(OBJS) 
	 $(CXX) $(OBJS) $(LIBS) -o $(EXEC)   
         

$(OBJS): $(INCL) 

.PHONY : clean

clean:
	 rm -f $(OBJS) $(EXEC)

