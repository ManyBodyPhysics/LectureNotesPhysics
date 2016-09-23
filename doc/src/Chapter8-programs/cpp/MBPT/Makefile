# http://mrbook.org/blog/tutorials/make/
# http://www.cs.colby.edu/maxwell/courses/tutorials/maketutor/

CFLAGS = -Wall -fopenmp
CC=g++
ifdef DEBUG
	CFLAGS += -g -fbounds-check -ftrapv -pedantic-errors
else
	CFLAGS += -O2
endif

ifdef PROFILE
	CFLAGS += -pg
endif

ifdef OPTIMIZE
	CFLAGS += -O3
endif

LIBS= -lopenblas
EXECUTABLES = MBPT

O_FILES = $(patsubst %.cpp, %.o, $(wildcard *.cpp))

$(EXECUTABLES): $(O_FILES)
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

%.o: %.cpp  
	$(CC) $(CFLAGS) -c -o $@ $< 

clean:
	rm *.o *~ *# 
