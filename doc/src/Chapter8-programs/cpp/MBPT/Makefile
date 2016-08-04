# http://mrbook.org/blog/tutorials/make/
# http://www.cs.colby.edu/maxwell/courses/tutorials/maketutor/

CC=g++
#CFLAGS= -Wall -g -fbounds-check -ftrapv -pedantic-errors #-o3 #-Wextra -pedantic -std=c++11
CFLAGS= -Wall -O3 -pg -static
CXXFLAGS= -fopenmp
#LIBS= -llapack -lblas -lgfortran -lm
LIBS= -lopenblas# -lm
EXECUTABLES = MBPT

O_FILES = $(patsubst %.cpp, %.o, $(wildcard *.cpp))

$(EXECUTABLES): $(O_FILES)
	$(CC) $(CXXFLAGS) $(CFLAGS) -o $@ $^ $(LIBS)

%.o: %.cpp  
	$(CC) $(CXXFLAGS) $(CFLAGS) -c -o $@ $< 

clean:
	rm *.o *~ *# 
