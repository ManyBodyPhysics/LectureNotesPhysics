CC= g++
CFLAGS= -Wall -Wextra -pedantic -std=c++11 -fopenmp 
LIBS= -lblas -lm 
OBJ_DIR= ../obj/
EXECUTABLES = main.exe

O_FILES = $(patsubst %.cpp, $(OBJ_DIR)%.o, $(wildcard *.cpp))

$(EXECUTABLES): $(O_FILES)
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

$(OBJ_DIR)%.o: %.cpp  
	$(CC) $(CFLAGS) -c -o $@ $< 

clean:
	rm $(OBJ_DIR)*.o *~ *# $(EXECUTABLES) 
 
