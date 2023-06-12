CC = clang++ # or g++
 
CFLAGS  = -g -Wall
 
  
all: permute
 
permute: permute.cc
	$(CC) $(CFLAGS) -pthread -o permute permute.cc
 
clean:
	-rm) permute