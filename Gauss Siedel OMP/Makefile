CC = gcc
FFLAGS = -Wall -std=c99 -Xpreprocessor -fopenmp
LFLAGS = -lgomp -lm
OBJECTS = gs_omp.o



main.exe: $(OBJECTS)
	$(CC) $(LFLAGS) $(OBJECTS) -o main.exe

%.o: %.c
	$(CC) $(FFLAGS) -c $<

clean:
	rm -f $(OBJECTS) *.exe *.out *.txt 
