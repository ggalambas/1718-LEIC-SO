CFLAGS = -g -Wall -pedantic
CC = gcc
ARGS = 6 10 10 0 0 100 3 0 "ola" 1
# N tEsq tSup tDir tInf iter trab maxD fichS periodoS

heatSim: main.o matrix2d.o util.o
	$(CC) $(CFLAGS) -pthread -o heatSim main.o util.o matrix2d.o

main.o : main.c
	$(CC) $(CFLAGS) -c main.c

util.o: util.c util.h
	$(CC) $(CFLAGS) -c util.c

matrix2d.o: matrix2d.c matrix2d.h
	$(CC) $(CFLAGS) -c matrix2d.c

clean:
	rm -f *.o heatSim

zip:
	zip heatSim.zip main.c util.c util.h matrix2d.c matrix2d.h Makefile

run:
	./heatSim $(ARGS)

gdb:
	gdb --args ./heatSim $(ARGS)

val:
	valgrind --leak-check=yes ./heatSim $(ARGS)
