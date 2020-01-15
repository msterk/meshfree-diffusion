CC = mpiCC
# MD: added optimization options:
CCOPT = -c -DNDEBUG -O3 -funroll-loops -I ../common_lib    -march=native

mygsl: myblas.o mygsl.o
	ar rc libmygsl.a myblas.o mygsl.o
	rm -f *.o
	ranlib libmygsl.a
	
myblas.o: myblas.c gsl_linalg.h
	${CC} ${CCOPT} myblas.c

mygsl.o: mygsl.c gsl_linalg.h
	${CC} ${CCOPT} mygsl.c
	
clean:
	-rm -f *.o
	-rm -f libmygsl.a
