CC = mpiCC
CCOPT = -c -DNDEBUG -O3 -funroll-loops -I ../common_lib

libmfree: mfreediffusion.o nodes.o sparsematrix.o
	ar rc liblibmfree.a mfreediffusion.o nodes.o sparsematrix.o
	rm -f *.o
	ranlib liblibmfree.a
	
mfreediffusion.o: mfreediffusion.h mfreediffusion.cpp base.h libmfree.h nodes.h sparsematrix.h
	${CC} ${CCOPT} mfreediffusion.cpp

nodes.o: nodes.cpp nodes.h base.h libmfree.h
	${CC} ${CCOPT} nodes.cpp

sparsematrix.o: sparsematrix.cpp sparsematrix.h base.h libmfree.h nodes.h
	${CC} ${CCOPT} sparsematrix.cpp
	
clean:
	-rm -f *.o
	-rm -f liblibmfree.a
