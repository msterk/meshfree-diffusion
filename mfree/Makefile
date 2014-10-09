CC = mpiCC.lam
CCOPT = -c -DNDEBUG -O3 -funroll-loops -I ../common_lib
CCOPT_DEBUG = -c -I ../common_lib
LINKOPT = -L ../common_lib -L ../libmfree -L ../mygsl -L /usr/lib/mpich/lib
LIBRARIES = -lcommon_lib -llibmfree -lmpich -lm -lmygsl

mfree: main.cpp
	echo WARNING: dependencies not resolved
	${CC} ${CCOPT} main.cpp
	g++-4.4 ${LINKOPT} -o mfre main.o ${LIBRARIES}
clean:
	-rm -f *.o
	-rm -f mfre
