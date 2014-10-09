CC = mpiCC.lam
CCOPT = -c -DNDEBUG -O3 -funroll-loops -I ../common_lib
LINKOPT = -L ../common_lib -L ../libmfree -L /usr/lib/mpich/lib
LIBRARIES = -lcommon_lib -llibmfree -lmpich -lm

fem: main.cpp
	echo WARNING: dependencies not resolved
	${CC} ${CCOPT} main.cpp
	g++-4.4 ${LINKOPT} -o fe main.o ${LIBRARIES}
clean:
	-rm -f *.o
	-rm -f fe
