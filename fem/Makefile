CC = mpiCC
CCOPT = -c -DNDEBUG -O3 -funroll-loops -I ../common_lib
LINKOPT = -L ../common_lib -L ../libmfree
LIBRARIES = -lcommon_lib -llibmfree  -lm

fem: main.cpp
	echo WARNING: dependencies not resolved
	${CC} ${CCOPT} main.cpp
	${CC} ${LINKOPT} -o fe main.o ${LIBRARIES}
clean:
	-rm -f *.o
	-rm -f fe
