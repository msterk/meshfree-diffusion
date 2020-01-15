CC = mpiCC
# MD added some optimization options
CCOPT = -c -DNDEBUG -O3 -funroll-loops -I ../common_lib    -march=native
CCOPT_DEBUG = -c -I ../common_lib
LINKOPT = -L ../common_lib -L ../libmfree -L ../mygsl
LIBRARIES = -lcommon_lib -llibmfree -lm -lmygsl

mfre: main.cpp ../libmfree/liblibmfree.a
	echo WARNING: dependencies not resolved
	${CC} ${CCOPT} main.cpp
	${CC} ${LINKOPT} -o mfre main.o ${LIBRARIES}
clean:
	-rm -f *.o
	-rm -f mfre
