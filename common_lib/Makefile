CC = mpiCC -I .
# MD: added optimization options
CCOPT = -c -O3 -funroll-loops   -march=native

common_lib: inifiler.o averager.o common_lib.o table.o
	ar rc libcommon_lib.a averager.o common_lib.o inifiler.o table.o
	rm -f *.o
	ranlib libcommon_lib.a

inifiler.o: inifiler.cpp inifiler.h
	${CC} ${CCOPT} inifiler.cpp
	 
averager.o: averager.cpp averager.h
	${CC} ${CCOPT} averager.cpp

common_lib.o: common_lib.cpp common_lib.h vector.h
	${CC} ${CCOPT} common_lib.cpp	

table.o: table.cpp table.h loaddistribution.h
	${CC} ${CCOPT} table.cpp

clean:
	rm -f *.o
	rm -f libcommon_lib.a
