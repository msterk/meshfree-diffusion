#include "sparsematrix.h"

namespace sparse
{

void CSparseVector::reallocTable(int newSpace)
//changes the allocated space for dynamic tables
//NOTE: new space can be larger or smaller (or the same) than old
{
	assert(newSpace > 0 && newSpace >= nnz);
	if (space == 0)
	{
		buf = (T*)malloc(newSpace*sizeof(T));
		indices = (int*)malloc(newSpace*sizeof(int));
	} else {
		buf = (T*)realloc(buf, newSpace*sizeof(T));
		indices = (int*)realloc(indices, newSpace*sizeof(int));
	}
	assert(buf);
	assert(indices);
	space = newSpace;
}

void CSparseMatrix::gatherFrom(const CSparseMatrix &A, int numProcs)
//gathers all rows on this processor into a full matrix
//USAGE:
//if (myRank==0)
//		wholeA.gatherFrom(A);
//else
//		A.sendToRoot(0);
//contents of root's A are simply copied to wholeA
{
	MPI_Status status;
	int numRecvdRows = 0;
	//copy root's rows
	for (int row = 0; row < A.getNumRows(); row++)
	{
		access(row) = A(row);
		numRecvdRows++;
	}

	for (int proc = 1; proc < numProcs; proc++)
	{
		int numRowsProc;
		MPI_Recv((void*)&numRowsProc, 1, MPI_INT, proc, 180, MPI_COMM_WORLD, &status);
		for (int rowProc = 0; rowProc < numRowsProc; rowProc++)
			access(numRecvdRows++).recvFrom(proc);
	}
}

}; //end namespace sparse

