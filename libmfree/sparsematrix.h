#ifndef __SPARSEMATRIX_H_INCLUDED__
#define __SPARSEMATRIX_H_INCLUDED__

#include <common_lib.h>
#include <ostream>
#include "base.h"
#include <mpi.h>

namespace sparse
{

class CSparseVector //non-distributable sparse vector,
					//implemented as dynamic arrays of indices and values, sorted by indices
{
public:
	typedef double T; //in case we will later templatize this, we use T everywhere instead of double

protected:

	int length;		//total length of vector; a full vector would have nnz == length
	int rowIndex;	//tells which row of total matrix this vector is

	int nnz;		//number of nonzeros (i.e., number of elements whose values are not implicitly 0)
	int space;		//allocated space in dynamic arrays buf and indices

	T *buf;			//array of nonzero element values, sorted by their indices
	int *indices;	//sorted array of nonzero element indices (zero-based)

	int findElement(int index) const
	//searches for element index in array indices, using bisection
	//if found, its position is returned; -1 otherwise
	{
		assert(index >= 0 && index < length);
		int lower = 0, upper = nnz-1, middle;
		while (lower <= upper)
		{
			middle = (lower+upper)/2;
			if (index == indices[middle])
				return middle;
			else if (index < indices[middle])
				upper = middle-1;
			else
				lower = middle+1;
		}
		return -1;
	}

	T &access(int index)
	//private function used to find and access the value of element index
	//adds store for element if it is currently implicit zero
	{
		assert(index >= 0 && index < length);
		int position = findElement(index);
		if (position != -1)
			//return reference to existing value
			return buf[position];
		else
		{
			//add store for element and initialize it to 0

			if (nnz == space) //increase allocated space if necessary
				reallocTable(int(space*1.5));

			int lower = 0, upper = nnz-1, middle;
			while (lower <= upper) //find place to insert element
			{
				middle = (lower+upper)/2;
				if (index < indices[middle])
					upper = middle-1;
				else
					lower = middle+1;
			}
			if (nnz == 0)
			{
				buf[0] = 0;
				indices[0] = index;
				nnz = 1;
				return buf[0];
			}
			assert(lower <= nnz);
			assert(lower >= 0);
			assert(indices[lower-1] < index);
			assert(lower == nnz || indices[lower] > index);
			for (int i = nnz-1; i >= lower; i--) //shift elements with higher indices to make room
			{
				buf[i+1] = buf[i];
				indices[i+1] = indices[i];
			}

			//insert a 0 value & return reference to it
			nnz++;
			buf[lower] = 0;
			indices[lower] = index;
			return buf[lower];
		}
	}

	T accessConst(const int index) const
	//private function used to find and access the value of element index
	//if element is implicit 0, just returns 0 instead of adding store
	//not for left-side of assignment; is a const function
	{
		assert(index >= 0 && index < length);
		int position = findElement(index);
		if (position != -1)
			return buf[position];
		else
			return 0;
	}

public:
	CSparseVector()
	//default constructor creates an invalid object
	{
		length = -1;
		rowIndex = -1;
		nnz = space = 0;
		buf = NULL;
		indices = NULL;
	}

	CSparseVector(int _length, int _space = 10)
	{
		buf = NULL;
		indices = NULL;
		initialize(_length, _space);
	}

	virtual ~CSparseVector()
	{
		if (buf) delete []buf;
		if (indices) delete []indices;
	}

	void initialize(int _length, int _space = 10)
	//initializes fields and allocates _space elements in dynamic tables
	{
		length = _length;
		nnz = space = 0;
		rowIndex = -1;
		buf = NULL;
		indices = NULL;
		reallocTable(_space);
	}

	void reallocTable(int newSpace);
	//changes the allocated space for dynamic tables
	//NOTE: new space can be larger or smaller (or the same) than old

	int getLength() const
	{
		return length;
	}

	int getNnz() const
	{
		return nnz;
	}

	bool isImplicitZero(int index) const
	{
		return (findElement(index) == -1);
	}

	T operator()(int index) const // () returns value; [] returns reference and can thus be on left of =
	{
		return accessConst(index);
	}

	T &operator[](int index) // () returns value; [] returns reference and can thus be on left of =
	{
		return access(index);
	}

	void operator=(const CSparseVector &other)
	{
		reallocTable(other.space);
		length = other.length;
		rowIndex = other.rowIndex;
		nnz = other.nnz;
		for (int i = 0; i < nnz; i++)
		{
			buf[i] = other.buf[i];
			indices[i] = other.indices[i];
		}
	}

	int getRowIndex() const
	{
		return rowIndex;
	}

	void setRowIndex(int _rowIndex)
	{
		assert(_rowIndex >= 0);
		rowIndex = _rowIndex;
	}

	void setElemToZero(int index)
	//makes element index an implicit zero (that is, removes its store from dynamic tables)
	{
		int position = findElement(index);
		if (position != -1)
		{
			nnz--;
			for (int i = position; i < nnz; i++)
			{			
				buf[i] = buf[i+1];
				indices[i] = indices[i+1];
			}
		}
	}

	void appendValue(int index, T value)
	//adds value to element index
	//constant time if currently the vector ends either with index or with something smaller
	//used to add n pre-sorted contributions in O(n) time
	{
		if (nnz > 0 && index < indices[nnz-1])
			access(index) += value;
		else
		{
			if (nnz == space)
				reallocTable(int(space*1.5));
			if (nnz == 0 || index > indices[nnz-1])
			{
				nnz++;
				buf[nnz-1] = 0;
				indices[nnz-1] = index;
			}
			buf[nnz-1] += value;
		}
	}

	void sendTo(int dest) const
	//sends the vector to processor dest
	//the protocol consists of the following messages:
	//	a. nnz (1 int, tag 182)
	//	b. buf (nnz doubles, tag 183)
	//	c. indices (nnz ints, tag 184)
	//	d. rowIndex (1 int, tag 185)
	{
		MPI_Send((void*)&nnz, 1, MPI_INT, dest, 182, MPI_COMM_WORLD);
		if (nnz > 0)
		{
			MPI_Send(buf, nnz, MPI_DOUBLE, dest, 183, MPI_COMM_WORLD);
			MPI_Send(indices, nnz, MPI_INT, dest, 184, MPI_COMM_WORLD);
		}
		MPI_Send((void*)&rowIndex, 1, MPI_INT, dest, 185, MPI_COMM_WORLD);
	}

	void recvFrom(int source)
	//recieves a vector from processor source, discarding the old content
	//for protocol, see sentTo()
	{
		MPI_Status status;
		MPI_Recv(&nnz, 1, MPI_INT, source, 182, MPI_COMM_WORLD, &status);
		if (space < nnz)
			reallocTable(nnz);
		if (nnz > 0)
		{
			MPI_Recv(buf, nnz, MPI_DOUBLE, source, 183, MPI_COMM_WORLD, &status);
			MPI_Recv(indices, nnz, MPI_INT, source, 184, MPI_COMM_WORLD, &status);
		}
		MPI_Recv(&rowIndex, 1, MPI_INT, source, 185, MPI_COMM_WORLD, &status);
	}

	void writeToMFile(FILE *f, char *varName, int row) const
	//writes vector to Matlab script in format:
	//  varName(row,col)=value;
	{
		for (int i = 0; i < nnz; i++)
			fprintf(f, "%s(%d,%d)=%.15lf;\n", varName, row+1, indices[i]+1, buf[i]);//+1 because matlab indices start with 1
	}
}; //end class CSparseVector


class CSparseVectorInput
//contribution to a vector
{
public:
	int col;		//column index
	double value;	//contribution value; values with the same col are added together

	CSparseVectorInput() {}
	CSparseVectorInput(int _col, double _value)
		: col(_col), value(_value) {}
}; //end class CSparseVectorInput


class CSparseMatrixInput
//contribution to a matrix
//includes info to which physical row (row) this goes
//as well as which row index (rowIndex) it corresponds to
{
public:
	int row;		//physical row in an object of class CMatrix that this contribution should contribute to
	int rowIndex;	//index of row in a global matrix that row corresponds to; different CSparseMatrixInputs with same row should also have same rowIndex
	int col;		//column index
	double value;	//contribution value; values with the same row and col are added together

	CSparseMatrixInput() {}
	CSparseMatrixInput(int _row, int _rowIndex, int _col, double _value)
		: row(_row), rowIndex(_rowIndex), col(_col), value(_value) {}
}; //end class CSparseMatrixInput


inline int compareMatrixInputsByCol(const void *input1, const void *input2)
//callback for qsort
{
	CSparseMatrixInput *in1 = (CSparseMatrixInput*)input1;
	CSparseMatrixInput *in2 = (CSparseMatrixInput*)input2;
	if (in1->col == in2->col)
		return 0;
	else if (in1->col < in2->col)
		return -1;
	else
		return 1;
}

inline int compareVectorInputsByCol(const void *input1, const void *input2)
//callback for qsort
{
	CSparseVectorInput *in1 = (CSparseVectorInput*)input1;
	CSparseVectorInput *in2 = (CSparseVectorInput*)input2;
	if (in1->col == in2->col)
		return 0;
	else if (in1->col < in2->col)
		return -1;
	else
		return 1;
}

class CSparseMatrix //distributable sparse matrix,
					//implemented as a full vector (static array) of CSparseVectors
					//distributed matrix means that a CSparseMatrix object only stores some rows of total matrix
{
	typedef CSparseVector::T T; //in case we will later templatize this, we use T everywhere instead of double

	int numRows;			//number of rows stored locally (i.e. length of array rows)
	int numCols;			//length of each row
	CSparseVector *rows;	//array of local rows

	CSparseVector &access(int row)
	//private function used to find and access an individual row
	{
		assert(row >= 0 && row < numRows);
		return rows[row];
	}

	const CSparseVector &accessConst(int row) const
	//private function used to find and access an individual row
	//not for left-side of assignment; is a const function
	{
		assert(row >= 0 && row < numRows);
		return rows[row];
	}

	int getIndexOf(int row) const
	//returns rowIndex of the vector stored in position row of table *rows
	{
		assert(row >= 0 && row < numRows);
		return accessConst(row).getRowIndex();
	}

	void setIndexOf(int row, int index)
	//sets rowIndex of the vector stored in position row of table *rows
	{
		assert(row >= 0 && row < numRows);
		assert(getIndexOf(row) == -1 || getIndexOf(row) == index);
		access(row).setRowIndex(index);
	}

public:
	CSparseMatrix(int _numRows, int _numCols, int spacePerRow)
	{
		rows = NULL;
		initialize(_numRows, _numCols, spacePerRow);
	}

	CSparseMatrix()
	//default constructor creates an invalid object
	{
		numRows = numCols = -1;
		rows = NULL;
	}

	void initialize(int _numRows, int _numCols, int spacePerRow)
	//initializes fields and allocates _numRows elements in table *rows,
	//then initializes all rows and allocates space for spacePerRow nonzeros in each
	{
		if (rows) delete []rows;
		numRows = _numRows;
		numCols = _numCols;
		rows = new CSparseVector[numRows];
		for (int i = 0; i < numRows; i++)
			rows[i].initialize(numCols, spacePerRow);
	}

	virtual ~CSparseMatrix()
	{
		if (rows)
			delete []rows;
	}

	int getNumRows() const
	{
		return numRows;
	}

	int getNumCols() const
	{
		return numCols;
	}

	/* usage of operators () and []:
	//() returns value
	double x = A(1,3);
	int nnz1 = A(1).getNnz();
	A(1,3) = 5; //ERROR!
	//[] returns reference and can thus be on left of =
	A[1][3] = 5;
	A[1].setRowIndex(3);
	//NOTE!
	double x = A[1][4]; //if A(1,4) is implicit zero, this will turn it to explicit zero (i.e., allocate store for it)
	*/
	T operator()(int row, int col) const
	{
		return accessConst(row)(col);
	}

	const CSparseVector &operator()(int row) const
	{
		return accessConst(row);
	}

	CSparseVector &operator[](int row)
	{
		return access(row);
	}

	void addInputs(CSparseMatrixInput *inputs, int numInputs)
	//adds a series of contributions to matrix
	//if inputs are sorted by columns, time complexity is linear
	//if there remain any empty rows, should be followed by removeEmptyRows()
	{ 
		assert(inputs);
		for (int i = 0; i < numInputs; i++)
		{
			setIndexOf(inputs[i].row, inputs[i].rowIndex);
			access(inputs[i].row).appendValue(inputs[i].col, inputs[i].value);
		}
	}

	void addInputsSort(CSparseMatrixInput *inputs, int numInputs)
	//sorts inputs by columns, then cals addInputs()
	//if there remain any empty rows, should be followed by removeEmptyRows()
	{ 
		assert(inputs);
		qsort(inputs, numInputs, sizeof(CSparseMatrixInput), compareMatrixInputsByCol);
		addInputs(inputs, numInputs);
	}

	void addInputsSort(CSparseVectorInput *inputs, int numInputs, int row, int rowIndex)
	//adds a series of vector contributions to row
	//pre-sorts them by columns to ensure linear time
	//if there remain any empty rows, should be followed by removeEmptyRows()
	{ 
		assert(inputs);
		assert(row >= 0 && row < numRows);
		setIndexOf(row, rowIndex);
		qsort(inputs, numInputs, sizeof(CSparseVectorInput), compareVectorInputsByCol);
		for (int i = 0; i < numInputs; i++)
			access(row).appendValue(inputs[i].col, inputs[i].value);
	}

	void removeEmptyRows()
	//removes any trailing empty (i.e. all zero) rows
	{
		for (int row = numRows-1; row >= 0; row--)
			if (getIndexOf(row) == -1)
				numRows--;
		assert(numRows > 0);
		for (int row = 0; row < numRows; row++)
			assert(getIndexOf(row) != -1);
	}

	void sendToRoot(int root) const
	//sends all rows stored in this object to processor root (which usually gathers them by calling gatherFrom()
	//the protocol consists of the following messages:
	//  a. numRows (1 int, tag 180)
	//	B. sending rows one-by-one using protocol described in CSparseVector::sendTo()
	{
		MPI_Send((void*)&numRows, 1, MPI_INT, root, 180, MPI_COMM_WORLD);
		for (int row = 0; row < numRows; row++)
			accessConst(row).sendTo(root);
	}

	void gatherFrom(const CSparseMatrix &A, int numProcs);
	//gathers all rows on this processor into a full matrix
	//USAGE:
	//if (myRank==0)
	//		wholeA.gatherFrom(A);
	//else
	//		A.sendToRoot(0);
	//contents of root's A are simply copied to wholeA

	void writeToMFile(FILE *f, char *varName) const
	//writes matrix to Matlab script in format:
	//	varName=spalloc(numRows, numCols, 1);
	//	varName(row,col)=value;
	{
		fprintf(f, "%s=spalloc(%d, %d, 1);\n", varName, numRows, numCols);
		for (int row = 0; row < numRows; row++)
			accessConst(row).writeToMFile(f, varName, getIndexOf(row));
	}
}; //end class CSparseMatrix

}; //end namespace sparse

#endif


