#ifndef __BASE_H_INCLUDED__
#define __BASE_H_INCLUDED__

#include <string.h>
#include <stdio.h>
#include <assert.h>
#include <common_lib.h>
#include <stdlib.h>

//various replacements for MPI_Barrier; use for debugging; select by uncommenting one version

//barrier that prints its location in the program
#define MPI_BARRIER MPI_Barrier(MPI_COMM_WORLD); if (numProcs > 1 && myRank == 0) { std::cout << "Barrier in line " << __LINE__ << " of file " << __FILE__ << std::endl; } MPI_Barrier(MPI_COMM_WORLD);

//barrier that prints its location in the program & waits for user to press Enter
//#define MPI_BARRIER MPI_Barrier(MPI_COMM_WORLD); if (numProcs > 1 && myRank == 0) { std::cout << "Barrier in line " << __LINE__ << " of file " << __FILE__ << std::endl; getchar(); } MPI_Barrier(MPI_COMM_WORLD);

//just barrier
//#define MPI_BARRIER MPI_Barrier(MPI_COMM_WORLD);

//no barrier
//#define MPI_BARRIER ;

//executes the next statement only on processor with rank R
#define IFRANK(R) int __myRank__; MPI_Comm_rank(MPI_COMM_WORLD, &__myRank__); if (__myRank__==(R))

//prints DATA only on processor with rank R
#define PRINTIFRANK(R, DATA) {int __myRank__; MPI_Comm_rank(MPI_COMM_WORLD, &__myRank__); if (__myRank__==(R)) std::cout << DATA << std::endl;}

using namespace CommonLib;

namespace mfree
{

class MFreeException
//usage: throw new MFreeException(__LINE__, "Class::Method: A bloody unbelievable error occurred.");
{
public:
	int errorNum;
	char *errorText;

	MFreeException(const int _errorNum, const char *_errorText)
	{
		errorNum = _errorNum;
		errorText = new char[strlen(_errorText)+1];
		strcpy(errorText, _errorText);
	}
	
	virtual ~MFreeException()
	{
		delete []errorText;
	}
};

template <class T, int reallocSpace>
class C1DTable
//template class of a 1-D auto-allocate table
//only for simple types T that need not be explicitly constructed
//key method: operator[] returns reference to an element
//accessing the first non-existing element automatically adds one element to the table
//reallocSpace = number of elements allocated when out of space
{
protected:
	static const int maxSpace = 200*1024*1024/sizeof(T); //upper size limit in num of elements
	int allocedSpace, numElems; //space currently allocated, elements actually contained
	T* data; //element store (dynamic array)

	T& access(const int index)
	//private function used to access the value of an element in position index
	//if next non-present element (ie, index == numElems) is demanded, it is added to table
	//(and realloc is done if necessary)
	{
		assert(index >= 0 && index <= numElems);
		if (index == numElems)
		{
			numElems++;
			if (numElems > allocedSpace)
				reallocTable(allocedSpace + reallocSpace);
		}
		return data[index];
	}

	T& accessConst(const int index) const //not for left-side of assignment; is a const function
	//private function used to access the value of an element in position index
	//only existing elements may be accessed
	//not for left-side of assignment; is a const function
	{
		assert(index >= 0 && index < numElems);
		return data[index];
	}

public:
	C1DTable()
	{
		data = NULL;
		allocedSpace = 0;
		numElems = 0;
		reallocTable(reallocSpace);
	}

	C1DTable(const C1DTable &other)
	{
		data = NULL;
		allocedSpace = 0;
		numElems = 0;
		*this = other;
	}

	C1DTable(const int spaceToAlloc)
	{
		data = NULL;
		allocedSpace = 0;
		numElems = 0;
		assert(spaceToAlloc >= 0 && spaceToAlloc <= maxSpace);
		if (spaceToAlloc > 0)
			reallocTable(spaceToAlloc);
	}

	virtual ~C1DTable()
	{
		if (data)
		{
			free(data);
			data = NULL;
		}
	}

	void emptyTable()
	{
		numElems = 0;
	}

	void reallocTable(const int newSpace)
	//changes the allocated space for dynamic array
	//NOTE: new space can be larger or smaller (or the same) than old
	{
		assert(newSpace > 0 && newSpace < maxSpace && newSpace >= numElems);
		data = (T*)realloc(data, newSpace*sizeof(T));
		if (!data)
			throw new MFreeException(__LINE__, "CTable::reallocTable: out of memory");
		allocedSpace = newSpace;
	}

	int getNumElems() const
	{
		return numElems;
	}

	T *getRawPointer()
	//returns pointer to array of elements
	//USE WITH CAUTION
	{
		return data;
	}

	T& operator[](const int index) const
	//returns a const reference to element: this is called by default
	{
		return accessConst(index);
	}

	T& operator[](const int index)
	//returns a non-const reference to element: this is called automatically when used eg. on the left of =
	{
		return access(index);
	}

	void deleteElems(int firstToDelete, int lastToDelete)
	//deletes a range of elements
	//firstToDelete < 0 means "delete from beginning"
	//lastToDelete >= numElems means "delete to end"
	{
		assert(firstToDelete <= lastToDelete);
		firstToDelete = CommonLib::mmax(firstToDelete, 0);
		lastToDelete = CommonLib::mmin(lastToDelete, numElems-1);
		for (int i = lastToDelete+1; i < numElems; i++)
			access(i-1-lastToDelete+firstToDelete) = access(i);
		numElems -= lastToDelete-firstToDelete+1;
		reallocTable(numElems + reallocSpace);
	}

	void insertElem(int insertBefore, const T &newElem)
	//inserts element into the table
	{
		assert(insertBefore >= 0 && insertBefore <= numElems);
		access(numElems);
		for (int i = numElems-1; i > insertBefore; i--)
			access(i) = access(i-1);
		access(insertBefore) = newElem;
	}

	void appendElem(const T &newElem)
	//appends element to the end of the table
	{
		assert(numElems < maxSpace);
		access(numElems) = newElem;
	}

	void writeToFile(FILE *file) const
	//writes all the elements to file in raw binary format
	//NOTE: files may not be compatible between different machines, eg. big endian and little endian
	{
		assert(file);
		if (fwrite(data, sizeof(T), numElems, file) != (size_t)numElems)
			throw new MFreeException(__LINE__, "C1DTable::writeToFile: writing data failed");
		if (!fprintf(file, "\n"))
			throw new MFreeException(__LINE__, "C1DTable::writeToFile: writing terminating newline failed");
	}

	void readFromFile(FILE *file)
	//reads numElems elements from raw binary file
	//NOTE: numElems should be set before readFromFile() is called
	{
		assert(file);
		reallocTable(mmax(numElems, 1));
		if ((int)fread(data, sizeof(T), numElems, file) != numElems)
			throw new MFreeException(__LINE__, "C1DTable::readFromFile: reading data failed");
		if (fgetc(file) != '\n')
			throw new MFreeException(__LINE__, "C1DTable::readFromFile: terminating newline not found");
	}

	C1DTable<T, reallocSpace> &operator=(const C1DTable<T, reallocSpace> &other)
	{
		numElems = 0;
		reallocTable(other.allocedSpace);
		for (int i = 0; i < other.numElems; i++)
			access(i) = other.accessConst(i);
		return *this;
	}

}; //end class CTable

struct CPoint;

struct CRect //rectangle spanning from xMin to xMax in x direction, similar in y
{
	double xMin, xMax;
	double yMin, yMax;

	CRect()
	//default constructor creates an invalid object
	{
		xMin = CommonLib::NaN;
		xMax = CommonLib::NaN;
		yMin = CommonLib::NaN;
		yMax = CommonLib::NaN;
	}

	CRect(const CRect &other)
	{
		xMin = other.xMin;
		xMax = other.xMax;
		yMin = other.yMin;
		yMax = other.yMax;
	}

	CRect(double _xMin, double _xMax, double _yMin, double _yMax)
	{
		xMin = _xMin;
		xMax = _xMax;
		yMin = _yMin;
		yMax = _yMax;
	}
	
	bool overlaps(const CRect &other) const;	//true iff two rectangles overlap

    CRect operator*(const CRect &other) const
	//returns intersection of two rectangles;
	//if they do not overlap, returns CRect(NaN,NaN,NaN,NaN)
	{
		CRect res;
		if (overlaps(other))
		{
			res.xMin = CommonLib::mmax(xMin, other.xMin);
			res.xMax = CommonLib::mmin(xMax, other.xMax);
			res.yMin = CommonLib::mmax(yMin, other.yMin);
			res.yMax = CommonLib::mmin(yMax, other.yMax);
		}
		return res;
	}

	CRect operator+(const CRect &other) const
	//returns the smallest rectangle containing given two rectangles
	{
		CRect res;
		res.xMin = CommonLib::mmin(xMin, other.xMin);
		res.xMax = CommonLib::mmax(xMax, other.xMax);
		res.yMin = CommonLib::mmin(yMin, other.yMin);
		res.yMax = CommonLib::mmax(yMax, other.yMax);
		return res;
	}

	CRect operator+(const CPoint &p) const; //smallest rectangle containing given rectangle and point

	CRect operator+(double boundary) const
	//rectangle enlarged by boundary on each side
	{
		CRect res;
		res.xMin = xMin-boundary;
		res.xMax = xMax+boundary;
		res.yMin = yMin-boundary;
		res.yMax = yMax+boundary;
		return res;
	}

	bool operator==(const CRect &other) const
	{
		return xMin == other.xMin && xMax == other.xMax
				&& yMin == other.yMin && yMax == other.yMax;
	}

	double getWidth() const
	{
		return xMax-xMin;
	}

	double getHeight() const
	{
		return yMax-yMin;
	}
}; //end class CRect

struct CPoint //2D point (x,y)
{
	double x, y;

	CPoint()
	//default constructor creates an invalid object
	{
		x = CommonLib::NaN;
		y = CommonLib::NaN;
	}

	CPoint(double _x, double _y)
	{
		x = _x;
		y = _y;
	}

	CPoint(const CPoint &other)
	{
		x = other.x;
		y = other.y;
	}

	CPoint &operator=(const CPoint &other)
	{
		x = other.x;
		y = other.y;
		return *this;
	}

	CPoint operator+(const CPoint &other) const
	//sum of 2 vectors represented by given points
	{
		return CPoint(x+other.x, y+other.y);
	}

	CPoint operator-(const CPoint &other) const
	{
		return CPoint(x-other.x, y-other.y);
	}

	bool operator==(const CPoint &other) const
	{
		return x == other.x && y == other.y;
	}

	bool operator!=(const CPoint &other) const
	{
		return !(*this == other);
	}

	CPoint operator*(double scalar) const
	{
		return CPoint(x*scalar, y*scalar);
	}

	CPoint operator/(double scalar) const
	{
		assert(scalar != 0);
		return CPoint(x/scalar, y/scalar);
	}

	double length() const
	//length of vector represented by a point
	{
		return sqrt(x*x + y*y);
	}

	bool isInRect(const CRect &rect) const
	//true iff point lies in rectangle or on its boundary
	{
		return (CommonLib::in(x, rect.xMin, rect.xMax) && CommonLib::in(y, rect.yMin, rect.yMax));
	}

	double distToClosest(const CRect &rect) const
	//distance from point to closest point of the rectangle
	{
		if (x > rect.xMax)
			if (y > rect.yMax)
				//closest point is upper-right corner
				return (*this-CPoint(rect.xMax, rect.yMax)).length();
			else if (y < rect.yMin)
				//closest point is lower-right corner
				return (*this-CPoint(rect.xMax, rect.yMin)).length();
			else
				//closest point is on right border
				return x-rect.xMax;
		else if (x < rect.xMin)
			if (y > rect.yMax)
				//closest point is upper-left corner
				return (*this-CPoint(rect.xMin, rect.yMax)).length();
			else if (y < rect.yMin)
				//closest point is lower-left corner
				return (*this-CPoint(rect.xMin, rect.yMin)).length();
			else
				//closest point is on left border
				return rect.xMin-x;
		else
			if (y > rect.yMax)
				//closest point is on upper border
				return y-rect.yMax;
			else if (y < rect.yMin)
				//closest point is on lower border
				return rect.yMin-y;
			else 
				return 0; //point in rect
	}

	void writeToFile(FILE *f) const
	{
		assert(f);
		fprintf(f, "%lf %lf", x, y);
	}

}; //end class CPoint

inline CRect CRect::operator+(const CPoint &p) const
//smallest rectangle containing given rectangle and point
{
	CRect res;
	res.xMin = CommonLib::mmin(xMin, p.x);
	res.xMax = CommonLib::mmax(xMax, p.x);
	res.yMin = CommonLib::mmin(yMin, p.y);
	res.yMax = CommonLib::mmax(yMax, p.y);
	return res;
}

inline bool CRect::overlaps(const CRect &other) const
//true iff rects overlay
{
	//two rectangles overlay iff one corner of first rectangle lies in second rectangle
	return CPoint(xMin,yMin).isInRect(other) || CPoint(xMin,yMax).isInRect(other)
			|| CPoint(xMax,yMin).isInRect(other) || CPoint(xMax,yMax).isInRect(other);
}

inline std::ostream &operator<<(std::ostream &str, const mfree::CRect &rect)
//writes rectangle in format (xMin-xMax)x(yMin-yMax)
{
	str << "(" << rect.xMin << "-" << rect.xMax << ")x(" << rect.yMin << "-" << rect.yMax << ")";
	return str;
}

struct MFreeOptions //options & parameters defining everything for the MLPG1 method
{
	int suggestedNumNodes;	//actual number of nodes will be similar, but not always equal to suggestedNumNodes (because of limitations of our node generation algorithm)
	double irregularity;	//node position irregularity
	int randSeed;			//seed for random number generator; using the same seed ensures reproducibility of results

	char *nodesFile;		//file to read nodes from

	double rHole;			//hole diameter; hole boundary must not intersect unit square boundary
	CPoint xyHole;			//hole position

	double c;				//diffusivity constant; > 0

	int MLSdegree;			//monomial degree; 1 or 2
	int quadDegree;			//quadrature degree; 1, 3, 5, ... 11

	double timeStep;		

	double alphaSupport, alphaQuad; //dimensionless support and quadrature domain sizes
							//(if alphaSupport is NaN, constant number of support nodes is used)
	int nI;					//number of nodes in support domain (if <=0, alphaSupport/avgDistMesh is used)
	double betaQuad;		//ratio of quad domain size/sup domain size, for the case where constant suppord node count is used

	bool preFindQuadSupport; //iff true, isciKvad version of algorithm is used
	int verboseLevel;		//amount of data printed during execution
	int numRepetitions;		//number of times everything will be executed; use >1 for small node numbers to ensure accurate time measurements

	enum DataDistrType //data disribution types for parallel execution
	{
		DataDistrHierarchical,	//each subtree goes to one sub-hypercube
								//each proc gets its domain + border
								//avgDistMesh is evaluated before distribution
		DataDistr1D,			//1-D distribution on y
								//each proc gets its domain + estimated border
								//avgDistMesh is evaluated after distribution
		DataDistrUnstructured	//unstructured; currently implemented as 1-D distribution,
								//but each proc gets all nodes
								//avgDistMesh is evaluated after distribution
	};
	DataDistrType dataDistrType; //selected type of data distribution

	MFreeOptions()
	//default options
	{
		suggestedNumNodes = 49;
		randSeed = 17;
		nodesFile = NULL;
		rHole = 0;
		xyHole = CPoint(0.6, 0.35);
		c = 0.001;
		timeStep = 5/40.0;
		irregularity = 0.3;
		alphaSupport = NaN;
		alphaQuad = 0.8;
		nI = 13;
		betaQuad = 0.7;
		MLSdegree = 2;
		quadDegree = 5;
		preFindQuadSupport = false;
		verboseLevel = 0;
		numRepetitions = 1;
		dataDistrType = DataDistrHierarchical;
	}

	bool constSupNodeCount() const
	{
		return isnan(alphaSupport);
	}

	bool validate() const 
	//true iff none of the option values is nonsense
	{
		return CommonLib::in(suggestedNumNodes, 10, 100000)
			&& CommonLib::in(rHole, 0.0, 0.5)
			&& CommonLib::in(xyHole.x, rHole, 1-rHole)
			&& CommonLib::in(xyHole.y, rHole, 1-rHole)
			&& CommonLib::in(irregularity, 0.0, 1.0)
			&& CommonLib::in(MLSdegree, 1, 2)
			&& CommonLib::in(quadDegree, 1, 11)
			&& quadDegree % 2 == 1
			&& (CommonLib::in(alphaSupport, 0.4, 10.0) || nI >= 3)
			&& (isnan(alphaSupport) || nI <= 0)
			&& CommonLib::in(alphaQuad, 0.01, 10.0)
			&& CommonLib::in(betaQuad, 0.01, 10.0)
			&& CommonLib::in(verboseLevel, 0, 100)
			&& (nI <= 0 || !preFindQuadSupport);
	}

	void writeToMFile(FILE *f, char *varName, int numNodes) const
	//export to structure varName.mfreeParams using Matlab script format
	{
		fprintf(f, "clear mfreeParams; global mfreeParams; \n");
		fprintf(f, "mfreeParams.numNodes = %d;\n", numNodes);
		fprintf(f, "mfreeParams.randSeed = %d;\n", randSeed);
		fprintf(f, "domain.rHole = %.15f;\n", rHole);
		fprintf(f, "domain.xHole = %.15f;\n", xyHole.x);
		fprintf(f, "domain.yHole = %.15f;\n", xyHole.y);
		fprintf(f, "mfreeParams.C = %.15f;\n", c);
		fprintf(f, "mfreeParams.irregularity = %.15f;\n", irregularity);
		fprintf(f, "mfreeParams.baseFtype = 'mls';\n");
		fprintf(f, "mfreeParams.mlsDegree = %d;\n", MLSdegree);
		fprintf(f, "mfreeParams.quadDegree = %d;\n", quadDegree);
		fprintf(f, "mfreeParams.timeStep = %.15f;\n", timeStep);
		if (constSupNodeCount()) {
			fprintf(f, "mfreeParams.nI = %d;\n", nI);
			fprintf(f, "mfreeParams.betaQuad = %.15f;\n", betaQuad);
		} else {
			fprintf(f, "mfreeParams.alphaInterp = %.15f;\n", alphaSupport);
			fprintf(f, "mfreeParams.alphaQuad = %.15f;\n", alphaQuad);
		}
		fprintf(f, "mfreeParams.preFindQuadSupport = %d;\n", (preFindQuadSupport ? 1 : 0));
		fprintf(f, "mfreeParams.numRepetitions = %d;\n", numRepetitions);
		fprintf(f, "mfreeParams.dataDistrType = %d;\n",
					(dataDistrType == DataDistrHierarchical ? 0 : (dataDistrType == DataDistr1D ? 1 : 2)));
	}
};

}; //end namespace mfree

#endif
