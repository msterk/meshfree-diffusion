#ifndef __COMMONLIB_AUX_H__
#define __COMMONLIB_AUX_H__

#ifndef WIN32
	#define stricmp strcasecmp	//in Visual C, the function is named strcasecmp
								//we want to always use the name stricmp
#endif

#include <assert.h>
#include <mpi.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <float.h>

namespace CommonLib {

    const double pi = 3.1415926535897932384626433832795; 
	const double eps = 1e-15;	//don't really know the purpose of that, so "Don't touch it!"

//definitions to unify usage of Not-a-Number in Visual C++ and linux/gcc
	//NaN is available as __builtin_nan("") in gcc, but not available in VC; thus we create NaN for both cases
	//similarly, there is isnan() in gcc but _isnan() in VC; thus we create isNaN for both cases
	#ifdef WIN32
		const unsigned long ___nanConstruction[2]={0xffffffff, 0x7fffffff};
		const double NaN = *( double* )___nanConstruction;

		inline bool isNaN(const double x)
		{
			return (_isnan(x) ? true : false);
		}
	#else
		const double NaN = __builtin_nan("");
		inline bool isNaN(const double x)
		{
			return isnan(x);
		}
	#endif

//simple function templates that need no explanation
	template <class T> void swap(T &x, T &y)
	{
		T temp;
		temp = x;
		x = y;
		y = temp;
	}
	template <class T> inline T randFromRange(T min, T max) 
	{
		return (T)rand()/RAND_MAX*(max-min) + min;
	}
	template <class T> inline T mmin(T a, T b)
	{
		return a < b ? a : b;
	}
	template <class T> inline T mmin(T a, T b, T c)
	{
		return mmin(a,mmin(b,c));
	}
	template <class T> inline T mmin(T a, T b, T c, T d)
	{
		return mmin(mmin(a,b),mmin(c,d));
	}
	template <class T> inline T mmax(T a, T b)
	{
		return a > b ? a : b;
	}
	template <class T> inline T mmax(T a, T b, T c)
	{
		return mmax(a,mmax(b,c));
	}
	template <class T> inline T mmax(T a, T b, T c, T d)
	{
		return mmax(mmax(a,b),mmax(c,d));
	}
	template <class T> inline T sqr(T x)
	{
		return x*x;
	}
	template <class T> inline T avg(T x, T y) 
	{ 
		return (x+y)/2;
	}
	template <class T> inline T sign(T x) 
	{
		return (x > 0) ? 1 : ((x < 0) ? -1 : 0); 
	}
	inline int round(double x)
	{
		return int(x+0.5);
	}
	inline bool in(int x, int lowerBound, int upperBound) //true iff x in [lb, ub]
	{
		return (x >= lowerBound) && (x <= upperBound);
	}
	inline bool in(double x, double lowerBound, double upperBound) //true iff x in [lb, ub]
	{
		return (x >= lowerBound) && (x <= upperBound);
	}

//the graphics functions only work in windows
	#ifdef WIN32
		extern "C" void CreateBMPFile(char *pszFile, void *hBMP, void *hDC);
		//given a handle hBMP of a bitmap that is currently selected in the display context hDC,
		//creates a .bmp file named pszFile containing the bitmap

		unsigned long HLSToRGB(int hue, int lum, int sat);
		//converts a color given as (h,l,s) to RGB(r,g,b)

	#endif //WIN32

//combinatorics
	class CCombination
	//represents a combination of numEl different elements taken from interval [minEl..maxEl]
	{
		int numEl,			//number of elements comprising a combination
			minEl, maxEl;	//interval that element values are taken from
		int *el;			//array of currently selected elements, in ascending order
		int combNumber;		//index of the combination
							//combNumber = 0 stands for combination < minEl, minEl+1, ..minEl+(numEl-1) >
							//combNumber = binomial_symbol( (maxEl-minEl+1) , numEl )
							//	stands for combination < maxEl-(numEl-1), ... , maxEl-1, maxEl >

	public:
		CCombination(int _numEl, int _minEl, int _maxEl)
		//creates the combination with index 0
		{
			assert(_numEl <= _maxEl - _minEl +1);
			numEl = _numEl;
			minEl = _minEl;
			maxEl = _maxEl;
			el = new int[numEl];
			for (int i = 0; i < numEl; i++)
				el[i] = i+minEl;
			combNumber = 0;
		}

		~CCombination()
		{
			delete []el;
		}

		int &operator[](int i)
		//returns reference to i+th element
		{
			return el[i];
		}

		bool next()
		//if current combination is the last one possible, sets combNumber to -1 and returns false
		//otherwise, increments combNumber by 1, updates the combination accordingly and returns true
		{
			combNumber++;
			int hole; //index of last element whose value is smaller by 2 or more than the value of next element
					  //enlarging elements hole..numEl-1 by 1 gives the next combination

			//if value of last element is not maxEl, just enlarge it & done
			if (el[numEl-1] < maxEl)
			{
				el[numEl-1]++;
				return true;
			}
			//otherwise, find hole
            for (hole = numEl-2; ; hole--)
			{
				if (hole < 0)
				{
					combNumber = -1;
					return false;
				}
				if (el[hole]+1 < el[hole+1])
					break;
			}
			//enlarge el[hole..numEl-1] by one
			el[hole]++;
			for (int i = hole+1; i < numEl; i++)
				el[i] = el[i-1]+1;
			return true;
		}

		int getCombNumber()
		{
			return combNumber;
		}

		int countCombsApprox()
		//returns number of possible combinations, disregarding overflows & rounding errors
		//calculation is done in double; otherwise, overflow would happen even for small numEl
		{
			double c = 1;
			for (int i = 1; i <= numEl; i++)
				c *= (double)(maxEl-minEl+2-i)/i;
			return (int)c;
		}
	}; //end class CCombination

//simplified wrappers for some MPI functions
//all of these should be run on all processes that are part of MPI_COMM_WORLD
	inline void commBroadcast(double &x, int count = 1)
	//broadcasts 1 or 'count' doubles from process 0 to all processes
	{
		MPI_Bcast(&x, count, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	}

	inline void commBroadcast(int &x, int count = 1)
	//broadcasts 1 or 'count' ints from process 0 to all processes
	{
		MPI_Bcast(&x, count, MPI_INT, 0, MPI_COMM_WORLD);
	}

	template <typename T>
	inline void commBroadcast(T &x, int count = 1)
	//broadcasts 1 or 'count' Ts (sent as bytes) from process 0 to all processes
	{
		//std::cout << "T size " << sizeof(T) << ", count " << count << std::endl;
		MPI_Bcast(&x, count*sizeof(T), MPI_BYTE, 0, MPI_COMM_WORLD);
	}

	inline void commAllReduce(double &x, MPI_Op op)
	//Allreduces a double using op
	{
		double temp = x;
		MPI_Allreduce(&temp, &x, 1, MPI_DOUBLE, op, MPI_COMM_WORLD);
	}

	inline void commAllReduce(int &x, MPI_Op op)
	//Allreduces an int using op
	{
		int temp = x;
		MPI_Allreduce(&temp, &x, 1, MPI_INT, op, MPI_COMM_WORLD);
	}

	template <typename T>
	inline void commDbgOutput(const T &x, char *prolog = NULL)
	//on process 0, prints prolog and local values of x from all processes
	//on other processes, just sends local value to 0
	{
		T allX[100];
		int myRank, numProcs;
		MPI_Comm_size(MPI_COMM_WORLD,&numProcs);
		MPI_Comm_rank(MPI_COMM_WORLD,&myRank);
		MPI_Gather((void*)&x, sizeof(T), MPI_BYTE, allX, sizeof(T), MPI_BYTE, 0, MPI_COMM_WORLD);
		if (myRank == 0)
		{
			if (prolog)
				std::cout << prolog << std::endl;
			for (int i = 0; i < numProcs; i++)
				std::cout << "Proc " << i << ": " << allX[i] << std::endl;
			std::cout << std::flush;
		}
	}
}; //end namespace CommonLib

#endif

