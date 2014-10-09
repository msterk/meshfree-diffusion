#ifndef __LIBDATA_AUX_H_INCLUDED__
#define __LIBDATA_AUX_H_INCLUDED__

#include <mpi.h>
#include <assert.h>
#include <stdlib.h>
#include <iostream>
#include <stdio.h>

#include <commonlib_aux.h>

namespace CommonLib {

class ILoadDistribution
//interface to 1D load distributions for data indexed 1..size
//the purpose is to distribute non-border elements (2..size-1) between procs somehow
//index 1 implicitly belongs to process owning index 2; similarly, size belongs to the same process as size-1
{
protected:
	int size;		//size of the distributed data
	int myRank;		//rank of the process on which this object resides
	int numProcs;	//total number of processes

public:
	ILoadDistribution(int _size, int _rank, int _numProcs)
		: size(_size)
	{
		assert(size > 0);
		assert(_numProcs > 0);
		assert(_rank >= 0 && _rank < _numProcs);
		numProcs = _numProcs;
		myRank = _rank;
	}

	virtual ~ILoadDistribution()
	{}

	int getNumProcs() const
	{
		return numProcs;
	}

	int getRank() const
	{
		return myRank;
	}

	virtual int getDomainStart(int rank) const = 0;
	//abstract function that should return the start of subdomain owned by process 'rank'
	//	i.e. process rank owns data from getDomainStart(rank) to getDomainEnd(rank)
	//should return size if the process owns no data

	int getDomainStart() const
	//start of subdomain of process this->myRank
	{
		return getDomainStart(myRank);
	}

	int getDomainEnd(int rank) const
	//end of subdomain
	{
		assert(rank >= 0 && rank < numProcs);
		return (rank == numProcs-1 ? size-1 : getDomainStart(rank+1)-1);
	}

	int getDomainEnd() const
	//end of subdomain of process this->myRank
	{
		return getDomainEnd(myRank);
	}

	bool isDomainEmpty(int rank) const
	{
		return getDomainEnd(rank) < getDomainStart(rank);
	}

	bool isDomainEmpty() const
	{
		return isDomainEmpty(myRank);
	}

	int getDomainSize(int rank) const
	{
		return isDomainEmpty(rank) ? 0 : getDomainEnd(rank)-getDomainStart(rank)+1;
	}

	int getDomainSize() const
	{
		return getDomainSize(myRank);
	}

	int getRightNeighbour(int rank) const
	//returns the rank of process that owns the subdomain following the subdomain of given process
	//	(in straightforward 1D distribution, this is always rank+1)
	//returns -1 if not found
	{
		assert(rank >= 0 && rank < numProcs);
		for (int i = rank+1; i < numProcs; i++)
			if (getDomainStart(i) == getDomainEnd(rank)+1 && !isDomainEmpty(i))
				return i;
		return -1;
	}

	int getRightNeighbour() const
	{
		return getRightNeighbour(myRank);
	}

	int getLeftNeighbour(int rank) const
	//returns the rank of process that owns the subdomain preceding the subdomain of given process
	//	(in straightforward 1D distribution, this is always rank-1)
	//returns -1 if not found
	{
		assert(rank >= 0 && rank < numProcs);
		for (int i = rank-1; i >= 0; i--)
			if (getDomainEnd(i) == getDomainStart(rank)-1 && !isDomainEmpty(i))
				return i;
		return -1;
	}

	int getLeftNeighbour() const
	{
		return getLeftNeighbour(myRank);
	}

	int getOwner(int index) const
	//returns the rank of processor owning the data at given index
	//returns -1 if not found
	{
		assert(index >= 2 && index < size);
		for (int proc = 0; proc < numProcs; proc++)
			if (!isDomainEmpty(proc) && getDomainStart(proc) <= index && getDomainEnd(proc) >= index)
				return proc;
		return -1;
	}

	bool isDistributed()
	//returns false IFF everything is on root
	{
		return !(getDomainSize(0) > 0 && getRightNeighbour(0) == -1);
	}

	virtual bool coarsenDist(ILoadDistribution &fineDist) = 0;
	//abstract function that should make *this a distribution of data of half the size of fineDist
	//	used eg in restriction operator in parallel multigrid method; each subdomain should also be
	//	the coarsened version of corresponding subdomain in fineDist
	//should return true if successful;
	//	that is, if no process has a non-empty subdomain in fineDist but an empty one in the coarse distribution

}; //end class ILoadDistribution


class CNondistributedLoadDistribution : public ILoadDistribution
//implementation of ILoadDistribution that just puts everything on process 0
{
public:
	CNondistributedLoadDistribution(int _size, int _rank, int _numProcs)
		: ILoadDistribution(_size, _rank, _numProcs)
	{}
	
	virtual int getDomainStart(int rank) const
	{
		assert(rank >= 0 && rank < numProcs);
		return (rank == 0 ? 2 : size);
	}

	virtual bool isDistributed()
	{
		return false;
	}

	virtual bool coarsenDist(ILoadDistribution &fineDist)
	{
		return true;
	}
}; //end class CNondistributedLoadDistribution


class CUniformLoadDistribution : public ILoadDistribution
//implementation that distributes data as close to uniformly as possible
{
public:
	CUniformLoadDistribution(int _size, int _rank, int _numProcs)
		: ILoadDistribution(_size, _rank, _numProcs)
	{}

	virtual int getDomainStart(int rank) const
	{
		assert(rank >= 0 && rank < numProcs);
		int num = size-2;
		int firstFew = num % numProcs; //these will have one plane more
		int zPerProc = num / numProcs; //the rest will have this many planes
		if (rank <= firstFew)
			return (zPerProc+1) * rank + 2;
		else
			return (zPerProc+1) * firstFew + zPerProc * (rank-firstFew) + 2;
	}
	virtual bool coarsenDist(ILoadDistribution &fineDist)
	{
		assert("Cannot coarsen uniform distribution" == NULL);
		return false;
	}
}; //end class CUniformLoadDistribution


class CCustomLoadDistribution : public CUniformLoadDistribution
//implementation that can distribute data in any way
//on construction, uniform distribution is created
{
	int *domainStarts;		//start of each individual subdomain
	bool broadcastFinished;	//true iff array domainStarts has been successfully broadcast and is thus valid on all processes

public:
	CCustomLoadDistribution(int _size, int _rank, int _numProcs)
		: CUniformLoadDistribution(_size, _rank, _numProcs)
	//creates uniform distribution
	{
		domainStarts = new int[numProcs];
		assert(domainStarts);
		for (int i = 0; i < numProcs; i++)
			domainStarts[i] = CUniformLoadDistribution::getDomainStart(i);
		broadcastFinished = true;
	}

	virtual ~CCustomLoadDistribution()
	{
		if (domainStarts)
			delete domainStarts;
	}

	virtual int getDomainStart(int rank) const
	{
		assert(rank >= 0 && rank < numProcs);
		assert(broadcastFinished);
		return domainStarts[rank];
	}

	void setDomainStart(int rank, int start)
	//changes the distribution and marks the need to broadcast the new distribution
	{
		assert(rank >= 0 && rank < numProcs);
		assert(start > 1 && start <= size);
		domainStarts[rank] = start;
		broadcastFinished = false;
	}

	void broadcastDistribution()
	//broadcasts the distribution to all processes so that each has full knowledge of distribution
	{
		MPI_Bcast(domainStarts, numProcs, MPI_INT, 0, MPI_COMM_WORLD);
		broadcastFinished = true;
	}

	bool coarsenDist(ILoadDistribution &fineDist, int minDomainSize)
	//coarsen the distribution according to requirements of parallel multigrid
	{
		int coarseSizesOK;
		if (myRank == 0) //find a coarse distribution on root
		{
			//assign points to processor where the corresponding fine points were
			for (int proc = 0; proc < numProcs; proc++)
				setDomainStart(proc, (fineDist.getDomainStart(proc)+2)/2);
			//return false if any coarse domain is too small
			coarseSizesOK = 1;
			for (int proc = 0; proc < numProcs; proc++)
				if ((proc == numProcs-1 ? size : domainStarts[proc+1]) - domainStarts[proc] < minDomainSize
					&& fineDist.getDomainSize(proc) > 0)
				{
					coarseSizesOK = 0;
					break;
				}
		}
		broadcastDistribution(); //broadcast it to all
		commBroadcast(coarseSizesOK);
		return (coarseSizesOK ? true : false);
	}
}; //end class CCustomLoadDistribution

}; //end namespace CommonLib

#endif


