#ifndef __TABLE_H_INCLUDED__
#define __TABLE_H_INCLUDED__

#include "loaddistribution.h"
#include <iostream>

//uncomment the first to disable, second to enable tracing messages
//#define LOGMESSAGE(operation, other, numBytes) ;
#define LOGMESSAGE(operation, other, numBytes) if (loadDist->getRank() == 0) printf("Process %3d %13s %3d msg of %8d bytes.\n", loadDist->getRank(), operation, other, numBytes)

namespace CommonLib
{

template <class T>
class CTable
//3D table with 1D distribution
//(sliced in xy plane, so values <x,y,z> and <x,y,z+1> potentially reside on different processes
//
//each process also stores the neighbouring slice belonging to neighbour process
//	to facilitate finite difference computations
//
//the class was written to store boundary values in addition to values <1,1,1>..<sizeX,sizeY,sizeZ>,
//	i.e. all values <0,0,0>..<sizeX+1,sizeY+1,sizeZ+1> were to be stored; but this was never tested
{
	T *data;							//stores the data
	ILoadDistribution *loadDist;		//description of load distribution
	int domainStart, domainEnd;			//subdomain boundaries of current process
	MPI_Datatype checkerBoardSubsetType;//MPI type that facilitates sending just the values of black or red checkerboard fields
	int checkerBoardSubsetTypeSize;		//size of checkerBoardSubsetType

	int getIndex(int x, int y, int z)
	//returns index at which value <x,y,z> is stored in the array *data
	{
		assert(data != NULL);
		assert(x >= 0 && x <= intSizeX+1);
		assert(y >= 0 && y <= sizeY+1);
		assert(z >= domainStart-1 && z <= domainEnd+1);
		return ((z - domainStart + 1)*(sizeY + 2) + y) * (intSizeX + 2) + x;
	}

	T &access(int x, int y, int z)
	//returns reference to a value
	{
		return data[ getIndex(x, y, z) ];
	}

public:
	int sizeX, sizeY, sizeZ;//size of table in each dimension
	int	intSizeX;			//x-size used internally for data storage
				//intSizeX is either sizeX (natural storing) or sizeX+1 (additional empty memory) so that
				//intSizeX is always odd for easy checkerboard-pattern communication

	CTable()
	//creates invalid object
	{
		data = NULL;
		loadDist = NULL;
		sizeX = sizeY = sizeZ = -1;
	}

	void allocate(int _sizeX, int _sizeY, int _sizeZ, ILoadDistribution *_loadDist)
	//allocates space for table of given size and distribution
	//creates checkerBoardSubsetType
	{
		assert(_sizeX > 0 && _sizeY > 0 && _sizeZ > 0);
		assert(_loadDist);
		sizeX = _sizeX;
		intSizeX = (sizeX % 2 == 0 ? sizeX + 1 : sizeX);
		sizeY = _sizeY;
		sizeZ = _sizeZ;
		loadDist = _loadDist;
		domainStart = loadDist->getDomainStart();
		domainEnd = loadDist->getDomainEnd();

		data = new T[(intSizeX+2)*(sizeY+2)*(domainEnd - domainStart +3)];
		assert(data != NULL);

		if (loadDist->getNumProcs() > 1)
		{
			//checkerBoardSubsetType is a vector with stride 2 and with elements of sizeof(T) MPI_BYTEs each
			MPI_Type_vector( int((getIndex(intSizeX+1, sizeY+1, domainStart) - getIndex(0, 0, domainStart) + 1) /2),
								sizeof(T), 2*sizeof(T), MPI_BYTE, &checkerBoardSubsetType);
			MPI_Type_commit(&checkerBoardSubsetType);
			MPI_Type_size(checkerBoardSubsetType, &checkerBoardSubsetTypeSize);
		}
	}

	CTable(int _sizeX, int _sizeY, int _sizeZ, ILoadDistribution *_loadDist)
	{
		allocate(_sizeX, _sizeY, _sizeZ, _loadDist);
	}

	~CTable() 
	{
		if (data != NULL)
			delete []data;
	}

	void resetValues()
	//sets all values to 0
	{
		for (int z = domainStart-1; z <= domainEnd+1; z++)
			for (int y = 0; y <= sizeY+1; y++)
				for (int x = 0; x <= intSizeX+1; x++)
					access(x, y, z) = 0;
	}

	void exchangeValues(CTable<T> &t2)
	//swaps values between this and CTable t2
	{
		assert(sizeX ==t2.sizeX && sizeY == t2.sizeY && sizeZ == t2.sizeZ);
		T *temp;
		temp = data;
		data = t2.data;
		t2.data = temp;
	}

	T &operator()(int x, int y, int z)
	{
		return access(x, y, z);
	}

	enum CommSubset {
		subsetBlack,
		subsetRed,
		subsetAll
	};

	void commSend(int zFrom, int zTo, int dest)
	//send all data in slices [zFrom..zTo] to process dest
	{
		assert(zFrom >= domainStart-1);
		assert(zTo >= zFrom && zTo <= domainEnd+1);
		MPI_Send(&data[ getIndex(0, 0, zFrom) ],
			sizeof(T)*( getIndex(intSizeX+1, sizeY+1, zTo) - getIndex(0, 0, zFrom) + 1 ),
			MPI_BYTE, dest, zTo+zFrom, MPI_COMM_WORLD);
		LOGMESSAGE("sending to", dest, sizeof(T)*( getIndex(intSizeX+1, sizeY+1, zTo) - getIndex(0, 0, zFrom) + 1) );
	}

	void commISend(int zFrom, int zTo, int dest, MPI_Request &request, CommSubset subsetToSend = subsetAll)
	//Isend either all or a subset of all data in slices [zFrom..zTo] to process dest
	{
		assert(zFrom >= domainStart-1);
		assert(zTo >= zFrom && zTo <= domainEnd+1);

		switch (subsetToSend)
		{
		case subsetAll:
			MPI_Isend(&data[ getIndex(0, 0, zFrom) ],
				sizeof(T)*( getIndex(intSizeX+1, sizeY+1, zTo) - getIndex(0, 0, zFrom) + 1 ),
				MPI_BYTE, dest, zTo+zFrom, MPI_COMM_WORLD, &request);
			LOGMESSAGE("Isending to", dest, sizeof(T)*( getIndex(intSizeX+1, sizeY+1, zTo) - getIndex(0, 0, zFrom) + 1 ) );
			break;
		case subsetBlack:
			assert(zFrom == zTo);
			//black subset is checkerBoardSubsetType starting at x = zFrom % 2
			MPI_Isend(&data[ getIndex(zFrom % 2, 0, zFrom) ], 1, checkerBoardSubsetType, dest, zTo+zFrom, MPI_COMM_WORLD, &request);
			LOGMESSAGE("Isending to", dest, checkerBoardSubsetTypeSize );
			break;
		case subsetRed:
			assert(zFrom == zTo);
			//red subset is checkerBoardSubsetType starting at x = (zFrom+1) % 2
			MPI_Isend(&data[ getIndex((zFrom+1) % 2, 0, zFrom) ], 1, checkerBoardSubsetType, dest, zTo+zFrom, MPI_COMM_WORLD, &request);
			LOGMESSAGE("Isending to", dest, checkerBoardSubsetTypeSize );
			break;
		default:
			assert("Not yet implemented" == NULL);
			break;
		}
	}

	void commRecv(int zFrom, int zTo, int source)
	//recv all data in slices [zFrom..zTo] from process source
	{
		assert(zFrom >= domainStart-1);
		assert(zTo >= zFrom && zTo <= domainEnd+1);
		MPI_Status status;
		MPI_Recv(&data[ getIndex(0, 0, zFrom) ],
			sizeof(T)*( getIndex(intSizeX+1, sizeY+1, zTo) - getIndex(0, 0, zFrom) + 1 ),
			MPI_BYTE, source, zTo+zFrom, MPI_COMM_WORLD, &status);
		LOGMESSAGE("recving from", source, sizeof(T)*( getIndex(intSizeX+1, sizeY+1, zTo) - getIndex(0, 0, zFrom) + 1 ) );
	}

	void commIRecv(int zFrom, int zTo, int source, MPI_Request &request, CommSubset subsetToRecv = subsetAll)
	//Isrecv either all or a subset of all data in slices [zFrom..zTo] from process source
	{
		assert(zFrom >= domainStart-1);
		assert(zTo >= zFrom && zTo <= domainEnd+1);

		switch (subsetToRecv)
		{
		case subsetAll:
			MPI_Irecv(&data[ getIndex(0, 0, zFrom) ],
				sizeof(T)*( getIndex(intSizeX+1, sizeY+1, zTo) - getIndex(0, 0, zFrom) + 1 ),
				MPI_BYTE, source, zTo+zFrom, MPI_COMM_WORLD, &request);
			LOGMESSAGE("Irecving from", source, sizeof(T)*( getIndex(intSizeX+1, sizeY+1, zTo) - getIndex(0, 0, zFrom) + 1 ) );
			break;
		case subsetBlack:
			assert(zFrom == zTo);
			//black subset is checkerBoardSubsetType starting at x = zFrom % 2
			MPI_Irecv(&data[ getIndex(zFrom % 2, 0, zFrom) ], 1, checkerBoardSubsetType, source, zTo+zFrom, MPI_COMM_WORLD, &request);
			LOGMESSAGE("Irecving from", source, checkerBoardSubsetTypeSize );
			break;
		case subsetRed:
			assert(zFrom == zTo);
			//red subset is checkerBoardSubsetType starting at x = (zFrom+1) % 2
			MPI_Irecv(&data[ getIndex((zFrom + 1) % 2, 0, zFrom) ], 1, checkerBoardSubsetType, source, zTo+zFrom, MPI_COMM_WORLD, &request);
			LOGMESSAGE("Irecving from", source, checkerBoardSubsetTypeSize );
			break;
		default:
			assert("Not yet implemented" == NULL);
			break;
		}
	}

	void commExchangeBorders(CommSubset subsetToExchange = subsetAll)
	//send both subdomain borders to respective neighbour processes and
	//	receive their respective borders
	//does either complete borders or just a red or black subset of them
	{
		if (loadDist->isDomainEmpty())
			return;
		MPI_Request request[4];
		MPI_Status status[4];
		int commNum = 0;

		if (loadDist->getLeftNeighbour() != -1) //sendrecv with left neighbour
		{
			commISend(domainStart, domainStart,
						loadDist->getLeftNeighbour(), request[commNum++], subsetToExchange);
			commIRecv(domainStart-1, domainStart-1,
						loadDist->getLeftNeighbour(), request[commNum++], subsetToExchange);
		}
		if (loadDist->getRightNeighbour() != -1) //and with right n.
		{
			commISend(domainEnd, domainEnd,
						loadDist->getRightNeighbour(), request[commNum++], subsetToExchange);
			commIRecv(domainEnd+1, domainEnd+1,
						loadDist->getRightNeighbour(), request[commNum++], subsetToExchange);
		}
		MPI_Waitall(commNum, request, status);
	}

	void commGatherTo(CTable<T> &dest)
	//gathers the table *this (which is distributed according to this->loadDist)
	//into table dest on root proc
	{
		if (loadDist->getRank() == 0)
		{
			//receive other non-empty subdomains
			MPI_Request *request = new MPI_Request[loadDist->getNumProcs()];
			MPI_Status *status = new MPI_Status[loadDist->getNumProcs()];
			int commNum = 0;
			for (int proc = 1; proc < loadDist->getNumProcs(); proc++)
				if (loadDist->getDomainSize(proc) > 0)
					dest.commIRecv(loadDist->getLeftNeighbour(proc) == -1
									? loadDist->getDomainStart(proc)-1
									: loadDist->getDomainStart(proc), 
								loadDist->getRightNeighbour(proc) == -1
									? loadDist->getDomainEnd(proc)+1 
									: loadDist->getDomainEnd(proc),
									proc, request[commNum++]);

			//copy root's subdomain of *this to dest
			if (loadDist->getDomainSize() > 0)
			{
				for (int z = 1; z <= (loadDist->getRightNeighbour() == -1
									? loadDist->getDomainEnd()+1 
									: loadDist->getDomainEnd())
								; z++)
					for (int y = 1; y <= sizeY; y++)
						for (int x = 1; x <= sizeX; x++)
							dest(x, y, z) = access(x, y, z);
			}

			//wait for end of comm
			MPI_Waitall(commNum, request, status);
			delete request;
			delete status;
		} else if (loadDist->getDomainSize() > 0) //if subdomain non-empty, send it to root
				commSend(loadDist->getLeftNeighbour() == -1
							? domainStart-1 
							: domainStart, 
						loadDist->getRightNeighbour() == -1
							? domainEnd+1 
							: domainEnd, 0 );
	}

	void commScatterFrom(CTable<T> &src)
	//scatters table src, which resides on root proc,
	//into distributed table *this (which is distributed according to this->loadDist)
	//also sends border values so that commExchangeBorders is not needed afterwards
	{
		if (loadDist->getRank() == 0)
		{
			//send other non-empty subdomains including borders
			MPI_Request *request = new MPI_Request[loadDist->getNumProcs()];
			MPI_Status *status = new MPI_Status[loadDist->getNumProcs()];
			int commNum = 0;
			for (int proc = 1; proc < loadDist->getNumProcs(); proc++)
				if (loadDist->getDomainSize(proc) > 0)
					src.commISend(loadDist->getDomainStart(proc)-1,
									loadDist->getDomainEnd(proc)+1, 
									proc, request[commNum++]);

			//copy root's subdomain of *this from src
			if (loadDist->getDomainSize() > 0)
			{
				for (int z = 1; z <= domainEnd+1; z++)
					for (int y = 1; y <= sizeY; y++)
						for (int x = 1; x <= sizeX; x++)
							access(x, y, z) = src(x, y, z);
			}

			//wait for end of comm
			MPI_Waitall(commNum, request, status);
			delete request;
			delete status;
		} else if (loadDist->getDomainSize() > 0) //if subdomain non-empty, receive it from root
				commRecv(domainStart-1, 
						domainEnd+1, 0 );
	}

	void printGathered(std::ostream &str);
	//gathers table on process 0 and prints its contents to str
	
	void printDistributed(std::ostream &str);
	//prints table contents of all processes without gathering
	//	using a bunch of MPI_Barriers to ensure correct order

}; //end class CTable

}; //end namespace CommonLib

#endif

