#include "table.h"

namespace CommonLib
{

template<>
void CTable<double>::printGathered(std::ostream &str)
//gathers table on process 0 and prints its contents to str
{
	CTable<double> wholeTable;
	if (loadDist->getRank() == 0)
	{
		CNondistributedLoadDistribution nonDistr(sizeZ, loadDist->getRank(), loadDist->getNumProcs());
		wholeTable.allocate(sizeX, sizeY, sizeZ, &nonDistr);
		commGatherTo(wholeTable);
		char buf[1000];
		sprintf(buf, "CTable<double> %d x %d x %d\n", sizeX, sizeY, sizeZ);
		str << buf;
		for (int z = 1; z <= sizeZ; z++)
		{
			str << "Plane " << z << std::endl;
			for (int y = 1; y <= sizeY; y++)
			{
				for (int x = 1; x <= sizeX; x++)
				{
					sprintf(buf, "%14.7lf ", wholeTable(x, y, z));
					str << buf;
				}
				str << std::endl;
			}
			str << std::endl;
		}
	} else {
		commGatherTo(wholeTable);
	}
}

template <>
void CTable<double>::printDistributed(std::ostream &str)
//prints table contents of all processes without gathering
//	using a bunch of MPI_Barriers to ensure correct order
{
	char buf[1000];
	for (int i = 0; i < loadDist->getNumProcs(); i++)
	{
		MPI_Barrier(MPI_COMM_WORLD);
		if (loadDist->getRank() == i)
		{
			sprintf(buf, "CTable<double> %d x %d x %d; %d-%d on proc %d\n",
					sizeX, sizeY, sizeZ,
					loadDist->getDomainStart()-1, loadDist->getDomainEnd()+1, i);
			str << buf;
			for (int z = loadDist->getDomainStart()-1; z <= loadDist->getDomainEnd()+1; z++)
			{
				str << "Plane " << z << std::endl;
				for (int y = 1; y <= sizeY; y++)
				{
					for (int x = 1; x <= sizeX; x++)
					{
						sprintf(buf, "%14.7lf ", access(x, y, z));
						str << buf;
					}
					str << std::endl;
				}
				str << std::endl;
			}
		}
		str.flush();
	}
}

}; //end namespace CommonLib


