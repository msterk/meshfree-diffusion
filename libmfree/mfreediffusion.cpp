#include "mfreediffusion.h"

#include "../mygsl/gsl_linalg.h"
#include <mpi.h>
#include <iostream>

#include <vector>
#include <algorithm>

namespace mfree
{

void CMFreeDiffusion::importNodes(const CNodeTable &_nodesTable)
//copies nodes from given table and sorts them by index
{
	nodesTable = _nodesTable;
	nodesTable.sortNodes(CNodeTable::sortOrderByIndex);
}

void CMFreeDiffusion::createNodes()
//calls CNodeTable::createNodes
//reads all parameters from options
{
	totalNumNodes = nodesTable.createOrReadNodes(options);
	subdomain = CRect(0,1,0,1);
}

void CMFreeDiffusion::prepareNodeTree()
//builds kd tree from nodes stored in nodesTable
{
	//make separate copies of nodesTable for sorting by x and by y
	//	because the original node ordering in nodesTable must not be changed
	CNodeTable nodesX, nodesY;
	nodesX = nodesTable;
	nodesY = nodesTable;
	nodesX.sortNodes(CNodeTable::sortOrderByX);
	nodesY.sortNodes(CNodeTable::sortOrderByY);
	int height = nodesTree.buildKDTree(nodesTable.getNumElems(), nodesX.getRawPointer(), nodesY.getRawPointer(), 0);
	//printf("tree height = %d (theor. minimum is %lf)\n", height, log(nodesTable.getNumElems()+1.0)/log(2.0)-1);
}

void CMFreeDiffusion::exportNodes(char *filename) const
//exports nodes to a struct named filename in a new Matlab script named s_filename.m
//creates file, writes preamble (clear filename; global filename;) and node matrix
{
	char filename_m[1000], varName[100];
	sprintf(filename_m, "s_%s.m", filename);
	FILE *f = fopen(filename_m, "w");
	const char *domainStruct = "domain";
	fprintf(f, "clear %s;\nglobal %s;\n", domainStruct, domainStruct);
	sprintf(varName, "%s.nodes", domainStruct);
	nodesTable.writeToMFile(f, varName);
	fclose(f);
}

void CMFreeDiffusion::distributeNodes()
//distributes the nodes among processes
//kd tree containing subdomain nodes is also prepared on all processes because building it can be overlapped with communication
//for hierarchical distribution, avgDistMesh should be initialized before and kd tree built
{
	//on root:
	//for each other processor: calculate its subdomain and ISend it its nodes
	//for own subdomain: mark which nodes are outside it and should thus be deleted
	//wait for others to recv data
	//delete marked nodes
	//prepare kd tree from remaining nodes
	//
	//on others: just recv nodes & create kd tree from them

	int numRunningComm;		//total number of communications started by root (which should later be waited for)
	MPI_Request *iRequest;	//requests for all communications started by root

	CNodeTable *subdomainNodes; //nodes, rects & num of nodes for subdomain of each process
	CRect *subdomains;			//subdomainNodes is only used with hierarchical distribution because nodes belonging
	int *subdomainCounts;		//to a single subdomain in 1D distribution can be collected by sorting all nodes by y-coordinate

	int nodesPerProc,	//avg number of nodes per proc for 1D distribution
		firstNodeInSubdomain, lastNodeInSubdomain,	//...of some  process for 1D distribution
		firstNodeToSend, lastNodeToSend;			//nodes in subdomain + border for 1D distribution

	int firstNodeToDelete = -1; //which nodes to delete from root in case of 1D distribution
								//if >1, this marks that nodes from toDelete to end should be deleted from nodesTable

	//tell all procs the total num of nodes
	CommonLib::commBroadcast(totalNumNodes);

	if (myRank == 0)
	{
		//alloc tables & calculate border size (d_max)
		subdomainNodes = new CNodeTable[numProcs];
		subdomains = new CRect[numProcs];
		subdomainCounts = new int[numProcs];
		iRequest = new MPI_Request[numProcs*3];
		numRunningComm = 0;
		double borderSize = (options.dataDistrType == MFreeOptions::DataDistrHierarchical
			? avgDistMesh.getMaxAvgDist()*(options.alphaQuad/2+options.alphaSupport/2)
			: options.dataDistrType == MFreeOptions::DataDistr1D 
				? options.alphaSupport/sqrt(totalNumNodes*0.5) : 1);

		if (options.dataDistrType != MFreeOptions::DataDistrHierarchical)
		{
			nodesTable.sortNodes(CNodeTable::sortOrderByY);
			nodesPerProc = totalNumNodes/numProcs;
		}

		//calc subdomain for each non-root process & send its nodes
		for (int proc = numProcs-1; proc > 0; proc--)
		{
			if (options.dataDistrType == MFreeOptions::DataDistrHierarchical)
			{
				subdomains[proc] = nodesTree.getSubdomainRect(proc, numProcs);
				subdomainCounts[proc] = 
						nodesTree.findNodesInRect(subdomains[proc] + borderSize, subdomainNodes[proc]);
			/*	std::cout << "to p" << proc << ": " << subdomainCounts[proc] << "nodes";
				CNodeTable lj; int bu;
				bu = nodesTree.findNodesInRect(subdomains[proc], lj);
				std::cout << "(" << bu << " real)\n";*/
			} else {
				firstNodeInSubdomain = proc*nodesPerProc;
				lastNodeInSubdomain = (proc == numProcs-1 ? totalNumNodes-1 : (proc+1)*nodesPerProc-1);
				subdomains[proc] = CRect(0,1,nodesTable[firstNodeInSubdomain].y,
											nodesTable[lastNodeInSubdomain].y);
				if (options.dataDistrType == MFreeOptions::DataDistr1D)
				{
					firstNodeToSend = firstNodeInSubdomain;
					while (firstNodeToSend > 0
							&& nodesTable[firstNodeToSend-1].y > subdomains[proc].yMin-borderSize)
						firstNodeToSend--;
					lastNodeToSend = lastNodeInSubdomain;
					while (lastNodeToSend < totalNumNodes-1
							&& nodesTable[lastNodeToSend+1].y < subdomains[proc].yMax+borderSize)
						lastNodeToSend++;
				} else {
					firstNodeToSend = 0;
					lastNodeToSend = totalNumNodes-1;
				}
				subdomainCounts[proc] = lastNodeToSend-firstNodeToSend+1;
			/*	std::cout << "to p" << proc << ": " << subdomainCounts[proc] << "nodes";
				std::cout << "(" << lastNodeInSubdomain-firstNodeInSubdomain+1 << " real)\n";*/
			}

			if (options.verboseLevel >= 5)
				//print subdomain of each process
				std::cout << "to p" << proc << ": subdomain " << subdomains[proc] << std::endl;

			//first send subdomain rectangle & num of nodes that will be sent
			MPI_Isend(&subdomains[proc], sizeof(subdomains[proc]), MPI_BYTE, proc,
						42, MPI_COMM_WORLD, &iRequest[numRunningComm++]);
			MPI_Isend(&subdomainCounts[proc], 1, MPI_INT, proc,
						43, MPI_COMM_WORLD, &iRequest[numRunningComm++]);

			//send nodes themselves
			if (options.dataDistrType == MFreeOptions::DataDistrHierarchical)
				MPI_Isend(subdomainNodes[proc].getRawPointer(), sizeof(CNode)*subdomainCounts[proc], MPI_BYTE, proc,
						44, MPI_COMM_WORLD, &iRequest[numRunningComm++]);
			else
				MPI_Isend(nodesTable.getRawPointer()+firstNodeToSend, sizeof(CNode)*subdomainCounts[proc], MPI_BYTE, proc,
						44, MPI_COMM_WORLD, &iRequest[numRunningComm++]);
		}

		//"send" nodes from root to root
		if (options.dataDistrType == MFreeOptions::DataDistrHierarchical)
		{
			//construct new tree from nodes in subdomain+border, delete other nodes
			subdomain = nodesTree.getSubdomainRect(0, numProcs);
			nodesTable.emptyTable();
			nodesTree.findNodesInRect(subdomain + borderSize, nodesTable);
			nodesTree.emptyTree();
			prepareNodeTree();
		} else {
			//mark which nodes are to be deleted
			firstNodeInSubdomain = 0;
			lastNodeInSubdomain = nodesPerProc-1;
			subdomain = CRect(0,1,nodesTable[firstNodeInSubdomain].y,
									nodesTable[lastNodeInSubdomain].y);				
			if (options.dataDistrType == MFreeOptions::DataDistr1D)
			{
				lastNodeToSend = lastNodeInSubdomain;
				while (lastNodeToSend < totalNumNodes-1
						&& nodesTable[lastNodeToSend+1].y < subdomain.yMax+borderSize)
					lastNodeToSend++;
				if (lastNodeToSend < totalNumNodes-1)
					firstNodeToDelete = lastNodeToSend+1;
			}
			if (firstNodeToDelete == -1)
				prepareNodeTree();
		}
	} else {
		//on non-root procs: receive subdomain and nodes, prepare tree
		MPI_Status status;
		MPI_Recv(&subdomain, sizeof(subdomain), MPI_BYTE, 0, 42, MPI_COMM_WORLD, &status);
		nodesTable.recvNodes(0);
		prepareNodeTree();
	}

//	CommonLib::commDbgOutput(subdomain, "subdomain");

	if (myRank == 0)
	{
		//on root: wait for all communication to finish, delete temporary tables...
		MPI_Status *iStatus;
		iStatus = new MPI_Status[numProcs*3];
		MPI_Waitall(numRunningComm, iRequest, iStatus);
		delete []iStatus;
		delete []iRequest;
		delete []subdomainNodes;
		delete []subdomains;
		delete []subdomainCounts;
		if (firstNodeToDelete != -1)
		{
			//delete marked nodes
			nodesTable.deleteElems(firstNodeToDelete, nodesTable.getNumElems()-1);
			//only in this case must new tree be prepared at the end;
			//in all other cases we prepare it as soon as possible so construction overlaps with communication
			prepareNodeTree();
		}
	}
}

bool CMFreeDiffusion::testTree() const
//tests tree for correctness
//tests finding in various rects & circles
//tests division to subdomains
{
	nodesTree.testTreeCorrectness();
	for (int i = 1; i < 10; i++)
	{
		//test findNodesInRect() with random rectangles
		CRect rect;
		rect.xMin = CommonLib::randFromRange(0.0,1.0);
		rect.xMax = CommonLib::randFromRange(0.0,1.0);
		if (rect.xMin > rect.xMax)
			CommonLib::swap(rect.xMin, rect.xMax);
		rect.yMin = CommonLib::randFromRange(0.0,1.0);
		rect.yMax = CommonLib::randFromRange(0.0,1.0);
		if (rect.yMin > rect.yMax)
			CommonLib::swap(rect.yMin, rect.yMax);
		CNodeTable nodesAllT, nodesAllF;
		nodesTree.findNodesInRect(rect, nodesAllT);
		for (int j = 0; j < nodesTable.getNumElems(); j++)
			if (nodesTable[j].isInRect(rect))
				nodesAllF.appendElem(nodesTable[j]);
		nodesAllT.sortNodes(mfree::CNodeTable::sortOrderByIndex);
		nodesAllF.sortNodes(mfree::CNodeTable::sortOrderByIndex);
		if (nodesAllT.getNumElems() != nodesAllF.getNumElems())
		{
			printf("Different num of nodes!\n");
			return false;
		}
		for (int i = 0; i < nodesAllT.getNumElems(); i++)
			if (nodesAllT[i] != nodesAllF[i])
			{
				printf("Different nodes!\n");
				return false;
			}
	}

	for (int i = 1; i < 10; i++)
	{
		//test findNodesInCircle() with random circles
		CPoint point;
		double diameter;
		point.x = CommonLib::randFromRange(0.0,1.0);
		point.y = CommonLib::randFromRange(0.0,1.0);
		diameter = CommonLib::randFromRange(0.0,1.42);
		CNodeTable nodesAllT, nodesAllF;
		nodesTree.findNodesInCircle(point, diameter, nodesAllT);
		for (int j = 0; j < nodesTable.getNumElems(); j++)
			if ((nodesTable[j]-point).length() <= diameter/2)
				nodesAllF.appendElem(nodesTable[j]);
		nodesAllT.sortNodes(mfree::CNodeTable::sortOrderByIndex);
		nodesAllF.sortNodes(mfree::CNodeTable::sortOrderByIndex);
		if (nodesAllT.getNumElems() != nodesAllF.getNumElems())
		{
			printf("Different num of nodes!\n");
			return false;
		}
		for (int i = 0; i < nodesAllT.getNumElems(); i++)
			if (nodesAllT[i].index != nodesAllF[i].index)
			{
				printf("Different nodes!\n");
				return false;
			}
	}

	if (numProcs == 1 && nodesTable.getNumElems() > 64)
	{
		//test dividing to subdomains
		const int testNumProcs = 16;
		CRect subdomains[testNumProcs];

		CNodeTable allNodes;
		for (int proc = 0; proc < testNumProcs; proc++)
		{
			CNodeTable lj;
			subdomains[proc] = nodesTree.getSubdomainRect(proc, testNumProcs);
			nodesTree.findNodesInRect(subdomains[proc], lj);
			if (options.verboseLevel >= 4)
				std::cout << "proc " << proc << " has domain " << subdomains[proc]
							<< "\n\twith " << lj.getNumElems() << "nodes: ";
			nodesTree.findNodesInRect(subdomains[proc], allNodes);
			if (options.verboseLevel >= 6)
			{
				lj.sortNodes(CNodeTable::sortOrderByIndex);
				for (int i = 0; i < lj.getNumElems(); i++)
					std::cout << lj[i].index << " ";
				std::cout << std::endl;
			}
		}
		//assert(allNodes.getNumElems() == nodesTable.getNumElems());
		allNodes.sortNodes(CNodeTable::sortOrderByIndex);
		for (int i = 1; i < allNodes.getNumElems(); i++)
			if (allNodes[i].index != allNodes[i-1].index+1)
				std::cout << allNodes[i].index << " incorrectly scattered ";

		assert(allNodes.getNumElems() == nodesTable.getNumElems());
		for (int i = 0; i < allNodes.getNumElems(); i++)
			assert(allNodes[i].index == i);
	}

	return true;
}

class CSparseVectorInputTable: public C1DTable<sparse::CSparseVectorInput, 100>
//table of CSparseVectorInputs
{
public:
	void append(int col, double value)
	{
		appendElem(sparse::CSparseVectorInput(col, value));
	}
}; //end class CSparseVectorInputTable

//test function and its derivatives
double Wquad(CPoint node, CPoint quadPoint, double quadDomainSize)
{
	CPoint dist = (quadPoint - node)/(quadDomainSize/2);
	if (fabs(dist.x) < 1 && fabs(dist.y) < 1)
		return (1-dist.x*dist.x) * (1-dist.y*dist.y);
	else
		return 0;
}
double Wquad_x(CPoint node, CPoint quadPoint, double quadDomainSize)
{
	CPoint dist = (quadPoint - node)/(quadDomainSize/2);
	if (fabs(dist.x) < 1 && fabs(dist.y) < 1)
		return -2*dist.x/(quadDomainSize/2) * (1-dist.y*dist.y);
	else
		return 0;
}
double Wquad_y(CPoint node, CPoint quadPoint, double quadDomainSize)
{
	CPoint dist = (quadPoint - node)/(quadDomainSize/2);
	if (fabs(dist.x) < 1 && fabs(dist.y) < 1)
		return (1 - dist.x*dist.x) * (-2)*dist.y/(quadDomainSize/2);
	else
		return 0;
}

int findIntersections(const CRect &rect, const CPoint &c, double r, double *angles)
//finds intersections of rectangle and circle
//writes angles of intersections (as seen from center) and returns no of intersections
//(there can be up to 8)
//*angles should thus have space for at least 8 numbers
{
	double intX, intY;
	int numIntersections = 0;
	if (CommonLib::in(rect.yMin, c.y-r, c.y+r)) //lower edge
	{
		//solve second degree equation to obtain up to two intersections
		intX = c.x - sqrt(r*r - (rect.yMin-c.y)*(rect.yMin-c.y));
		if (CommonLib::in(intX, rect.xMin, rect.xMax))
			angles[numIntersections++] = atan2(rect.yMin-c.y, intX-c.x);
		intX = c.x + sqrt(r*r - (rect.yMin-c.y)*(rect.yMin-c.y));
		if (CommonLib::in(intX, rect.xMin, rect.xMax))
			angles[numIntersections++] = atan2(rect.yMin-c.y, intX-c.x);
	}
	if (CommonLib::in(rect.yMax, c.y-r, c.y+r)) //upper edge
	{
		intX = c.x - sqrt(r*r - (rect.yMax-c.y)*(rect.yMax-c.y));
		if (CommonLib::in(intX, rect.xMin, rect.xMax))
			angles[numIntersections++] = atan2(rect.yMax-c.y, intX-c.x);
		intX = c.x + sqrt(r*r - (rect.yMax-c.y)*(rect.yMax-c.y));
		if (CommonLib::in(intX, rect.xMin, rect.xMax))
			angles[numIntersections++] = atan2(rect.yMax-c.y, intX-c.x);
	}
	if (CommonLib::in(rect.xMin, c.x-r, c.x+r)) //left edge
	{
		intY = c.y - sqrt(r*r - (rect.xMin-c.x)*(rect.xMin-c.x));
		if (CommonLib::in(intY, rect.yMin, rect.yMax))
			angles[numIntersections++] = atan2(intY-c.y, rect.xMin-c.x);
		intY = c.y + sqrt(r*r - (rect.xMin-c.x)*(rect.xMin-c.x));
		if (CommonLib::in(intY, rect.yMin, rect.yMax))
			angles[numIntersections++] = atan2(intY-c.y, rect.xMin-c.x);
	}
	if (CommonLib::in(rect.xMax, c.x-r, c.x+r)) //right edge
	{
		intY = c.y - sqrt(r*r - (rect.xMax-c.x)*(rect.xMax-c.x));
		if (CommonLib::in(intY, rect.yMin, rect.yMax))
			angles[numIntersections++] = atan2(intY-c.y, rect.xMax-c.x);
		intY = c.y + sqrt(r*r - (rect.xMax-c.x)*(rect.xMax-c.x));
		if (CommonLib::in(intY, rect.yMin, rect.yMax))
			angles[numIntersections++] = atan2(intY-c.y, rect.xMax-c.x);
	}
	return numIntersections;
}

#define MAX_SUP_DOMAIN_COUNT 100 //max number of nodes in support domain of any point

class CQuadPointsWeights
//stores abscissae (points) and weights for 1D Gaussian integration over interval (0,1)
//the required degree (1,3,...,13 for single-interval integration and 2*1, 2*3, 2*5 for double-interval integration)
//is passed to constructor, which initializes the object with the right constants
//points and weights can then be directly read from *quadPoints and *quadWeights
{
protected:
	static double quadPoints1 [], quadWeights1 [];
	static double quadPoints2 [], quadWeights2 [];
	static double quadPoints3 [], quadWeights3 [];
	static double quadPoints6 [], quadWeights6 [];
	static double quadPoints5 [], quadWeights5 [];
	static double quadPoints10[], quadWeights10[];
	static double quadPoints7 [], quadWeights7 [];
	static double quadPoints9 [], quadWeights9 [];
	static double quadPoints11[], quadWeights11[];
	static double quadPoints13[], quadWeights13[];

public:
	static const int MaxNumQuadPoints = 10;
	double *quadPoints, *quadWeights;	//points & weights to use; access directly
	int numQuadPoints;					//num of quad points used

	CQuadPointsWeights(int gaussDegree)
	//copies the right constants to quadPoints & quadWeights
	{
		numQuadPoints = (gaussDegree % 2 == 1) ? (gaussDegree+1)/2 : gaussDegree/2+1;
		switch (gaussDegree)
		{
			case 1:
				quadPoints = quadPoints1;
				quadWeights = quadWeights1;
				break;
			case 3:
				quadPoints = quadPoints3;
				quadWeights = quadWeights3;
				break;
			case 5:
				quadPoints = quadPoints5;
				quadWeights = quadWeights5;
				break;
			case 7:
				quadPoints = quadPoints7;
				quadWeights = quadWeights7;
				break;
			case 9:
				quadPoints = quadPoints9;
				quadWeights = quadWeights9;
				break;
			case 11:
				quadPoints = quadPoints11;
				quadWeights = quadWeights11;
				break;
			case 13:
				quadPoints = quadPoints13;
				quadWeights = quadWeights13;
				break;
			case 2:
				quadPoints = quadPoints2;
				quadWeights = quadWeights2;
				break;
			case 6:
				quadPoints = quadPoints6;
				quadWeights = quadWeights6;
				break;
			case 10:
				quadPoints = quadPoints10;
				quadWeights = quadWeights10;
				break;
			default:
				throw new MFreeException(__LINE__, "QuadPointsWeights::fillData: invalid gauss degree");
				break;
		}
	}
};

//constants for CQuadPointsWeights; error-free (that is, checked in 3 different sources)
double CQuadPointsWeights::quadPoints1[]  = {0.5};
double CQuadPointsWeights::quadWeights1[] = {1};

double CQuadPointsWeights::quadPoints2[]  = {0.25, 0.75};
double CQuadPointsWeights::quadWeights2[] = {0.5, 0.5};

double CQuadPointsWeights::quadPoints3[]  = {0.21132486540519, 0.78867513459481};
double CQuadPointsWeights::quadWeights3[] = {0.5, 0.5};

double CQuadPointsWeights::quadPoints6[]  = {0.10566243270260, 0.39433756729740, 0.60566243270260, 0.89433756729740};
double CQuadPointsWeights::quadWeights6[] = {0.25, 0.25, 0.25, 0.25};

double CQuadPointsWeights::quadPoints5[]  = {0.11270166537926, 0.50000000000000, 0.88729833462074};
double CQuadPointsWeights::quadWeights5[] = {0.27777777777778, 0.44444444444444, 0.27777777777778};

double CQuadPointsWeights::quadPoints10[]  = {0.05635083268963, 0.25000000000000, 0.44364916731037, 0.55635083268963, 0.75000000000000, 0.94364916731037};
double CQuadPointsWeights::quadWeights10[] = {0.13888888888889, 0.22222222222222, 0.13888888888889, 0.13888888888889, 0.22222222222222, 0.13888888888889};

double CQuadPointsWeights::quadPoints7[]  = {0.06943184420297, 0.33000947820757, 0.66999052179243, 0.93056815579703};
double CQuadPointsWeights::quadWeights7[] = {0.17392742256873, 0.32607257743127, 0.32607257743127, 0.17392742256873};

double CQuadPointsWeights::quadPoints9[]  = {0.04691007703067, 0.23076534494716, 0.50000000000000, 0.76923465505284, 0.95308992296933};
double CQuadPointsWeights::quadWeights9[] = {0.11846344252809, 0.23931433524968, 0.28444444444444, 0.23931433524968, 0.11846344252809};

double CQuadPointsWeights::quadPoints11[]  = {0.03376524289842, 0.16939530676687, 0.38069040695840, 0.61930959304160, 0.83060469323313, 0.96623475710158};
double CQuadPointsWeights::quadWeights11[] = {0.08566224618958, 0.18038078652407, 0.23395696728635, 0.23395696728635, 0.18038078652407, 0.08566224618958};

double CQuadPointsWeights::quadPoints13[]  = {0.02544604382862, 0.12923440720030, 0.29707742431130, 0.50000000000000, 0.70292257568870, 0.87076559279970, 0.97455395617138};
double CQuadPointsWeights::quadWeights13[] = {0.06474248307573, 0.13985269574744, 0.19091502525638, 0.20897959184091, 0.19091502525638, 0.13985269574744, 0.06474248307573};


void coutTable(CNodeTable& nodes, bool parentheses = true, bool commas = true) {
	if (parentheses)
		std::cout << "[";
	
	std::vector<int> nodeIndices(nodes.getNumElems());
	for (size_t i = 0; i < nodes.getNumElems(); ++i) 
		nodeIndices[i] = (nodes[i].index);
	std::sort(nodeIndices.begin(), nodeIndices.end());
	
	for (size_t i = 0; i < nodeIndices.size(); ++i) {
		std::cout << (i > 0 ? (commas ? ", " : " ") : "") << nodeIndices[i];
	}
	if (parentheses)
		std::cout << "]";
}


void coutSupDomain(const char* supDomainType, int node, int qIndex, CNodeTable& supDomain, double domainSize) {
	//std::cout << supDomainType << " support domain of size " << supDomain.getNumElems() << ", radius=" << domainSize << "[";
	std::cout << supDomainType << " node " << node << " quad point " << qIndex << ": ";
	coutTable(supDomain);
	std::cout << "\n";
}

void coutSupDomain_m(const char* supDomainType, int node, int qIndex, CNodeTable& supDomain, double domainSize) {
	//std::cout << supDomainType << ", ";
	std::cout << node << ", " << qIndex << ", " << domainSize << ", ";
	coutTable(supDomain, false, true);
	std::cout << "\n";
}


void CMFreeDiffusion::constructSystem()
//main part of program: constructs system matrices (A, B) & RHS vector (fstar)
//
//type of boundary conditions (Dirichlet or Neumann) are hard-coded in this function
//	boundary conditions values (whether Dirichlet or Neumann) are taken from node.u
//
//after sequential execution, A, B and fstar will containt the whole system
//after parallel construction of system, each process will have its part of matrices A and B
//during parallel construction of system, each process will allocate space for full
//	vector fstar, but will only put values into components of fstar corresponding to
//	its own nodes; components of fstar corresponding to other nodes will be NaN
{
	int numNodes = nodesTable.getNumElems();

	//allocate space for full vector fstar
	fstar = new double[totalNumNodes];
	for (int i = 0; i < totalNumNodes; i++)
		fstar[i] = CommonLib::NaN;

	CQuadPointsWeights gaussQuad(options.quadDegree);	//1D Gaussian weights & abscissae
	CSparseVectorInputTable AsparseInput, BsparseInput;	//contributions into the matrix line corresponding to current node

	CNodeTable quadSupDomain,	//support nodes of whole quadrature domain
				supDomain;		//support nodes of a single quadrature point
	CRect quadDomain;			//quadrature domain
	double supDomainDiameter, quadDomainSize; //diameters of (circular) support and (square) quadrature domain

	double phi[MAX_SUP_DOMAIN_COUNT], phi_x[MAX_SUP_DOMAIN_COUNT], phi_y[MAX_SUP_DOMAIN_COUNT];
				//values of base functions and their first derivatives

	int matrixARow = 0, matrixBRow = 0; //number of rows in each matrix constucted so far

	A.initialize(numNodes, totalNumNodes, 40); //total matrix size is totalNumNodes x totalNumNodes;
	B.initialize(numNodes, totalNumNodes, 40); //of those, numNodes lines will be allocated locally

	// MD: debugging
	supDomains.resize(numNodes);

	for (int node = 0; node < numNodes; node++)
	{
		// MD: debugging, create 9 quadrature points for each node
		supDomains[node].resize(gaussQuad.numQuadPoints);
		
		//for each node that is stored in this process:

		if (nodesTable[node].distToClosest(subdomain) > 0)
			//skip over nodes outside subdomain
			//(which are only in the table because they might be in sup. domain of some quad point,
			//	i.e. their distance from subdomain is less than d_MAX)
			continue;

		if (myRank == 0 && options.verboseLevel >= 5)
			std::cout << "Proc " << myRank << ", node = " << node << " index " << nodesTable[node].index << std::endl;


		bool internalNodeEquation = true,	//true iff this node's equation is (3.22)
				dirichletEquation;			//true iff !internalNodeEquation AND node is in Dirichlet boundary

		//select the version to be used (1-3)
		#define BOUNDARY_VERSION 1

		#if BOUNDARY_VERSION == 1 //dissertation version: Dirichlet conditions on all boundaries
			if (nodesTable[node].boundary)
			{
				internalNodeEquation = false;
				dirichletEquation = true;
			}
		#elif BOUNDARY_VERSION == 2 //Neumman version 1: Neumann S & W boundaries, implemented with collocation in boundary nodes
			if (nodesTable[node].boundary)
			{
				internalNodeEquation = false;
				dirichletEquation = (nodesTable[node].x == 1 || nodesTable[node].y == 1);
			}
		#elif BOUNDARY_VERSION == 3 //Neumann version 2: Neumann S & W boundaries, implemented with zero heat flow
			if (nodesTable[node].boundary && (nodesTable[node].x == 1 || nodesTable[node].y == 1))
			{
				internalNodeEquation = false;
				dirichletEquation = true;
			}
		#else
			#error Invalid BOUNDARY_VERSION
		#endif

		if (!internalNodeEquation)
		{
			//collocation in boundary nodes
			fstar[nodesTable[node].index] = nodesTable[node].u; //value of BC is taken from inital value of node.u
			supDomainDiameter = getSupDomain(nodesTable[node], supDomain);
			//coutSupDomain("border", node, -1, supDomain, supDomainDiameter);
			assert(supDomain.getNumElems() <= MAX_SUP_DOMAIN_COUNT);

			//depending on type of BC, put the right value as contribution to matrix A
			//matrix B is empty if !internalNodeEquation

			if (dirichletEquation)
			{
				//equation says that value of trial function should be equal to something (eq. 3.18)
				calculatePhi(nodesTable[node], supDomain, supDomainDiameter, false, phi, NULL, NULL);
				for (int supNode = 0; supNode < supDomain.getNumElems(); supNode++)
					AsparseInput.append(supDomain[supNode].index, phi[supNode]);
			} else { 
				//equation says that normal derivative of trial function should be equal to something
				calculatePhi(nodesTable[node], supDomain, supDomainDiameter, true, phi, phi_x, phi_y);
				if (nodesTable[node].x == 0 && nodesTable[node].y > 0)
				{
					//left boundary; phi_x = sth
					for (int supNode = 0; supNode < supDomain.getNumElems(); supNode++)
						AsparseInput.append(supDomain[supNode].index, phi_x[supNode]);
				} else if (nodesTable[node].x > 0 && nodesTable[node].y == 0)
				{
					//lower boundery; phi_y = sth
					for (int supNode = 0; supNode < supDomain.getNumElems(); supNode++)
						AsparseInput.append(supDomain[supNode].index, phi_y[supNode]);
				} else {
					//Neumann corner: phi_x + phi_y = sth
					//	(this is only a guess; something else might work better)
					for (int supNode = 0; supNode < supDomain.getNumElems(); supNode++)
						AsparseInput.append(supDomain[supNode].index, phi_x[supNode]+phi_y[supNode]);
				}
			}

		} else {
			//internal nodes: equation (3.22)
			fstar[nodesTable[node].index] = 0; //no heat sources
			
			// support domain is found but is reduced to its radius in the following code
			//calculate quadrature domain, no part of which can lie outside the global domain
			if (options.constSupNodeCount()) {
				CNodeTable supNodes;
				quadDomainSize = options.betaQuad * getSupDomain(nodesTable[node], supNodes);
				//coutSupDomain("constant size internal", node, -1, supNodes, quadDomainSize);
			} else {
				quadDomainSize = options.alphaQuad * avgDistMesh.interpolateAvgDist(nodesTable[node]);
				//coutSupDomain("interpolated internal", node, -1, nodesTable, quadDomainSize);
			}
			quadDomain.xMin = CommonLib::mmax(0.0, nodesTable[node].x - quadDomainSize/2);
			quadDomain.xMax = CommonLib::mmin(1.0, nodesTable[node].x + quadDomainSize/2);
			quadDomain.yMin = CommonLib::mmax(0.0, nodesTable[node].y - quadDomainSize/2);
			quadDomain.yMax = CommonLib::mmin(1.0, nodesTable[node].y + quadDomainSize/2);

			if (options.preFindQuadSupport)
			{
				//find support nodes of all quadrature domain
				quadSupDomain.emptyTable();
				std::cerr << "pre-find support\n";
				nodesTree.findNodesInCircle(nodesTable[node],
						(options.alphaSupport + sqrt(2.0)*options.alphaQuad)*avgDistMesh.getMaxAvgDist(quadDomain),
						quadSupDomain);
			}

			int quadraturePointIndex1 = 0;
			int quadraturePointIndex2 = 0;
			//integrate K and C within rectangle (but not in hole) by looping over quadrature points
			//put corresponding values as contributions to A and B 
			for (int ix = 0; ix < gaussQuad.numQuadPoints; ix++) 
			{
				CPoint quadPoint;
				quadPoint.x = quadDomain.xMin + gaussQuad.quadPoints[ix]*quadDomain.getWidth();

				//divide vertical line through quadPoint.x
				//(that is, line from (quadPoint.x, quadDomain.yMin) to (quadPoint.x, quadDomain.yMax))
				//into two lines - below and above hole
				//lines will go from yMin1 to yMax1 and from yMin2 to yMax2
				//NaN means a line does not exist
				double yMin1 = CommonLib::NaN, yMax1 = CommonLib::NaN,
						yMin2 = CommonLib::NaN, yMax2 = CommonLib::NaN;
				if (options.rHole > 0 && CommonLib::in(quadPoint.x, options.xyHole.x-options.rHole, options.xyHole.x+options.rHole))
				{
					//find intersections
					double yInters1 = options.xyHole.y - sqrt(options.rHole*options.rHole
															- CommonLib::sqr(quadPoint.x-options.xyHole.x));
					double yInters2 = options.xyHole.y + sqrt(options.rHole*options.rHole
															- CommonLib::sqr(quadPoint.x-options.xyHole.x));
					if (yInters1 > quadDomain.yMin && yInters2 < quadDomain.yMax)
					{
						//quadDomain includes both intersecions - define two lines
						yMin1 = quadDomain.yMin;
						yMax1 = yInters1;
						yMin2 = yInters2;
						yMax2 = quadDomain.yMax;
					} else if (yInters1 <= quadDomain.yMin && yInters2 > quadDomain.yMin && yInters2 < quadDomain.yMax)
					{
						//quadDomain includes only inters2 - define line from inters2 to yMax
						yMin1 = yInters2;
						yMax1 = quadDomain.yMax;
					} else if (yInters1 > quadDomain.yMin && yInters1 < quadDomain.yMax && yInters2 >= quadDomain.yMax)
					{
						//quadDomain includes only inters1 - define line from yMin to inters1
						yMin1 = quadDomain.yMin;
						yMax1 = yInters1;
					} else {
						//quadDomain includes no intersection - leave original line
						//	(this is only true if quadDomain is smaller than hole; otherwise there are
						//	additional possibilities that whole original line is inside hole, but
						//	holes smaller than any quad domain are not supported)
						yMin1 = quadDomain.yMin;
						yMax1 = quadDomain.yMax;
					}
				} else {
					//quadPoint.x is left or right of hole, or hole does not exist - leave original line
					yMin1 = quadDomain.yMin;
					yMax1 = quadDomain.yMax;
				}

				//evaluate line 1
				if (!CommonLib::isNaN(yMin1))
					for (int iy = 0; iy < gaussQuad.numQuadPoints; iy++) 
					{
						quadraturePointIndex1++;
						quadPoint.y = yMin1 + gaussQuad.quadPoints[iy]*(yMax1-yMin1);
						//now quadPoint is fully defined

						//find support domain diameter & support nodes of quadPoint
						if (options.preFindQuadSupport)
							supDomainDiameter = getSupDomain(quadPoint, supDomain, quadSupDomain);
						else
							supDomainDiameter = getSupDomain(quadPoint, supDomain);
						//coutSupDomain_m("line1", node, quadraturePointIndex1, supDomain, supDomainDiameter);
						supDomains[node][iy] = supDomain;
						
						assert(supDomain.getNumElems() <= MAX_SUP_DOMAIN_COUNT);

						//calculate values of base f. and their derivatives
						calculatePhi(quadPoint, supDomain, supDomainDiameter, true, phi, phi_x, phi_y);
						double gaussWeight = gaussQuad.quadWeights[ix]*gaussQuad.quadWeights[iy]
										* quadDomain.getWidth()*(yMax1-yMin1);

						for (int supNode = 0; supNode < supDomain.getNumElems(); supNode++)
						{
							//for each support node of quadPoint:
							//contribution of this supNode for this quadPoint is:
							//gaussWeight*(2*C +/- c*dt*K) for matrices A and B, respectively
							//C and K are from eq. (3.22) but without integral over Gamma_Qi in K

							double K, C;
							C = Wquad(nodesTable[node], quadPoint, quadDomainSize) * phi[supNode];
							K = Wquad_x(nodesTable[node], quadPoint, quadDomainSize) * phi_x[supNode]
								+ Wquad_y(nodesTable[node], quadPoint, quadDomainSize) * phi_y[supNode];

							AsparseInput.append(supDomain[supNode].index,
													gaussWeight * (2*C + options.c*options.timeStep*K));
							BsparseInput.append(supDomain[supNode].index,
													gaussWeight * (2*C - options.c*options.timeStep*K));
						}
					}

				//evaluate line 2 in the same way
				if (!CommonLib::isNaN(yMin2))
					for (int iy = 0; iy < gaussQuad.numQuadPoints; iy++) 
					{
						quadraturePointIndex2++;
						quadPoint.y = yMin2 + gaussQuad.quadPoints[iy]*(yMax2-yMin2);
						if (options.preFindQuadSupport)
							supDomainDiameter = getSupDomain(quadPoint, supDomain, quadSupDomain);
						else
							supDomainDiameter = getSupDomain(quadPoint, supDomain);
							
						//coutSupDomain_m("line2", node, quadraturePointIndex2, supDomain, supDomainDiameter);

						assert(supDomain.getNumElems() <= MAX_SUP_DOMAIN_COUNT);
						calculatePhi(quadPoint, supDomain, supDomainDiameter, true, phi, phi_x, phi_y);
						double gaussWeight = gaussQuad.quadWeights[ix]*gaussQuad.quadWeights[iy]
										* quadDomain.getWidth()*(yMax2-yMin2);

						for (int supNode = 0; supNode < supDomain.getNumElems(); supNode++)
						{
							double K, C;
							C = Wquad(nodesTable[node], quadPoint, quadDomainSize) * phi[supNode];
							K = Wquad_x(nodesTable[node], quadPoint, quadDomainSize) * phi_x[supNode]
								+ Wquad_y(nodesTable[node], quadPoint, quadDomainSize) * phi_y[supNode];

							AsparseInput.append(supDomain[supNode].index,
													gaussWeight * (2*C + options.c*options.timeStep*K));
							BsparseInput.append(supDomain[supNode].index,
													gaussWeight * (2*C - options.c*options.timeStep*K));
						}
					}
			}

			//integrate K on intersecion of quad domain with problem domain boundary
			//	(integral over Gamma_Qi in eq. (3.22)

			//first, build a list of all boundary quad points
			//	together with their weights and outer normal vectors (n_x,n_y)
			CPoint boundQuadPoints[5*gaussQuad.MaxNumQuadPoints];
			double boundGaussWeights[5*gaussQuad.MaxNumQuadPoints],
					n_x[5*gaussQuad.MaxNumQuadPoints], n_y[5*gaussQuad.MaxNumQuadPoints];
			int numBoundaryQuadPoints = 0;

            if (nodesTable[node].x - quadDomainSize/2 < 0)
				//left edge of problem domain
				for (int iy = 0; iy < gaussQuad.numQuadPoints; iy++)
				{
					boundQuadPoints[numBoundaryQuadPoints]
						= CPoint(0, quadDomain.yMin + gaussQuad.quadPoints[iy]*quadDomain.getHeight());
					boundGaussWeights[numBoundaryQuadPoints] = gaussQuad.quadWeights[iy]*quadDomain.getHeight();
					n_x[numBoundaryQuadPoints] = -1;
					n_y[numBoundaryQuadPoints] = 0;
					numBoundaryQuadPoints++;
				}

            if (nodesTable[node].x + quadDomainSize/2 > 1)
				//right edge
				for (int iy = 0; iy < gaussQuad.numQuadPoints; iy++)
				{
					boundQuadPoints[numBoundaryQuadPoints]
						= CPoint(1, quadDomain.yMin + gaussQuad.quadPoints[iy]*quadDomain.getHeight());
					boundGaussWeights[numBoundaryQuadPoints] = gaussQuad.quadWeights[iy]*quadDomain.getHeight();
					n_x[numBoundaryQuadPoints] = 1;
					n_y[numBoundaryQuadPoints] = 0;
					numBoundaryQuadPoints++;
				}

            if (nodesTable[node].y - quadDomainSize/2 < 0)
				//bottom edge
				for (int ix = 0; ix < gaussQuad.numQuadPoints; ix++)
				{
					boundQuadPoints[numBoundaryQuadPoints]
						= CPoint(quadDomain.xMin + gaussQuad.quadPoints[ix]*quadDomain.getWidth(), 0);
					boundGaussWeights[numBoundaryQuadPoints] = gaussQuad.quadWeights[ix]*quadDomain.getWidth();
					n_x[numBoundaryQuadPoints] = 0;
					n_y[numBoundaryQuadPoints] = -1;
					numBoundaryQuadPoints++;
				}

            if (nodesTable[node].y + quadDomainSize/2 > 1)
				//top edge
				for (int ix = 0; ix < gaussQuad.numQuadPoints; ix++)
				{
					boundQuadPoints[numBoundaryQuadPoints]
						= CPoint(quadDomain.xMin + gaussQuad.quadPoints[ix]*quadDomain.getWidth(), 1);
					boundGaussWeights[numBoundaryQuadPoints] = gaussQuad.quadWeights[ix]*quadDomain.getWidth();
					n_x[numBoundaryQuadPoints] = 0;
					n_y[numBoundaryQuadPoints] = 1;
					numBoundaryQuadPoints++;
				}

			if (options.xyHole.distToClosest(quadDomain) < options.rHole)
			{
				//hole edge
				
				//find intersections of quad domain boundary and the hole
				//if quad domain is smaller than the hole, there should be 0 or 2 intersections
				double intersAngles[8];
				int numInters = findIntersections(quadDomain, options.xyHole, options.rHole, intersAngles);
				assert(numInters == 0 || numInters == 2);

				if (numInters == 2)
				{
					//arrange intersection angles so that first < second
					if (intersAngles[0] > intersAngles[1])
						CommonLib::swap(intersAngles[0], intersAngles[1]);
					if (intersAngles[1] - intersAngles[0] > CommonLib::pi)
					{
						intersAngles[0] += 2*CommonLib::pi;
						CommonLib::swap(intersAngles[0], intersAngles[1]);
					}
					//add points on hole edge to list of boundary quad points
					double angleRange = intersAngles[1]-intersAngles[0];
					for (int ic = 0; ic < gaussQuad.numQuadPoints; ic++)
					{
						double angle = intersAngles[0] + gaussQuad.quadPoints[ic]*angleRange;
						boundQuadPoints[numBoundaryQuadPoints]
							= options.xyHole + CPoint(cos(angle), sin(angle)) * options.rHole;
						boundGaussWeights[numBoundaryQuadPoints] = gaussQuad.quadWeights[ic]*options.rHole*angleRange;
						n_x[numBoundaryQuadPoints] = -cos(angle);
						n_y[numBoundaryQuadPoints] = -sin(angle);
						numBoundaryQuadPoints++;
					}
				}							
			}

			//iterate over all boundary quad points
			for (int ib = 0; ib < numBoundaryQuadPoints; ib++)
			{
				//find support domain diameter & support nodes of quadPoint
				if (options.preFindQuadSupport)
					supDomainDiameter = getSupDomain(boundQuadPoints[ib], supDomain, quadSupDomain);
				else
					supDomainDiameter = getSupDomain(boundQuadPoints[ib], supDomain);
				assert(supDomain.getNumElems() <= MAX_SUP_DOMAIN_COUNT);

				//calculate values of base f. and their derivatives
				calculatePhi(boundQuadPoints[ib], supDomain, supDomainDiameter, true, phi, phi_x, phi_y);
				
				for (int supNode = 0; supNode < supDomain.getNumElems(); supNode++)
				{
					//for each support node of this boundQuadPoint:
					//contribution of this supNode for this boundQuadPoint is:
					//-/+ gaussWeight*K*c*dt for matrices A and B, respectively
					//K is integral over Gamma_Qi from eq. (3.22)
					double K = - boundGaussWeights[ib] //minus that is before path integral in K_ij
							* Wquad(nodesTable[node], boundQuadPoints[ib], quadDomainSize) *(
								n_x[ib] * phi_x[supNode] + n_y[ib] * phi_y[supNode] ); 
					AsparseInput.append(supDomain[supNode].index,
											options.c*options.timeStep*K);
					BsparseInput.append(supDomain[supNode].index,
											-options.c*options.timeStep*K);
				}
			}
		}
		if (AsparseInput.getNumElems() > 0)
		{
			//if there are any contributions to A (there should always be some), add them to the matrix
			A.addInputsSort(AsparseInput.getRawPointer(), AsparseInput.getNumElems(),
							matrixARow, nodesTable[node].index);

			//prepare the AsparseInput table for construction of next row
			AsparseInput.emptyTable();
			matrixARow++;
		}
		if (BsparseInput.getNumElems() > 0)
		{
			//if there are any contributions to B (there will be none for boundary collocation nodes), add them to the matrix
			B.addInputsSort(BsparseInput.getRawPointer(), BsparseInput.getNumElems(),
							matrixBRow, nodesTable[node].index);

			//prepare the BsparseInput table for construction of next row
			BsparseInput.emptyTable();
			matrixBRow++;
		}
	}

	//remove any empty rows (although there should be none) from matrices; otherwise, gathering of matrices will fail
	A.removeEmptyRows();
	B.removeEmptyRows();
}

void MPI_Op_join( void *invec, void *inoutvec, int *len, MPI_Datatype *datatype)
//for each component of given two vectors, outputs the component that is not NaN
//also checks that at least one is NaN
{
	double *inV = (double*)invec, *inoutV = (double*)inoutvec;
	for (int i = 0; i < *len; i++)
	{
		assert(CommonLib::isNaN(inV[i]) || CommonLib::isNaN(inoutV[i]));
		if (!CommonLib::isNaN(inV[i]))
			inoutV[i] = inV[i];
	}
}

void CMFreeDiffusion::exportSystem(char *filename) const
//gathers A, B, fstar on root and exports them to an existing Matlab script named s_filename.m
//also exports options
{
	MPI_Op opJoin;
	MPI_Op_create(MPI_Op_join, true, &opJoin);

	if (myRank == 0)
	{
		//gather all data into full matrices wholeA and wholeB
		sparse::CSparseMatrix wholeA(totalNumNodes, totalNumNodes, 40);
		wholeA.gatherFrom(A, numProcs);

		sparse::CSparseMatrix wholeB(totalNumNodes, totalNumNodes, 40);
		wholeB.gatherFrom(B, numProcs);
		double *wholefstar = new double[totalNumNodes];

		//gather fstar
		//for each i, fstar(i) should be non-NaN on exactly one processor before reduce
		//after reduce, wholefstar is a vector of non-NaN values
		MPI_Reduce(fstar, wholefstar, totalNumNodes, MPI_DOUBLE, opJoin, 0, MPI_COMM_WORLD);

		//export options, A, B, and f
		char filename_m[1000], varName[100];
		sprintf(filename_m, "s_%s.m", filename);
		FILE *f = fopen(filename_m, "a+");

		options.writeToMFile(f, filename, totalNumNodes);

		sprintf(varName, "%s.A", filename);
		wholeA.writeToMFile(f, varName);

		sprintf(varName, "%s.B", filename);
		wholeB.writeToMFile(f, varName);

		fprintf(f, "%s.fstar=[];\n", filename);
		for (int i = 0; i < totalNumNodes; i++)
			fprintf(f, "%s.fstar(%d)=%.15lf;\n", filename, i+1, wholefstar[i]);
			
		for (int i = 0; i < totalNumNodes; i++)
			fprintf(f, "%s.sup_domain_debug(%d)=%d;\n", filename, i+1, (int)wholefstar[i]);
		fclose(f);

		//free mem
		delete []wholefstar;
	} else {
		//send data to root
		A.sendToRoot(0);
		B.sendToRoot(0);
		MPI_Reduce(fstar, NULL, totalNumNodes, MPI_DOUBLE, opJoin, 0, MPI_COMM_WORLD);
	}

	MPI_Op_free(&opJoin);
}

void CAvgDistMesh::initialize(const CKDTree &nodesTree, const MFreeOptions &options)
//evaluates average nodal distance between nodes given in nodesTree on the regular mesh
//not parallelism-aware: if you want only one process to evaluate distances,
//	you should call initialize() just on that on, then call CAvgDistMesh::broadcast()
{
	//default size: sqrt(N)/5, but at least 5 (if no hole) or at least 10 (if hole)
	allocate( CommonLib::mmax(options.rHole > 0 ? 10 : 5, CommonLib::round(sqrt(nodesTree.getNumNodesInSubtree()*1.0)/5.0)) );

	//size as in Matlab version (for comparison of results): same density as the nodes
	//allocate( CommonLib::round(sqrt( nodesTree.getNumNodesInSubtree() / (1 - CommonLib::pi*options.rHole*options.rHole) )) );

	//calc nI using eq. (4.3)
	int reqSupDomainCount = int(ceil(1 + sqrt(CommonLib::pi)*options.alphaSupport
								+ CommonLib::pi/4*options.alphaSupport*options.alphaSupport));

	//for every point in regular mesh, evaluate & store avgDist
	CPoint p;
	maxAvgDist = -1;
	for (int ix = 0; ix < size; ix++)
		for (int iy = 0; iy < size; iy++)
		{
			p = CPoint(ix*dx, iy*dx);
			double estSupDomSize = dx*options.alphaSupport;
			CNodeTable supDomain;

			do {
				supDomain.emptyTable();
				nodesTree.findNodesInCircle(p, estSupDomSize, supDomain);
				estSupDomSize *= 1.5;
			} while (supDomain.getNumElems() < reqSupDomainCount+1);

			for (int i = 0; i < supDomain.getNumElems(); i++)
				supDomain[i].dist = (supDomain[i] - p).length();
			supDomain.sortNodes(CNodeTable::sortOrderByDist);
			double avgDistNearP = sqrt(CommonLib::pi)
							* (supDomain[reqSupDomainCount-1].dist + supDomain[reqSupDomainCount].dist)/2
							/ (sqrt((double)reqSupDomainCount)-1);
			(*this)(p) = avgDistNearP;
			maxAvgDist = CommonLib::mmax(maxAvgDist, avgDistNearP);
		}
}

void CAvgDistMesh::broadcast(int myRank)
//broadcasts the whole table of evaluated avg distances to all processors
{
	//broadcast root's size to all; allocate table with new size on all but root
	int _size;
	if (myRank == 0)
		_size = size;
	CommonLib::commBroadcast(_size); 
	if (myRank != 0)
		allocate(_size);

	//broadcast table contents line-by-line
	CommonLib::commBroadcast(maxAvgDist);
	for (int ix = 0; ix < size; ix++)
       	CommonLib::commBroadcast(avgDist[ix][0], size);
}

inline double Wmls(double dist, double supDomainDiameter)
//MLS weight of a node with given dist from point of interest, whose supDomainDiam must also be known
{
	double d = dist / (supDomainDiameter/2);
	return d >= 1 ? 0 : 1 - 6*d*d + 8*d*d*d - 3*d*d*d*d;
}

inline double Wmls(CPoint node, CPoint point, double supDomainDiameter)
//MLS weight of a node as a support node of a given point of interest, whose supDomainDiam must also be known
{
	return Wmls((node-point).length(), supDomainDiameter);
}

inline double Wmls_x(CPoint node, CPoint point, double supDomainDiameter)
//derivative of Wmls() on point.x
{
	double d = (node-point).length() / (supDomainDiameter/2);
	return d >= 1 ? 0 : (-12 + 24*d - 12*d*d) * (point.x-node.x)/(supDomainDiameter/2);
}

inline double Wmls_y(CPoint node, CPoint point, double supDomainDiameter)
//derivative of Wmls() on point.x
{
	double d = (node-point).length() / (supDomainDiameter/2);
	return d >= 1 ? 0 : (-12 + 24*d - 12*d*d) * (point.y-node.y)/(supDomainDiameter/2);
}

template <int MLSdegree> inline void pmls(CPoint point, double *p);
//calculates all monomes of given degree; p must have enough space

template<>
inline void pmls<1>(CPoint point, double *p)
//implementation for degree 1
{
	p[0] = 1;
	p[1] = point.x;
	p[2] = point.y;
}
template<>
inline void pmls<2>(CPoint point, double *p)
//implementation for degree 2
{
	p[0] = 1;
	p[1] = point.x;
	p[2] = point.y;
	p[3] = point.x*point.x;
	p[4] = point.x*point.y;
	p[5] = point.y*point.y;
}

template <int MLSdegree> inline void pmls_xy(CPoint point, double *p_x, double *p_y);
//calculates derivatives of all monomes; p_x and p_y must have enough space

template<>
inline void pmls_xy<1>(CPoint point, double *p_x, double *p_y)
//implementation for degree 1
{
	p_x[0] = 0;			p_y[0] = 0;
	p_x[1] = 1;			p_y[1] = 0;
	p_x[2] = 0;			p_y[2] = 1;
}
template<>
inline void pmls_xy<2>(CPoint point, double *p_x, double *p_y)
//implementation for degree 2
{
	p_x[0] = 0;			p_y[0] = 0;
	p_x[1] = 1;			p_y[1] = 0;
	p_x[2] = 0;			p_y[2] = 1;
	p_x[3] = 2*point.x;	p_y[3] = 0;
	p_x[4] = point.y;	p_y[4] = point.x;
	p_x[5] = 0;			p_y[5] = 2*point.y;
}

void printfgslmatrix(FILE *f, gsl_matrix *A)
//for test output
{
	for (int i = 0; i < (int)A->size1; i++)
		for (int j = 0; j < (int)A->size2; j++)
			fprintf(f, "A(%d,%d) = %f\n", i, j, gsl_matrix_get(A, i, j));
	fprintf(f, "\n");
}

void printfgslvector(FILE *f, gsl_vector *x)
//for test output
{
	for (int i = 0; i < (int)x->size; i++)
		fprintf(f, "x(%d) = %f\n", i, gsl_vector_get(x, i));
	fprintf(f, "\n");
}

void CMFreeDiffusion::calculatePhi(CPoint &point, const CNodeTable &supDomain, double supDomainDiameter,
							bool calculateDerivatives, double *phi, double *phi_x, double *phi_y)
//calculates values of base functions of a given point
//	with given support nodes and support domain diameter
//if specified, the derivaties are also calculated
//the output tables should have space for at least supDomain.getNumElems() numbers
{
	const int MaxNumMonomials = 6;
	int numMonomials = (options.MLSdegree == 1 ? 3 : 6);

	//gsl matrices & vectors are static variables and thus only allocated on first call of calcualtePhi()
	//they are only freed when the whole program terminates
	static gsl_matrix *A = gsl_matrix_calloc(numMonomials, numMonomials);
	static gsl_vector *a = gsl_vector_alloc(numMonomials);
	double p[MaxNumMonomials], //\vec{p}(some node from supDom)
		p_point[MaxNumMonomials]; //\vec{p}(point)
	double *w = new double[supDomain.getNumElems()]; //MLS weights of all nodes

	//calculate p(point) and MLS weights of all nodes
	if (options.MLSdegree == 1)
		pmls<1>(point, p_point);
	else
		pmls<2>(point, p_point);
	for (int k = 0; k < supDomain.getNumElems(); k++)
		w[k] = Wmls(supDomain[k], point, supDomainDiameter);

	//construct A
	gsl_matrix_set_zero(A);
	for (int k = 0; k < supDomain.getNumElems(); k++)
	{
		if (options.MLSdegree == 1)
			pmls<1>(supDomain[k], p);
		else
			pmls<2>(supDomain[k], p);
		for (int i = 0; i < numMonomials; i++)
			for (int j = 0; j < numMonomials; j++)
				gsl_matrix_set(A, i, j, gsl_matrix_get(A, i, j) + w[k]*p[i]*p[j]);
	}
	//decompose A to (LU,p)
	int signum;
	static gsl_permutation *perm = gsl_permutation_alloc(numMonomials);
	gsl_linalg_LU_decomp(A, perm, &signum);

	//similar things for derivatives (which will only be allocated & evaluated if needed)
	static gsl_matrix *A_x = gsl_matrix_calloc(numMonomials, numMonomials), 
					*A_y = gsl_matrix_calloc(numMonomials, numMonomials); //derivatives of A & B
	static gsl_vector *gamma = gsl_vector_alloc(numMonomials);
	static double p_point_x[MaxNumMonomials], p_point_y[MaxNumMonomials]; //derivatives of \vec{p}(point)
	static double delta_x[MaxNumMonomials], delta_y[MaxNumMonomials]; //delta_x = p_x^T - g^T*A_x
	double *w_x = NULL, *w_y = NULL; //derivatives of MLS weights of nodes

	if (calculateDerivatives)
	{
		//calculate derivatives of MLS weights & of monomials
		w_x = new double[supDomain.getNumElems()];
		w_y = new double[supDomain.getNumElems()];
		if (options.MLSdegree == 1)
			pmls_xy<1>(point, p_point_x, p_point_y);
		else
			pmls_xy<2>(point, p_point_x, p_point_y);
		for (int k = 0; k < supDomain.getNumElems(); k++)
		{
			w_x[k] = Wmls_x(supDomain[k], point, supDomainDiameter);
			w_y[k] = Wmls_y(supDomain[k], point, supDomainDiameter);
		}

		//construct A_x, A_y
		gsl_matrix_set_zero(A_x);
		gsl_matrix_set_zero(A_y);
		for (int k = 0; k < supDomain.getNumElems(); k++)
		{
			if (options.MLSdegree == 1)
				pmls<1>(supDomain[k], p);
			else
				pmls<2>(supDomain[k], p);
			for (int i = 0; i < numMonomials; i++)
				for (int j = 0; j < numMonomials; j++)
				{
					gsl_matrix_set(A_x, i, j, gsl_matrix_get(A_x, i, j) + w_x[k]*p[i]*p[j]);
					gsl_matrix_set(A_y, i, j, gsl_matrix_get(A_y, i, j) + w_y[k]*p[i]*p[j]);
				}
		}

		//solve A*gamma = p(x)
		for (int i = 0; i < numMonomials; i++)
			gsl_vector_set(gamma, i, p_point[i]);
		gsl_linalg_LU_svx(A, perm, gamma);

		//calc delta_x = p_x^T - g^T*A_x
		for (int i = 0; i < numMonomials; i++)
		{
			delta_x[i] = p_point_x[i];
			delta_y[i] = p_point_y[i];
			for (int j = 0; j < numMonomials; j++)
			{
				delta_x[i] -= gsl_vector_get(gamma, j) * gsl_matrix_get(A_x, j, i);
				delta_y[i] -= gsl_vector_get(gamma, j) * gsl_matrix_get(A_y, j, i);
			}
		}
	}

	//evaluate phi, phi_x, phi_y
	for (int j = 0; j < supDomain.getNumElems(); j++)
	{
		//construct RHS and solve MLS system for a
		if (options.MLSdegree == 1)
			pmls<1>(supDomain[j], p);
		else
			pmls<2>(supDomain[j], p);
		for (int i = 0; i < numMonomials; i++)
			gsl_vector_set(a, i, w[j]*p[i]);
		gsl_linalg_LU_svx(A, perm, a);

		//mult. p_point*a to get phi
		phi[j] = 0;
		for (int i = 0; i < numMonomials; i++)
			phi[j] += p_point[i] * gsl_vector_get(a, i);

		//evaluate phi_x, phi_y
		if (calculateDerivatives)
		{
			phi_x[j] = 0;
			phi_y[j] = 0;
			for (int i = 0; i < numMonomials; i++)
			{
				phi_x[j] += delta_x[i]*gsl_vector_get(a, i) + gsl_vector_get(gamma, i)*w_x[j]*p[i];
				phi_y[j] += delta_y[i]*gsl_vector_get(a, i) + gsl_vector_get(gamma, i)*w_y[j]*p[i];
			}
		}
	}

	//delete temporary arrays
	delete w;
	if (calculateDerivatives)
	{
		delete w_x;
		delete w_y;
	}
}

}; //end namespace mfree
