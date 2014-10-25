#ifndef __MFREEDIFFUSION_H_INCLUDED__
#define __MFREEDIFFUSION_H_INCLUDED__

#include "base.h"
#include "nodes.h"
#include "sparsematrix.h"
#include <math.h>
#include <common_lib.h>

namespace mfree
{

class CAvgDistMesh
//average nodal distance class
//-evaluates avg nodal distances on regular mesh & stores them
//-interpolates avg nodal distances in any point
{
private:
	int size;			//mesh size in each dimension
	double **avgDist;	//2D table (array of pointers to 1D arrays)
	double maxAvgDist;	//max value in table avgDist
	double dx;			//distance between 2 mesh points in each direction (=1/(size-1))

	void allocate(int _size)
	{
		if (avgDist)
		{
			//delete old table
			for (int i = 0; i < size; i++)
				delete avgDist[i];
			delete avgDist;
		}

		//allocate new table
		size = _size;
		dx = 1.0/(size-1);
		avgDist = new double*[size];
		for (int i = 0; i < size; i++)
			avgDist[i] = new double[size];
	}

public:
	CAvgDistMesh()
	{
		size = -1;
		dx = CommonLib::NaN;
		avgDist = NULL;
		maxAvgDist = CommonLib::NaN;
	}
	
	virtual ~CAvgDistMesh()
	{
		if (avgDist)
		{
			for (int i = 0; i < size; i++)
				delete avgDist[i];
			delete avgDist;
		}
	}

	void initialize(const CKDTree &nodesTree, const MFreeOptions &options);
	//evaluates average nodal distance between nodes given in nodesTree on the regular mesh
	//not parallelism-aware: if you want only one process to evaluate distances,
	//	you should call initialize() just on that on, then call CAvgDistMesh::broadcast()

	void broadcast(int myRank);
	//broadcasts the whole table of evaluated avg distances to all processors

	double &operator()(const CPoint &point)
	//returns reference to avg distance in regular mesh point closest to given point
	{
		return avgDist[int(floor(point.x/dx+0.5))][int(floor(point.y/dx+0.5))];
	}

	double interpolateAvgDist(const CPoint &point) const
	//returns the interpolated value of average nodal distance in given point
	{
		double ixExact = point.x/dx; //horizontal "index" of point that does not necessarily exist in regular mesh
		int ixFloor, ixCeil; //horizontal indices of closest points left and right from given point, respectively
		if (fabs(ixExact - CommonLib::round(ixExact)) < 1e-10)
			ixFloor = ixCeil = CommonLib::round(ixExact);
		else {
			ixFloor = int(floor(ixExact));
			ixCeil = ixFloor+1;
		}

		double iyExact = point.y/dx; //same for vertical indices
		int iyFloor, iyCeil; 
		if (fabs(iyExact - CommonLib::round(iyExact)) < 1e-10)
			iyFloor = iyCeil = CommonLib::round(iyExact);
		else {
			iyFloor = int(floor(iyExact));
			iyCeil = iyFloor+1;
		}

		//interpolate the value from four closest mesh points 
		double t = (point.x-ixFloor*dx) / dx; //normalized distance from ixFloor to ixExact
		double u = (point.y-iyFloor*dx) / dx; //normalized distance from iyFloor to iyExact
        return (1-t)*(1-u)*avgDist[ixFloor][iyFloor]
				+ t*(1-u)*avgDist[ixCeil][iyFloor]
				+ (1-t)*u*avgDist[ixFloor][iyCeil]
				+ t*u*avgDist[ixCeil][iyCeil];
	}

	double getMaxAvgDist() const
	{
		return maxAvgDist;
	}

	double getMaxAvgDist(const CRect &rect) const
	//returns max average distance in given rectangle
	{
		//avg distance in a rectangle has maximum value in one of the following points:
		double max = 0;
		//1. all four corners of rect
		max = CommonLib::mmax(max, interpolateAvgDist(CPoint(rect.xMin, rect.yMin)));
		max = CommonLib::mmax(max, interpolateAvgDist(CPoint(rect.xMax, rect.yMin)));
		max = CommonLib::mmax(max, interpolateAvgDist(CPoint(rect.xMin, rect.yMax)));
		max = CommonLib::mmax(max, interpolateAvgDist(CPoint(rect.xMax, rect.yMax)));

		int ixMinCeil = (int)ceil(rect.xMin/dx), //horizontal index of left-most regular mesh point in rectangle
			ixMaxFloor = (int)floor(rect.xMax/dx),
			iyMinCeil = (int)ceil(rect.yMin/dx),
			iyMaxFloor = (int)floor(rect.yMax/dx);
		for (int ix = ixMinCeil; ix <= ixMaxFloor; ix++)
		{
			//2. intersections of vertical mesh lines with horizontal rect borders
			max = CommonLib::mmax(max, interpolateAvgDist(CPoint(dx*ix, rect.yMin)));
			max = CommonLib::mmax(max, interpolateAvgDist(CPoint(dx*ix, rect.yMax)));
			//3. mesh points within rect
			for (int iy = iyMinCeil; iy <= iyMaxFloor; iy++)
				max = CommonLib::mmax(max, interpolateAvgDist(CPoint(dx*ix, dx*iy)));
		}
		for (int iy = iyMinCeil; iy <= iyMaxFloor; iy++)
		{	
			//4. intersections of horizontal mesh lines with vertical rect borders
			max = CommonLib::mmax(max, interpolateAvgDist(CPoint(rect.xMin, dx*iy)));
			max = CommonLib::mmax(max, interpolateAvgDist(CPoint(rect.xMax, dx*iy)));
		}
		return max;
	}

	CAvgDistMesh &operator=(const CAvgDistMesh &other)
	{
		allocate(other.size);
		maxAvgDist = other.maxAvgDist;
		for (int ix = 0; ix < size; ix++)
			for (int iy = 0; iy < size; iy++)
				avgDist[ix][iy] = other.avgDist[ix][iy];
		return *this;
	}
}; //end class CAvgDistMesh

class CMFreeDiffusion
//class responsible for constructing the linear system that should solve the diffusion equation on unit square
//typical invocation sequence for a serial MFree program is:
//1.	createNodes()
//2.	exportNodes()
//3.	prepareNodeTree()
//4.	initAvgDistMesh()
//5.	constructSystem()
//6.	exportSystem()
{
protected: 
	MFreeOptions options;

	int totalNumNodes;  //number of nodes on ALL processors
	CRect subdomain;	//subdomain of this process;
						//this process is responsible for constructing rows of A and B corresponding to nodes in subdomain
	int myRank, numProcs; //rank of this process and total number of processes

	CNodeTable nodesTable;	//table of all nodes sorted by their index
	CKDTree nodesTree;		//kD-tree used for node queries

	CAvgDistMesh avgDistMesh;//regular mesh with pre-evaluated average distances

	sparse::CSparseMatrix A, B; //system matrices; system is A*x_new = B*x_old + fstar;
								//after parallel construction of system, each process will have its part of matrices A and B
	double *fstar;				//RHS of system
								//during parallel construction of system, each process will allocate space for full
								//vector fstar, but will only put values into components of fstar corresponding to
								//its own nodes; components of fstar corresponding to other nodes will be NaN

	double getSupDomain(const CPoint &point, CNodeTable &supDomain) const
	//finds nodes in the support domain of given point among all nodes in nodesTree
	//-calculates supDomainDiameter
	//-empties table supDomain
	//-puts all nodes at most supDomainDiameter/2 from point to table supDomain
	//-returns supDomainDiameter
	//nodes returned also have distances from point calculated
	{
		static double maxR = 0;
		if (options.constSupNodeCount()) {
			//constant num of nodes in sup domain
			supDomain.emptyTable();
			nodesTree.findKNearestNodes(point, options.nI+1, supDomain);
			int lastIn = supDomain.getNumElems()-2;
			double supDomainDiameter = (supDomain[lastIn].dist + supDomain[lastIn+1].dist);
			supDomain.deleteElems(lastIn+1, lastIn+1);
			maxR = mmax(maxR, supDomainDiameter);
			return supDomainDiameter;
		} else {
			//interpolated sup domain radius
			double supDomainDiameter = options.alphaSupport * avgDistMesh.interpolateAvgDist(point);
			supDomain.emptyTable();
			nodesTree.findNodesInCircle(point, supDomainDiameter, supDomain);
			maxR = mmax(maxR, supDomainDiameter);
			return supDomainDiameter;
		}
	}

	double getSupDomain(const CPoint &point, CNodeTable &supDomain, const CNodeTable &candidates) const
	//finds nodes in the support domain of given point among nodes given in table candidates
	//-calculates supDomainDiameter
	//-empties table supDomain
	//-puts all nodes at most supDomainDiameter/2 from point to table supDomain
	//-returns supDomainDiameter
	//nodes returned also have distances from point calculated
	{
		double supDomainDiameter = options.alphaSupport * avgDistMesh.interpolateAvgDist(point);
		supDomain.emptyTable();
		candidates.findNodesInCircle(point, supDomainDiameter, supDomain);
		return supDomainDiameter;
	}

	void calculatePhi(CPoint &point, const CNodeTable &supDomain, double supDomainDiameter,
						bool calculateDerivatives, double *phi, double *phi_x, double *phi_y);
	//calculates values of base functions of a given point
	//	with given support nodes and support domain diameter
	//if specified, the derivaties are also calculated
	//the output tables should have space for at least supDomain.getNumElems() numbers

public:
	CMFreeDiffusion(const MFreeOptions &_options, int _myRank, int _numProcs)
	{	
		assert(options.validate());
		options = _options;
		fstar = NULL;
		myRank = _myRank;
		numProcs = _numProcs;
	}

	virtual ~CMFreeDiffusion()
	{
		if (fstar)
			delete fstar;
	}

	int getNumNodes()
	{
		return nodesTable.getNumElems();
	}

	void importNodes(const CNodeTable &_nodesTable);
	//copies nodes from given table and sorts them by index

	void createNodes();
	//calls CNodeTable::createNodes
	//reads all parameters from options

	void prepareNodeTree();
	//builds kd tree from nodes stored in nodesTable

	void initAvgDistMesh()
	{
		avgDistMesh.initialize(nodesTree, options);
	}

	void distributeNodes();
	//distributes the nodes among processes
	//kd tree containing subdomain nodes is also prepared on all processes because building it can be overlapped with communication
	//for hierarchical distribution, avgDistMesh should be initialized before and kd tree built
	
	bool testTree() const;
	//tests tree for correctness
	//tests finding in various rects & circles
	//tests division to subdomains

	void constructSystem();
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

	void exportNodes(char *filename) const;
	//exports nodes to a struct named filename in a new Matlab script named s_filename.m
	//creates file, writes preamble (clear filename; global filename;) and node matrix
	
	void exportSystem(char *filename) const;
	//gathers A, B, fstar on root and exports them to an existing Matlab script named s_filename.m
	//also exports options

	void deleteAllButAvgDistMesh()
	{
		nodesTable.emptyTable();
		nodesTree.emptyTree();
	}
	void copyAvgDistMeshFrom(const CMFreeDiffusion &other)
	{
		avgDistMesh = other.avgDistMesh;
	}
}; //end class CMFreeDiffusion



}; //end namespace mfree


#endif
