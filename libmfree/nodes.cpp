#include "nodes.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <common_lib.h>

namespace mfree
{

double assureFromBoundary(double x, double minDist)
//returns the x' that is closest to x and: x' >= minDist, x' <= 1-minDist
{
	return CommonLib::mmax(minDist, CommonLib::mmin(1-minDist, x));
}

int CNodeTable::createOrReadNodes(const MFreeOptions &options)
//generates approximately suggestedNumNodes nodes using algorithm given in dissertation
//returns actual number of nodes generated
//nodes are randomly moved from inital positions (ie those on a regular mesh)
//	by at most dx*irregularity
//given seed is used for the random generator
//also sets initial conditions: 0 on boundaries, 1 inside of domain
{
	if (options.nodesFile == NULL)
		return createNodes(options);
	else
		return readNodes(options.nodesFile);
}

int CNodeTable::readNodes(char *filename)
{
	FILE *f = fopen(filename, "r");
	if (f == NULL)
	{
		std::cerr << "File not found: " << filename << std::endl;
		return -1;
	}
	emptyTable();
	char line[1024];
	while (fgets(line, sizeof line, f) != NULL)
	{
		double x, y, isBoundary, initialCond;
		sscanf(line, "%lf %lf %lf %lf", &x, &y, &isBoundary, &initialCond);
		appendElem(CNode(getNumElems(), x, y, initialCond, isBoundary > 0.1));
	}
	fclose(f);
	return getNumElems();
}

int CNodeTable::createNodes(const MFreeOptions &options)
{
	assert(options.suggestedNumNodes >= 9);
	assert(CommonLib::in(options.irregularity, 0.0, 0.5));
	assert(options.rHole >= 0);
	
	numNodesX = CommonLib::round(sqrt( (double)options.suggestedNumNodes / (1 - CommonLib::pi*options.rHole*options.rHole) ));
	double dx = 1.0/(numNodesX-1);
	double minDistFromBoundary = CommonLib::mmax(1-options.irregularity, 0.1)*dx;
	double rnd1, rnd2;

	emptyTable();
	reallocTable(numNodesX*numNodesX);
	srand(options.randSeed);

	//4 corners
	appendElem(CNode(getNumElems(), 0, 0, CommonLib::NaN, true));
	appendElem(CNode(getNumElems(), 0, 1, CommonLib::NaN, true));
	appendElem(CNode(getNumElems(), 1, 0, CommonLib::NaN, true));
	appendElem(CNode(getNumElems(), 1, 1, CommonLib::NaN, true));
	//lower & upper boundaries
	for (int ix = 1; ix < numNodesX-1; ix++)
	{
		rnd1 = CommonLib::randFromRange(-options.irregularity, options.irregularity);
		rnd2= CommonLib::randFromRange(-options.irregularity, options.irregularity);
		appendElem(CNode(getNumElems(), assureFromBoundary(dx*(ix + rnd1), minDistFromBoundary), 0, CommonLib::NaN, true));
		appendElem(CNode(getNumElems(), assureFromBoundary(dx*(ix + rnd2), minDistFromBoundary), 1, CommonLib::NaN, true));
	}
	//left & right boundaries
	for (int iy = 1; iy < numNodesX-1; iy++)
	{
		rnd1 = CommonLib::randFromRange(-options.irregularity, options.irregularity);
		rnd2 = CommonLib::randFromRange(-options.irregularity, options.irregularity);
		appendElem(CNode(getNumElems(), 0, assureFromBoundary(dx*(iy + rnd1), minDistFromBoundary), CommonLib::NaN, true));
		appendElem(CNode(getNumElems(), 1, assureFromBoundary(dx*(iy + rnd2), minDistFromBoundary), CommonLib::NaN, true));
	}
	//inner nodes
	for (int ix = 1; ix < numNodesX-1; ix++)
		for (int iy = 1; iy < numNodesX-1; iy++)
		{
			//put node if it is at least ~1% probability it will not be in the hole nor on hole boundary
			if (options.rHole == 0 || (CPoint(dx*ix, dx*iy) - options.xyHole).length()
						> options.rHole + minDistFromBoundary/1.4 - dx*options.irregularity*0.9)
			{
				CNode p;
				do {
					//keep moving it until it is out of hole+boundary
					rnd1 = CommonLib::randFromRange(-options.irregularity, options.irregularity);
					rnd2 = CommonLib::randFromRange(-options.irregularity, options.irregularity);
					p = CNode(getNumElems(), assureFromBoundary(dx*(ix + rnd1), minDistFromBoundary),
												assureFromBoundary(dx*(iy + rnd2), minDistFromBoundary), CommonLib::NaN, false);
				} while (options.rHole > 0 && (p - CNode(-1, options.xyHole.x, options.xyHole.y)).length()
													< options.rHole + minDistFromBoundary/1.4);
				appendElem(p);
			}
		}
	//hole boundary & vicinity
	if (options.rHole > 0)
	{
		const double densenAroundHole = 1.4; //nodes in vicinity of hole will be densened by this factor;
											//done for comparability of results with Femlab, which also densens nodes

		//put nodes on hole boundary the same way as on other boundaries
		int numNodesHole = int(densenAroundHole * ceil(2*CommonLib::pi*options.rHole / dx));
		double dphi = 2*CommonLib::pi / numNodesHole;
		for (int iphi = 0; iphi < numNodesHole; iphi++)
		{
			double angle = dphi*(iphi + CommonLib::randFromRange(-options.irregularity, options.irregularity));
			appendElem(CNode(getNumElems(), options.xyHole.x + options.rHole*cos(angle),
											options.xyHole.y + options.rHole*sin(angle), CommonLib::NaN, true));
		}

		//denser nodes around hole
		for (int i = 0; i < 5*options.suggestedNumNodes; i++)
		{
			//generate a random node in vicinity of hole
			double phi = CommonLib::randFromRange(0.0, 2*CommonLib::pi);
			double rPlus = CommonLib::randFromRange(minDistFromBoundary/densenAroundHole,
														(densenAroundHole-1)*options.rHole);
			double densityFactor = densenAroundHole - rPlus/options.rHole;
			CPoint potentialNode = options.xyHole + CPoint(cos(phi), sin(phi))*(options.rHole+rPlus);

			//only add the new node if it is not too close to an existing node
			double minDist = 1e+100;
			for (int j = 0; j < getNumElems(); j++)
				minDist = CommonLib::mmin(minDist, (potentialNode-access(j)).length());
			if (minDist > dx/densityFactor)
			{
				appendElem(CNode(getNumElems(), potentialNode.x, potentialNode.y, CommonLib::NaN, false));
				if (options.verboseLevel >= 6)
					std::cout << "Node in vicinity of hole added: " << potentialNode.x << "," << potentialNode.y << std::endl;
			}
		}
	}

	//set initial conditions: 0 on boundaries, 1 inside of domain
	for (int i = 0; i < getNumElems(); i++)
		access(i).u = (access(i).boundary ? 0 : 1);

	//return actual number of generated nodes
	return getNumElems();
}

int compareNodesByX(const void *pn1, const void *pn2)
{
	CNode *n1 = (CNode*)pn1;
	CNode *n2 = (CNode*)pn2;
	if (n1->x < n2->x)
		return -1;
	else if (n1->x == n2->x)
		return 0;
	else
		return 1;
}

int compareNodesByY(const void *pn1, const void *pn2)
{
	CNode *n1 = (CNode*)pn1;
	CNode *n2 = (CNode*)pn2;
	if (n1->y < n2->y)
		return -1;
	else if (n1->y == n2->y)
		return 0;
	else
		return 1;
}

int compareNodesByIndex(const void *pn1, const void *pn2)
{
	CNode *n1 = (CNode*)pn1;
	CNode *n2 = (CNode*)pn2;
	if (n1->index < n2->index)
		return -1;
	else if (n1->index == n2->index)
		return 0;
	else
		return 1;
}

int compareNodesByDist(const void *pn1, const void *pn2)
{
	CNode *n1 = (CNode*)pn1;
	CNode *n2 = (CNode*)pn2;
	if (n1->dist < n2->dist)
		return -1;
	else if (n1->dist == n2->dist)
		return 0;
	else
		return 1;
}

void CNodeTable::sortNodes(SortOrder sortOrder)
//sorts nodes by given key in ascending order using qsort
{
	if (sortOrder == sortOrderByX)
		qsort(data, getNumElems(), sizeof(CNode), compareNodesByX);
	else if (sortOrder == sortOrderByY)
		qsort(data, getNumElems(), sizeof(CNode), compareNodesByY);
	else if (sortOrder == sortOrderByIndex)
		qsort(data, getNumElems(), sizeof(CNode), compareNodesByIndex);
	else if (sortOrder == sortOrderByDist)
		qsort(data, getNumElems(), sizeof(CNode), compareNodesByDist);
	else
		throw new MFreeException(__LINE__, "CNodeTable::sortNodes: invalid sort order");
}

void CNodeTable::writeToMFile(FILE *f, char *varName) const
//exports nodes to matrix varName using Matlab script format
{
	assert(f);
	fprintf(f, "%s = [\n", varName);
	for (int i = 0; i < getNumElems(); i++)
	{
		accessConst(i).writeToFile(f);
		fprintf(f, ";\n");
	}
	fprintf(f, "];\n");
}

void CNodeTable::findNodesInCircle(const mfree::CPoint &center, double diameter, mfree::CNodeTable &res) const
//appends all nodes in the given circle to res
//(any nodes previously in res will also remain there)
//also writes distance from center to the appended nodes into their dist field
{
	double dist;
	for (int i = 0; i < numElems; i++)
	{
		dist = (center-accessConst(i)).length();
		if (dist <= diameter/2)
		{
			res.appendElem(accessConst(i));
			res[res.getNumElems()-1].dist = dist;
		}
	}
}

int CKDTree::buildKDTree(int numNodes, CNode *nodesSortedX, CNode *nodesSortedY, int currDepth)
//recursively builds tree containing numNodes nodes, which must be given in two arrays:
//one sorted by x, the other containing the same nodes but sorted by y
//current depth should be 0 on first call
//algorithm is as given on Page 92
//returns maximal tree depth in final tree
//the arrays may be reordered during tree construction
{
	depth = currDepth;
	numNodesInSubtree = numNodes;
	coverRect.xMin = nodesSortedX[0].x;
	coverRect.xMax = nodesSortedX[numNodes-1].x;
	coverRect.yMin = nodesSortedY[0].y;
	coverRect.yMax = nodesSortedY[numNodes-1].y;

	int subtreeDepth, totalDepth = 0;
	CNode *nodesSortedLeft = NULL, *nodesSortedRight = NULL;
	int inLeft = 0, inRight = 0;

	if (currDepth % 2 == 0)
	{
		//this (root of tree) = first node whose x is equal to middleNode.x
		int rootIndex = numNodes/2;
		while (rootIndex > 0 && nodesSortedX[rootIndex-1].x == nodesSortedX[rootIndex].x)
			rootIndex--;
		CNode::operator=(nodesSortedX[rootIndex]);

		//divide nodesSortedY on left & right parts
		//nodes < than middle.x go to left subtree
		//nodes >= than middle.x go to right subtree
		if (rootIndex > 0)
			nodesSortedLeft = new CNode[rootIndex];
		if (rootIndex < numNodes-1)
            nodesSortedRight = new CNode[numNodes-rootIndex-1];

		for (int i = 0; i < numNodes; i++)
		{
			if (nodesSortedY[i].index == index)
				continue;
			if (nodesSortedY[i].x < x)
				nodesSortedLeft[inLeft++] = nodesSortedY[i];
			else
				nodesSortedRight[inRight++] = nodesSortedY[i];
		}

		//check
		assert(inLeft == rootIndex);
		assert(inLeft+inRight+1 == numNodes);

		//recursively build subtrees
		if (inLeft > 0)
		{
			left = new CKDTree;
			subtreeDepth = left->buildKDTree(inLeft, nodesSortedX, nodesSortedLeft, currDepth+1);
			totalDepth = CommonLib::mmax(totalDepth, subtreeDepth+1);
		}
		if (inRight > 0)
		{
			right = new CKDTree;
			subtreeDepth = right->buildKDTree(inRight, nodesSortedX+rootIndex+1, nodesSortedRight, currDepth+1);
			totalDepth = CommonLib::mmax(totalDepth, subtreeDepth+1);
		}
	} else {
		//this (root of tree) = first node whose y is equal to middleNode.y
		int rootIndex = numNodes/2;
		while (rootIndex > 0 && nodesSortedY[rootIndex-1].y == nodesSortedY[rootIndex].y)
			rootIndex--;
		CNode::operator=(nodesSortedY[rootIndex]);

		//divide nodesSortedX on left & right parts
		//nodes < than middle.y go to left subtree
		//nodes >= than middle.y go to right subtree
		if (rootIndex > 0)
			nodesSortedLeft = new CNode[rootIndex];
		if (rootIndex < numNodes-1)
            nodesSortedRight = new CNode[numNodes-rootIndex-1];

		for (int i = 0; i < numNodes; i++)
		{
			if (nodesSortedX[i].index == index)
				continue;
			if (nodesSortedX[i].y < y)
				nodesSortedLeft[inLeft++] = nodesSortedX[i];
			else
				nodesSortedRight[inRight++] = nodesSortedX[i];
		}

		//check
		assert(inLeft == rootIndex);
		assert(inLeft+inRight+1 == numNodes);

		//recursively build subtrees
		if (inLeft > 0)
		{
			left = new CKDTree;
			subtreeDepth = left->buildKDTree(inLeft, nodesSortedLeft, nodesSortedY, currDepth+1);
			totalDepth = CommonLib::mmax(totalDepth, subtreeDepth+1);
		}
		if (inRight > 0)
		{
			right = new CKDTree;
			subtreeDepth = right->buildKDTree(inRight, nodesSortedRight, nodesSortedY+rootIndex+1, currDepth+1);
			totalDepth = CommonLib::mmax(totalDepth, subtreeDepth+1);
		}
	}

	if (nodesSortedLeft)
		delete nodesSortedLeft;
	if (nodesSortedRight)
		delete nodesSortedRight;

	return totalDepth;
}

int CKDTree::findNodesInRect(const CRect &rect, CNodeTable &res) const
//appends found nodes to res & returns no of found
//can also be used to convert to CNodeTable (simply pass CRect(0,0,1,1))
{
	int numFound = 0;
	if (isInRect(rect))
	{
		res.appendElem(*this);
		numFound++;
	}

	//if left subtree's coverRect intersects rect, also search in left subtree
	if (left && ((depth % 2 == 0 && rect.xMin <= x)
				|| (depth % 2 == 1 && rect.yMin <= y)) )
		numFound += left->findNodesInRect(rect, res);

	//similar for right subtree
	if (right && ((depth % 2 == 0 && rect.xMax >= x)
				|| (depth % 2 == 1 && rect.yMax >= y)) )
		numFound += right->findNodesInRect(rect, res);

	return numFound;
}

int CKDTree::findNodesInCircle(const CPoint &center, double diameter, CNodeTable &res) const
//appends all nodes from the tree that lie in given circle to the list res
//(any nodes previously in res will also remain there)
//returns number of nodes found
//appended nodes will also have distances from center calculated
{
	int numFound = 0;
	double distToThis = (*this-center).length();
	if (distToThis <= diameter/2)
	{
		res.appendElem(*this);
		res[res.getNumElems()-1].dist = distToThis;
		numFound++;
	}

	//if left subtree's coverRect intersects circle, also search in left subtree
	if (left && center.distToClosest(left->coverRect) <= diameter/2)
		numFound += left->findNodesInCircle(center, diameter, res);

	//similar for right subtree
	if (right && center.distToClosest(right->coverRect) <= diameter/2)
		numFound += right->findNodesInCircle(center, diameter, res);

	return numFound;
}

void CKDTree::findKNearestNodes(const CPoint &center, int k, CNodeTable &res) const
{
	static double diameter = 0.1;
	do {
		res.emptyTable();
		findNodesInCircle(center, diameter, res);
		diameter *= sqrt(2);
	} while (res.getNumElems() < k);
	diameter /= sqrt(res.getNumElems()/k);
	res.sortNodes(CNodeTable::sortOrderByDist);
	if (res.getNumElems() > k)
		res.deleteElems(k, res.getNumElems()-1);
}

CRect CKDTree::getSubdomainRect(int rank, int numProcs, CNodeTable &extraNodes) const
//function called internally by the public version of getSubdomainRect; not meant to be called directly
//returns subdomain of process with a given rank (out of numProcs processes) using hierarchical distribution
//recursive algorithm as given on Page 120
//extraNodes is a list of all nodes that are higher than *this in the tree but
//	may belong to the subdomain of proc rank
//when initially calling, pass empty table as extraNodes
{
	assert(rank >= 0 && rank < numProcs);
	if (numProcs == 1)
	{
		//if there is only 1 process, the whole subtree + any remaining extraNodes are his
		CRect subdomain = coverRect;
		for (int i = 0; i < extraNodes.getNumElems(); i++)
			subdomain = subdomain + extraNodes[i];
		return subdomain;
	} else {
		//filter extraNodes; only pass on those that are on the correct side of this->x or this->y
		int i = 0;
		while (i < extraNodes.getNumElems())
		{
			if (depth % 2 == 0)
			{
				if ((rank % 2 == 0 && extraNodes[i].x >= x) || (rank % 2 == 1 && extraNodes[i].x < x))
					extraNodes.deleteElems(i, i);
				else
					i++;
			} else {
				if ((rank % 2 == 0 && extraNodes[i].y >= y) || (rank % 2 == 1 && extraNodes[i].y < y))
					extraNodes.deleteElems(i, i);
				else
					i++;
			}
		}
		//if last bit is 1, return right subtree; else left subtree
		if (rank % 2 == 1)
		{
			assert(right);
			extraNodes.appendElem(*this);
			return right->getSubdomainRect(rank >> 1, numProcs/2, extraNodes);
		} else {
			assert(left);
			return left->getSubdomainRect(rank >> 1, numProcs/2, extraNodes);
		}
	}
}

}; //end namespace mfree
