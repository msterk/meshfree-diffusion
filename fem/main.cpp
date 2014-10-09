#include <iostream>
#include <fstream>
#include <stdio.h>
#include <time.h>

#include "../libmfree/sparsematrix.h"

/*

fem - reads a FEM mesh from file and builds system matrices A and B for solution of the diffusion equation
		using Crank-Nicolson time integration scheme

nodes are stored in an array of CNode objects

elements are stored in an array of CElement objects

constructAB() actually builds the matrices

*/

class CNode
{
public:
	static int lastIndex; //index of previously generated node; used to give each node a unique index

	int index; //unique index
	double x, y;
	bool boundary;
	double value, //value of nodal parameter
			fFEM; //value of right hand side vector f

	friend std::istream &operator>>(std::istream &str, CNode &node);
};

std::istream &operator>>(std::istream &str, CNode &node)
{
	str >> node.x >> node.y;
	node.index = ++CNode::lastIndex;
	int bound;
	str >> bound;
	node.boundary = bound ? true : false;
	node.value = node.boundary ? 0 : 1;	//initial conditions: 0 on boundary, 1 inside
	node.fFEM = 0;
	return str;
}

class CElement
{
public:
	CNode *node[3]; //a triangular element is defined as pointers to three nodes

	void getPhi(double x, double y, double *phi)
	//value of 3 base functions of the element
	//does not check if (x,y) lies within the element
	{
		double area2 = 2*getArea();
		phi[0] = (node[1]->x * node[2]->y - node[2]->x * node[1]->y
				+ (node[1]->y - node[2]->y)*x + (node[2]->x - node[1]->x)*y) / area2;
		phi[1] = (node[2]->x * node[0]->y - node[0]->x * node[2]->y
				+ (node[2]->y - node[0]->y)*x + (node[0]->x - node[2]->x)*y) / area2;
		phi[2] = (node[0]->x * node[1]->y - node[1]->x * node[0]->y
				+ (node[0]->y - node[1]->y)*x + (node[1]->x - node[0]->x)*y) / area2;
	}
	double getArea()
	{
		return (node[1]->x*node[2]->y - node[1]->y*node[2]->x
				- (node[0]->x*node[2]->y - node[0]->y*node[2]->x)
				+ node[0]->x*node[1]->y - node[0]->y*node[1]->x) / 2;
	}

	friend std::istream &operator>>(std::istream &str, CElement &elem);
};

#define MAXNUMELEMS 400000

int CNode::lastIndex = -1;
CNode nodes[MAXNUMELEMS*3]; //global nodes table
CElement elems[MAXNUMELEMS]; //global elements table
int numNodes, numElems;

std::istream &operator>>(std::istream &str, CElement &elem)
{
	int n[3];
	str >> n[0] >> n[1] >> n[2];
	for (int i = 0; i < 3; i++)
		elem.node[i] = &nodes[n[i]-1];
	return str;
}

void readMesh(char *filename)
//reads mesh from file in format:
//  numNodes numElems foo foo
//  foo foo foo foo
//  node1x node1y node1boundary
//  node2x node2y node2boundary
//  ...
//  elem1node1index elem1node2index elem1node3index
//  ...
//(foo is data that is ignored; node indices start with 1 in input file, but with 0 in this program
{
	std::ifstream file;
	file.open(filename);
	
	int ibla;
	double dbla;
	file >> numNodes >> numElems >> ibla >> ibla;
	file >> dbla >> dbla >> ibla >> ibla;

	for (int i = 0; i < numNodes; i++)
		file >> nodes[i];
	for (int i = 0; i < numElems; i++)
		file >> elems[i];
}

sparse::CSparseMatrix A, B;
sparse::CSparseMatrixInput *AsparseInput, *BsparseInput;
int AspInSize, BspInSize;

void constructAB(double c, double dt)
{
	AspInSize = 0, BspInSize = 0;
	AsparseInput = new sparse::CSparseMatrixInput[numNodes*20];
	BsparseInput = new sparse::CSparseMatrixInput[numNodes*20];

	double x[3], y[3], //coordinates of nodes in an element
		dphix[3], dphiy[3]; //derivatives of the 3 base functions within element
	double area;
	double I_A[3][3], I_B[3][3]; //matrices C^L and K^L
	const double I_Ac[3][3] ={{2, 1, 1},
							{1, 2, 1},
							{1, 1, 2}};

	for (int iElem = 0; iElem < numElems; iElem++)
	{
		//for each element: calculate contribution of this element to its nodes
		for (int iNode = 0; iNode < 3; iNode++)
		{
			x[iNode] = elems[iElem].node[iNode]->x;
			y[iNode] = elems[iElem].node[iNode]->y;
		}
		area = elems[iElem].getArea();
		dphix[0] = y[1]-y[2]; dphix[1] = y[2]-y[0]; dphix[2] = y[0]-y[1];
		dphiy[0] = x[2]-x[1]; dphiy[1] = x[0]-x[2]; dphiy[2] = x[1]-x[0];

		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
			{
				I_A[i][j] = I_Ac[i][j] * area/12;
				I_B[i][j] = c*dt/2 * (dphix[i]*dphix[j] + dphiy[i]*dphiy[j]) /4/area;
			}

		for (int i = 0; i < 3; i++)
			if (!elems[iElem].node[i]->boundary)
				for (int j = 0; j < 3; j++)
				{
					//add the contributions to the input lists that will be used to construct sparse matrices A and B
					AsparseInput[AspInSize++] = sparse::CSparseMatrixInput(
						elems[iElem].node[i]->index, elems[iElem].node[i]->index,
						elems[iElem].node[j]->index, I_A[i][j] + I_B[i][j]);
					BsparseInput[BspInSize++] = sparse::CSparseMatrixInput(
						elems[iElem].node[i]->index, elems[iElem].node[i]->index, 
						elems[iElem].node[j]->index, I_A[i][j] - I_B[i][j]);
				}
	}

	//construct equations for boundary nodes
	for (int i = 0; i < numNodes; i++)
		if (nodes[i].boundary)
		{
			AsparseInput[AspInSize++] = sparse::CSparseMatrixInput(nodes[i].index, nodes[i].index, nodes[i].index, 1);
			BsparseInput[BspInSize++] = sparse::CSparseMatrixInput(nodes[i].index, nodes[i].index, nodes[i].index, 1);
		}

	//construct matrices from lists of contributions
	A.initialize(numNodes, numNodes, 30);
	B.initialize(numNodes, numNodes, 30);
	A.addInputs(AsparseInput, AspInSize);
	B.addInputs(BsparseInput, BspInSize);
	delete []AsparseInput;
	delete []BsparseInput;
}

int main(int argc, char **argv)
{
	if (argc == 1)
	{
		std::cerr << "Usage: " << argv[0] << " triangxFile\n" << std::endl
			<< "    reads mesh from file in format:" << std::endl
			<< "        numNodes numElems foo foo" << std::endl
			<< "        foo foo foo foo" << std::endl
			<< "        node1x node1y node1boundary" << std::endl
			<< "        node2x node2y node2boundary" << std::endl
			<< "        ..." << std::endl
			<< "        elem1node1index elem1node2index elem1node3index" << std::endl
			<< "        ..." << std::endl
			<< "    (foo is data that is ignored; node indices start with 1)" << std::endl << std::endl
			<< "    Then constructs system matrices A and B," << std::endl
			<< "    and outputs time required for construction" << std::endl
			<< "    (matrices themselves are not printed)";
		return -1;
	}
	readMesh(argv[1]);
	clock_t start = clock();
	double magic = 0; //something depending on matrix elements;
					//we calculate and print that to ensure that the compiler won't throw everything out during optimization
	int reps;
	double dt = 0.125, c = 0.001; //time-step and diffusivity constant
	for (reps = 0; clock()-start < CLOCKS_PER_SEC*5; reps++)
	{	
		constructAB(c, dt);
		magic += A(rand() % numNodes, rand() % numNodes) + B(rand() % numNodes, rand() % numNodes);
	}
	clock_t end = clock();
	double timeSetup = (end-start)/(double)CLOCKS_PER_SEC/reps;
	std::cout << "setup time " << timeSetup << ", magic number " << magic << std::endl;

	return 0;
}
