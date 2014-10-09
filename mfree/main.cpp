#include "../libmfree/libmfree.h"
#include <time.h>
#include <iostream>
#include <mpi.h>
#include <stdio.h>

#ifdef WIN32
	#include <windows.h>
#endif

/*

mfree - creates nodes and builds system matrices A and B and RHS vector f for solution of the diffusion equation
		using Crank-Nicolson time integration scheme

		exports everything to a Matlab script

*/


mfree::MFreeOptions options;

int main(int argc, char **argv)
{
	#ifdef WIN32
		SetPriorityClass( GetCurrentProcess(), IDLE_PRIORITY_CLASS );
	#endif 

	//MPI init
	int numProcs, myRank, namelen;
	char myName[200];
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&numProcs);
    MPI_Comm_rank(MPI_COMM_WORLD,&myRank);
    MPI_Get_processor_name(myName,&namelen);

	if (argc == 1)
	{
		//if no parameters given, just print help & exit
		if (myRank == 0)
		{
			std::cerr << "Usage: mfree outFile options\n\n"
				<< "Creates nodes, builds system matrices A,B and RHS vector f for MLPG1\n"
				<< "solution of diffusion equation; exports everything to Matlab script\n\n"
				<< "Options:\n"
				<< "n=I: suggestedNumNodes, default = " << options.suggestedNumNodes << std::endl
				<< "f=FILENAME: file to read nodes from, only allowed if n is not given" << std::endl
				<< "h=F: rHole, default = " << options.rHole << std::endl
				<< "x=F: xHole, default = " << options.xyHole.x << std::endl
				<< "y=F: yHole, default = " << options.xyHole.y << std::endl
				<< "c=F: c (diff. coef.), default = " << options.c << std::endl
				<< "i=F: irregularity, default = " << options.irregularity << std::endl
				<< "t=F: time step, default = " << options.timeStep << std::endl
				<< "s=F: alphaSupport (if float) or nI (if integer), default = " << options.nI << std::endl
				<< "q=F: alphaQuad (if alphaSupport is used), default = " << options.alphaQuad << std::endl
				<< "     betaQuad (if nI is used), default = " << options.betaQuad << std::endl
				<< "d=I: MLSdegree, default = " << options.MLSdegree << std::endl
				<< "g=I: Gaussian quadrature degree, default = " << options.quadDegree << std::endl
				<< "p=B: pre-find support of quad domain (t/f), default = "
						<< (options.preFindQuadSupport ? 't' : 'f') << std::endl
				<< "v=I: verbose level, default = " << options.verboseLevel << std::endl
				<< "r=I: randSeed, default = " << options.randSeed << std::endl
				<< "l=I: no. of repetitions (if >0) or seconds (if <0) to run, default = "
						<< options.numRepetitions << std::endl
				<< "z=I: dataDistrType (0=hiear., 1=1D, 2=all-to-all), default = "
						<< (options.dataDistrType == mfree::MFreeOptions::DataDistrHierarchical ? 0
								: (options.dataDistrType == mfree::MFreeOptions::DataDistr1D ? 1 : 2))
						<< std::endl;
		}
		MPI_Finalize();
		return -1;
	}

	if (myRank == 0)
	{
		//the root process parses the parameters
		for (int op = 2; op < argc; op++)
			switch(tolower(argv[op][0]))
			{
			case 'n':
				options.suggestedNumNodes = atoi(argv[op]+2);
				break;
			case 'f':
				options.nodesFile = new char[strlen(argv[op]+2)];
				strncpy(options.nodesFile, argv[op]+2, strlen(argv[op]+2)+1);
				break;
			case 'd':
				options.MLSdegree = atoi(argv[op]+2);
				break;
			case 'g':
				options.quadDegree = atoi(argv[op]+2);
				break;
			case 'v':
				options.verboseLevel = atoi(argv[op]+2);
				break;
			case 'h':
				options.rHole = atof(argv[op]+2);
				break;
			case 'x':
				options.xyHole.x = atof(argv[op]+2);
				break;
			case 'y':
				options.xyHole.y = atof(argv[op]+2);
				break;
			case 'c':
				options.c = atof(argv[op]+2);
				break;
			case 'i':
				options.irregularity = atof(argv[op]+2);
				break;
			case 't':
				options.timeStep = atof(argv[op]+2);
				break;
			case 's':
				if (strstr(argv[op]+2, ".") == NULL) {
					options.alphaSupport = NaN;
					options.nI = atoi(argv[op]+2);
				} else {
					options.alphaSupport = atof(argv[op]+2);
					options.nI = 0;
				}
				break;
			case 'q':
				options.alphaQuad = atof(argv[op]+2);
				options.betaQuad = options.alphaQuad;
				break;
			case 'r':
				options.randSeed = atoi(argv[op]+2);
				break;
			case 'p':
				options.preFindQuadSupport = (argv[op][2] == 't' ? true : false);
				break;
			case 'l':
				options.numRepetitions = atoi(argv[op]+2);
				break;
			case 'z':
				options.dataDistrType = (argv[op][2] == '0' ? mfree::MFreeOptions::DataDistrHierarchical
					: (argv[op][2] == '1' ? mfree::MFreeOptions::DataDistr1D : mfree::MFreeOptions::DataDistrUnstructured));
				break;
			default:
				std::cerr << "Invalid option.\n";
				break;
			}
		if (numProcs == 1 && options.dataDistrType != mfree::MFreeOptions::DataDistrHierarchical)
		{
			std::cerr << "Running on just 1 process requires hierarchical distribution!\n";
			MPI_Finalize();
			return -1;
		}
		if (options.preFindQuadSupport && options.dataDistrType != mfree::MFreeOptions::DataDistrHierarchical)
		{
			std::cerr << "Pre-find support of quad domain requires hierarchical distribution!\n";
			MPI_Finalize();
			return -1;
		}
	}

	CommonLib::commBroadcast(options); 
	if (options.verboseLevel >= 2)
	{
		fprintf(stderr,"Process %d on %s\n", myRank, myName);
		fflush(stderr);
	}

	mfree::CMFreeDiffusion *sampleMfreeDiffusion = NULL;
	int numNodes;
	if (!options.constSupNodeCount()) {
		//sampleMfreeDiffusion is created on each process and its avgDistMesh is initialized, which
		//	will later be copied into every created CMFreeDiffusion object (all CMFreeDiffusion objects
		//	will have the same nodes because they will be generated with the same random seed)
		//this is actually cheating (slightly), but I didn't have time to implement distribution
		//	of avgDistMesh from root to other processes
		sampleMfreeDiffusion = new mfree::CMFreeDiffusion(options, myRank, numProcs);
		sampleMfreeDiffusion->createNodes();
		sampleMfreeDiffusion->prepareNodeTree();
		sampleMfreeDiffusion->initAvgDistMesh();
		sampleMfreeDiffusion->deleteAllButAvgDistMesh();
	}

	MPI_Barrier(MPI_COMM_WORLD);
	double start = MPI_Wtime();
	int reps;
	for (reps = 0; (options.numRepetitions < 0 && (reps < 2 || MPI_Wtime()-start < (-options.numRepetitions)))
			|| (options.numRepetitions > 0 && reps < options.numRepetitions); reps++)
	{	
		mfree::CMFreeDiffusion mfreeDiffusion(options, myRank, numProcs);
		if (!options.constSupNodeCount())
			mfreeDiffusion.copyAvgDistMeshFrom(*sampleMfreeDiffusion); //the cheating part

		if (myRank == 0)
			mfreeDiffusion.createNodes();

		if (numProcs == 1)
			mfreeDiffusion.prepareNodeTree();
		else {
			if (options.dataDistrType == mfree::MFreeOptions::DataDistrHierarchical)
				mfreeDiffusion.prepareNodeTree();
			mfreeDiffusion.distributeNodes();
		}

		numNodes = mfreeDiffusion.getNumNodes();

		mfreeDiffusion.constructSystem();
		MPI_Barrier(MPI_COMM_WORLD);
		if (myRank == 0)
		{
			mfreeDiffusion.exportNodes(argv[1]);
			mfreeDiffusion.exportSystem(argv[1]);
			//exportSystem se sesuje z 1d distr., 400 nodi in 32 procesi!
		}
	}
	double end = MPI_Wtime();

	if (sampleMfreeDiffusion != NULL)
		delete sampleMfreeDiffusion;

	if (myRank == 0)
	{
		double timeSetup = (end-start)/reps;
		printf("Time %15f sec = %9f per n.(N=%10d, p=%10d)\n",
				timeSetup, timeSetup/numNodes, numNodes, numProcs);
	}

//	mfreeDiffusion.testTree();


	//append command: timeSetup = ... to end of file so that timing can be read into matlab
/*	char filename[1000];
	sprintf(filename, "s_%s.m", argv[1]);
	FILE *f = fopen(filename, "a+");
	fprintf(f, "timeSetup = %.20lf;\n", timeSetup);
	fclose(f);
*/
    MPI_Finalize();
	return 0;
}
