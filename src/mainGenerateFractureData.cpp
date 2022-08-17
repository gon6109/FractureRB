#include <fstream>
#include <tclap/CmdLine.h>

#include "FractureDataGenerator.h"

using namespace FractureSim;
using namespace std;

int main(int argc, char* argv[]) {
	setbuf(stdout, NULL);
	try {
		printf("\n%% parsing command line: \n%% %s ...\n%%  ", argv[0]);
		for (int i = 1; i < argc; ++i) printf(" %s", argv[i]);

		TCLAP::CmdLine cmd("Rigid Body and Fracture simulation", ' ', "using HyENA 2.0, Bullet 2.83 and OpenVDB 2.2 backend");
		TCLAP::UnlabeledValueArg<string> inputFile("input-file", "obj format input object", true, "default", "input-file", cmd);
		TCLAP::UnlabeledValueArg<string> paramFile("param-file", "parameter file (CSV format)", true, "default", "param-file", cmd);
		TCLAP::ValueArg<string> outDir("o", "out", "output directory and file prefix", false, "", "out-dir", cmd);
		TCLAP::ValueArg<double> rbTimestep("d", "time-step", "rigid body time step size", false, 1.0 / 60, "time-step", cmd);
		TCLAP::ValueArg<int> iterations("i", "iterations", "number of iterations", false, 64, "iterations", cmd);
		TCLAP::ValueArg<string> gridMin("", "gridMin", "minimam locate of grid AABB", false, "-1,-1,-1", "gridMin", cmd);
		TCLAP::ValueArg<string> gridMax("", "gridMax", "maximam locate of grid AABB", false, "1,1,1", "gridMax", cmd);
		TCLAP::ValueArg<int> gridNum("n", "gridNum", "number of stress grid cells", false, 64, "gridNum", cmd);
		TCLAP::ValueArg<int> crackGridNum("c", "crackGridNum", "number of crack SDF grid cells", false, 256, "gridNum", cmd);
		TCLAP::ValueArg<double> impulseSplit("s", "split-impulse", "depth in units of the minimal voxel size of all input objects where impulse splitting occurs in the RB solver (default disabled)", false, -1.0, "depth", cmd);
		TCLAP::SwitchArg        useDefaultRBsolver("", "default-solver", "use Bullet's default impulse solver instead of Danzig MLCP", cmd);
		TCLAP::ValueArg<int>    estFromElems("e", "est-elems", "if a fracture simulation contains more than the specified number of elements switch to faster estimated SIF evaluation (default -1=never) ", false, -1, "est-elems", cmd);

		cmd.parse(argc, argv);

		FractureDataGenerator generator(useDefaultRBsolver.getValue());
		generator.setOutputDir(outDir.getValue()); //ToDo: make this a constructor param. of BulletWrapper (?)
		printf("\n%% reading object %s\n", inputFile.getValue().c_str());

		printf("\n%% timestepping the rigid body sim ...\n");
		btVector3 gridMinVec = btVector3(-1, -1, -1), gridMaxVec = btVector3(1, 1, 1);
		double x, y, z;
		if (3 == sscanf(gridMin.getValue().c_str(), "%lf,%lf,%lf", &x, &y, &z)) gridMinVec = btVector3(x, y, z);
		if (3 == sscanf(gridMax.getValue().c_str(), "%lf,%lf,%lf", &x, &y, &z)) gridMaxVec = btVector3(x, y, z);
		generator.generate(
			inputFile.getValue(),
			paramFile.getValue(),
			impulseSplit.getValue(),
			estFromElems.getValue(),
			iterations.getValue(),
			gridNum.getValue(),
			crackGridNum.getValue(),
			gridMinVec,
			gridMaxVec,
			rbTimestep.getValue()
		);

		printf("\n%% simulation done!\n");
	}
	catch (TCLAP::ArgException& e) {
		std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
		return -1;
	}
	return 0;
}
