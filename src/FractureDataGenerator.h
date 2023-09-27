#pragma once

#include <random>
#include <vector>
#include <btBulletDynamicsCommon.h>
#include "types.h"

namespace FractureSim {
	class FractureRB;

	class FractureDataGenerator
	{
	public:
		FractureDataGenerator();
		virtual ~FractureDataGenerator();

		int generate(
			std::string inputFile,
			std::string paramFile,
			double splitImpulse = 3.0,
			int useEstSIFs = -1,
			int iterations = 0,
			int gridNum = 64,
			int crackGridNum = 256,
			btVector3 gridMin = btVector3(-1, -1, -1),
			btVector3 gridMax = btVector3(1, 1, 1),
			double dt = 1.0 / 250.0,
			int start = 0,
			int fileIdxOffset = 0,
			bool isOutputStress = true);

		void setOutputDir(std::string outDir);

	protected:
		std::mt19937 random;

		std::string outDir; // directory and file prefix for output files
		FractureRB* breakableRB;

		int loadObj(std::string filename, std::vector<btScalar>& vertices, std::vector<int>& indices);
		int saveObj(std::string filename, node_map& nodes, elem_map& elems);

		int loadParamFile(std::string filename, std::vector<std::string>& params);

		int initFractureRB(
			std::vector<std::string> params,
			std::vector<btScalar> vertices,
			std::vector<int> indices,
			double splitImpulse = 3.0,
			int useEstSIFs = -1);
	};
}