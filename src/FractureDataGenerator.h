#pragma once

#include <vector>
#include <btBulletDynamicsCommon.h>

namespace FractureSim {
	class FractureRB;

	class FractureDataGenerator
	{
	public:
		FractureDataGenerator(bool useDefaultSolver = false);
		virtual ~FractureDataGenerator();

		int init(
			std::string inputFile,
			std::string paramFile,
			double splitImpulse = 3.0,
			int useEstSIFs = -1,
			double dt = 1.0 / 250.0,
			double impulseThreshold = DBL_EPSILON,
			double forceThreshold = DBL_EPSILON
		);

		int generate(
			std::ostream* out = NULL,
			int iterations = 0,
			int gridNum = 64,
			btVector3 gridMin = btVector3(-1, -1, -1),
			btVector3 gridMax = btVector3(1, 1, 1));

		void setOutputDir(std::string outDir);

	protected:
		btBroadphaseInterface* broadphase;
		btDefaultCollisionConfiguration* collisionConfiguration;
		btCollisionDispatcher* dispatcher;
		btSequentialImpulseConstraintSolver* solver;
		btDiscreteDynamicsWorld* dynamicsWorld;

		std::string outDir; // directory and file prefix for output files
		FractureRB* breakableRB;

		int loadObj(std::string filename, std::vector<btScalar>& vertices, std::vector<int>& indices);

		int loadParamFile(std::string filename, std::vector<std::string>& params);
	};
}