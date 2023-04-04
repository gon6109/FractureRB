#include "FractureDataGenerator.h"
#include "FractureRB.h"
#include "ColliderData.h"
#include "FractureBEM.h"
#include "PostProcessor.h"
#include "FractureRB_config.h"
#include "VDBWrapper.h"

#include <BulletCollision/Gimpact/btGImpactCollisionAlgorithm.h>

#include <BulletDynamics/MLCPSolvers/btMLCPSolver.h>
#include <BulletDynamics/MLCPSolvers/btDantzigSolver.h>
#include <vcg/complex/complex.h>
#include <vcg/wrap/io_trimesh/import.h>
#include <vcg/wrap/io_trimesh/export_obj.h>

#define _SILENCE_EXPERIMENTAL_FILESYSTEM_DEPRECATION_WARNING
#include <experimental/filesystem>

using namespace vcg;
using namespace tri;

using namespace std;
namespace fs = std::experimental::filesystem;

namespace FractureSim
{
	namespace localReaderVCG {
		class MyVertex;
		class MyEdge;
		class MyFace;
		struct MyUsedTypes : public UsedTypes<Use<MyVertex>::AsVertexType, Use<MyEdge>::AsEdgeType, Use<MyFace>::AsFaceType> {};
		class MyVertex : public Vertex< MyUsedTypes,
			vertex::VFAdj,
			vertex::Coord3f,
			vertex::Normal3f,
			vertex::Mark,
			vertex::BitFlags  > {};
		class MyEdge : public Edge< MyUsedTypes> {};
		class MyFace : public Face< MyUsedTypes,
			face::VFAdj,
			face::VertexRef,
			face::BitFlags > {};
		class MyMesh : public vcg::tri::TriMesh<std::vector<MyVertex>, std::vector<MyFace> > {};
	}

	Eigen::Vector3d toEigen(const btVector3& in) {
		return Eigen::Vector3d(in[0], in[1], in[2]);
	}
	btVector3 fromEigen(const Eigen::Vector3d& in) {
		return btVector3(in[0], in[1], in[2]);
	}
	vdb::Vec3d toVDB(const btVector3& in) {
		return vdb::Vec3d(in[0], in[1], in[2]);
	}

	FractureDataGenerator::FractureDataGenerator()
	{
		std::random_device seedGen;
		random = mt19937(seedGen());
	}

	FractureDataGenerator::~FractureDataGenerator()
	{
		ColliderData::clearStoredData(false); // we've already delete the colliison shape for all the rigid bodies in the loop above
	}

	int FractureDataGenerator::generate(
		std::string inputFile,
		std::string paramFile,
		double splitImpulse,
		int useEstSIFs,
		int iterations,
		int gridNum,
		int crackGridNum,
		btVector3 gridMin,
		btVector3 gridMax,
		double dt,
		int start,
		int fileIdxOffset)
	{
		vector<string> params;
		loadParamFile(paramFile, params);

		vector<btScalar> vertices;
		vector<int> indices;
		loadObj(inputFile, vertices, indices);

		btVector3 gridSize = (gridMax - gridMin) / gridNum;
		vector<Eigen::Vector3d> stressEvalPos;
		{
			for (int z = 0; z < gridNum; z++)
			{
				for (int y = 0; y < gridNum; y++)
				{
					for (int x = 0; x < gridNum; x++)
					{
						Eigen::Vector3d evalPos =
							toEigen(btVector3(gridSize.x() * x, gridSize.y() * y, gridSize.z() * z) + gridSize / 2 + gridMin);
						stressEvalPos.push_back(evalPos);
					}
				}
			}
		}

		for (int i = start; i < iterations; i++)
		{
			printf("\niteration: %d\n", i + fileIdxOffset);

			string outputFileName = outDir + to_string(i + fileIdxOffset);
			if (fs::exists(outputFileName + "_stress.csv") && fs::exists(outputFileName + "_crack.csv"))
			{
				printf("exists crack & stress files.");
				continue;
			}

			initFractureRB(params, vertices, indices, splitImpulse, useEstSIFs);
			printf("\n");
			double totalImpulse = 0.0;
			std::vector<Eigen::Vector3d> contactPositions;
			std::vector<Eigen::Vector3d> contactDirs;
			std::vector<float> impulses;
			std::vector<float> durations;

			int contactNum = 1; // [1, 1]
			//for (int l = 0; l < contactNum; l++)
			{
				auto& elems = breakableRB->getBEM()->getElems();
				auto& nodes = breakableRB->getBEM()->getNodes();
				int elemId = next(elems.begin(), random() % elems.size())->first;
				//int elemId = next(elems.begin(), i % elems.size())->first;
				if (elemId == elems.size())
				{
					delete breakableRB;
					continue;
				}

				Eigen::Vector3d elemP0(nodes[elems[elemId][0]][0], nodes[elems[elemId][0]][1], nodes[elems[elemId][0]][2]);
				Eigen::Vector3d elemP1(nodes[elems[elemId][1]][0], nodes[elems[elemId][1]][1], nodes[elems[elemId][1]][2]);
				Eigen::Vector3d elemP2(nodes[elems[elemId][2]][0], nodes[elems[elemId][2]][1], nodes[elems[elemId][2]][2]);
				Eigen::Vector3d elemCenter = (elemP0 + elemP1 + elemP2) / 3;
				Eigen::Vector3d elemNormal = (elemP1 - elemP0).cross(elemP2 - elemP0).normalized();

				Eigen::Vector3d dir;
				//do
				//{
				//	velocity = Eigen::Vector3d(
				//		(random() % 200) / 100.0 - 1,
				//		(random() % 200) / 100.0 - 1,
				//		(random() % 200) / 100.0 - 1
				//	);
				//	velocity.normalize();
				//	velocity *= i / elems.size() % 5 / 5.0 * 0.9 + 0.1; // velocity [0.1, 1] m/s
				//} while (velocity.dot(elemNormal) > -cos(M_PI * 30.0f / 180.0f));

				int dirId = random() % 9;
				//int dirId = 4;
				Eigen::Matrix3d rot = (Eigen::AngleAxisd((dirId % 3 - 1) * 5.0f * M_PI / 180.0f, Eigen::Vector3d::UnitY())
					* Eigen::AngleAxisd((dirId / 3 - 1) * 5.0f * M_PI / 180.0f, Eigen::Vector3d::UnitZ())).matrix();

				//dir = rot * (-elemNormal);
				dir = rot * (-elemCenter);
				dir.normalize();
				//velocity *= i / elems.size() % 5 / 5.0 * 0.9 + 0.1; // velocity [0.1, 1] m/s

				double impulse = 100000 * dir.norm();
				//double duration = random() % (int)((dt - breakableRB->getBEM()->getTimeStep()) * 100000) / 100000.0 + breakableRB->getBEM()->getTimeStep();
				double duration = dt;

				breakableRB->addContact(elemId, dir, impulse, duration);
				printf("contact element id: %d(%lf,%lf,%lf), dir: (%lf,%lf,%lf), impulse: %lf, duration: %lf, mass: %lf\n",
					elemId,
					elemCenter.x(),
					elemCenter.y(),
					elemCenter.z(),
					dir.x(),
					dir.y(),
					dir.z(),
					impulse,
					duration,
					breakableRB->getRB()->getMass());
				contactPositions.push_back(elemCenter);
				contactDirs.push_back(dir);
				impulses.push_back(impulse);
				durations.push_back(duration);
				totalImpulse += impulse;
			}

			vect3d_map sigma, tau;
			vect3d_map unused;

			printf("fracture sim\n");
			int check = breakableRB->runFractureSim(dt, 0, [&]() {
				breakableRB->getBEM()->computeInteriorStresses(sigma, tau, unused, stressEvalPos);
				});

			if (breakableRB->getTotalContactForce() < 1e-6f)
			{
				printf("no fracture, ");
				printf("contact points: %d, total impulse: %lf, total force: %lf\n", contactNum, totalImpulse, breakableRB->getTotalContactForce());

				delete breakableRB;
				continue;
			}

			printf("contact points: %d, total impulse: %lf, total force: %lf\n", contactNum, totalImpulse, breakableRB->getTotalContactForce());
			breakableRB->getBEM()->writeVisualMesh(outputFileName, -1, false, false, true, true);
			breakableRB->getBEM()->getLevelSet().writeSdfCsv(outputFileName + "_crack.csv", toVDB(gridMin), toVDB(gridMax), crackGridNum);

			ofstream stressLog((outputFileName + "_stress.csv").c_str());
			for (int l = 0; l < sigma.size(); l++)
			{
				for (int m = 0; m < 3; m++)
				{
					stressLog << sigma[l][m] << ", ";
				}

				for (int m = 0; m < 3; m++)
				{
					stressLog << tau[l][m];
					if (m < 2)
						stressLog << ", ";
					else
						stressLog << endl;
				}
			}
			stressLog.close();

			ofstream contactLog((outputFileName + "_contact.csv").c_str());
			for (int l = 0; l < contactPositions.size(); l++)
			{
				contactLog << contactPositions[l].x() << ", " << contactPositions[l].y() << ", " << contactPositions[l].z() << ", ";
				contactLog << contactDirs[l].x() << ", " << contactDirs[l].y() << ", " << contactDirs[l].z() << ", ";
				contactLog << impulses[l] << ", " << durations[l] << endl;
			}
			contactLog.close();

			printf("save file: %s.obj, %s_stress.csv, %s_crack.csv\n\n", outputFileName.c_str(), outputFileName.c_str(), outputFileName.c_str());

			delete breakableRB;
		}
		return 0;
	}

	void FractureDataGenerator::setOutputDir(std::string outDir)
	{
		this->outDir = outDir;
	}

	int FractureDataGenerator::loadObj(std::string filename, std::vector<btScalar>& vertices, std::vector<int>& indices)
	{
		localReaderVCG::MyMesh mesh;
		vertices.clear();
		indices.clear();
		int r = io::Importer<localReaderVCG::MyMesh>::Open(mesh, filename.c_str());

		id_map vertexId; unsigned int k = 0;
		for (localReaderVCG::MyMesh::VertexIterator vi = mesh.vert.begin(); vi != mesh.vert.end(); ++vi) if (!vi->IsD()) {
			vertices.push_back(vi->P()[0]);
			vertices.push_back(vi->P()[1]);
			vertices.push_back(vi->P()[2]);
			vertexId[vi - mesh.vert.begin()] = k;
			++k;
		}
		for (localReaderVCG::MyMesh::FaceIterator fi = mesh.face.begin(); fi != mesh.face.end(); ++fi) if (!fi->IsD()) {
			indices.push_back(vertexId[tri::Index(mesh, fi->V(0))]);
			indices.push_back(vertexId[tri::Index(mesh, fi->V(1))]);
			indices.push_back(vertexId[tri::Index(mesh, fi->V(2))]);
		}

		return r;
	}

	int FractureDataGenerator::saveObj(std::string filename, node_map& nodes, elem_map& elems)
	{
		localReaderVCG::MyMesh mesh;

		unsigned int firstVert = mesh.vert.size();
		std::map<unsigned int, localReaderVCG::MyVertex*> vertex_map;
		localReaderVCG::MyMesh::VertexIterator vi = vcg::tri::Allocator<localReaderVCG::MyMesh>::AddVertices(mesh, nodes.size());
		for (auto node : nodes)
		{
			vi->P()[0] = node.second[0];
			vi->P()[1] = node.second[1];
			vi->P()[2] = node.second[2];
			vertex_map[node.first] = &*vi;
			vi++;
		}

		localReaderVCG::MyMesh::FaceIterator fi = vcg::tri::Allocator<localReaderVCG::MyMesh>::AddFaces(mesh, elems.size());
		for (auto elem : elems)
		{
			fi->V(0) = vertex_map[elem.second[0]];
			fi->V(1) = vertex_map[elem.second[1]];
			fi->V(2) = vertex_map[elem.second[2]];
			fi++;
		}

		vcg::tri::io::ExporterOBJ<localReaderVCG::MyMesh> writerObj;
		writerObj.Save(mesh, filename.c_str(), vcg::tri::io::Mask::IOM_NONE);

		return 0;
	}

	// string-trim helpers from http://www.cplusplus.com/faq/sequences/strings/trim/
	std::string& trim_right_inplace(std::string& s, const std::string& delimiters = " \f\n\r\t\v") {
		return s.erase(s.find_last_not_of(delimiters) + 1);
	}

	std::string& trim_left_inplace(std::string& s, const std::string& delimiters = " \f\n\r\t\v") {
		return s.erase(0, s.find_first_not_of(delimiters));
	}

	std::string& trim_inplace(std::string& s, const std::string& delimiters = " \f\n\r\t\v") {
		return trim_left_inplace(trim_right_inplace(s, delimiters), delimiters);
	}

	int FractureDataGenerator::loadParamFile(std::string filename, std::vector<std::string>& params)
	{
		params.clear();
		ifstream in(filename.c_str());
		if (!in.is_open())
			return -1;

		string line, token, name;
		bool first, comment;
		while (in.good()) {
			getline(in, line);
			if (line.empty() || line.at(0) == '%') { // comment -- ignore
			}
			else {
				stringstream tokenize(line);
				comment = false;
				while (tokenize.good() && !comment) {
					getline(tokenize, token, ',');
					trim_inplace(token);

					if (!token.empty()) if (token.at(0) == '%') comment = true; //allow commenting out rest of line
					if (!comment) params.push_back(token);
					//printf("added token \"%s\" to %s\n",token.c_str(), name.c_str());
				}

				in.close();
				return 0;
			}
		}
		in.close();

		return 0;
	}

	int FractureDataGenerator::initFractureRB(std::vector<std::string> params, std::vector<btScalar> vertices, std::vector<int> indices, double splitImpulse, int useEstSIFs)
	{
		double minVoxelSize = DBL_MAX;

		btTriangleIndexVertexArray* triArray
			= new btTriangleIndexVertexArray(indices.size() / 3, indices.data(), 3 * sizeof(int), vertices.size() / 3, vertices.data(), 3 * sizeof(btScalar));
		btGImpactMeshShape* shape = new btGImpactMeshShape(triArray);
		shape->updateBound();

		btRigidBody::btRigidBodyConstructionInfo rbci(
			0,
			new btDefaultMotionState(),
			shape,
			btVector3()
		);
		btRigidBody* rb = new btRigidBody(rbci);
		rb->setUserIndex(0);
		rb->setUserPointer(new string("main"));

		if (params.size() > 11) {
			double tmp, meshSize = 1.0, voxelSize = 0.1, E = 1.0, nu = 0.3, rho = 1.0, Sc = 1.0, Kc = 1.0, cf = 1.0, smplDens = 1.0;
			bool noSI = false; int remesh = 0, tmI, check, maxCracks = -1;
			string material = "default";
			FRACTURE_METHOD method = FULL_BEM;
			check = sscanf(params[0].c_str(), "%lf", &tmp); if (check == 1) meshSize = tmp;
			check = sscanf(params[1].c_str(), "%lf", &tmp); if (check == 1) voxelSize = tmp;
			check = sscanf(params[2].c_str(), "%lf", &tmp); if (check == 1) E = tmp;
			check = sscanf(params[3].c_str(), "%lf", &tmp); if (check == 1) nu = tmp;
			check = sscanf(params[4].c_str(), "%lf", &tmp); if (check == 1) rho = tmp;
			check = sscanf(params[5].c_str(), "%lf", &tmp); if (check == 1) Sc = tmp;
			check = sscanf(params[6].c_str(), "%lf", &tmp); if (check == 1) Kc = tmp;
			check = sscanf(params[7].c_str(), "%lf", &tmp); if (check == 1) cf = tmp;
			check = sscanf(params[8].c_str(), "%d", &tmI); if (check == 1) remesh = tmI;
			check = sscanf(params[9].c_str(), "%lf", &tmp); if (check == 1 && tmp > 0.0) noSI = true;
			material = params[10];
			// some optional params ...
			if (params.size() > 11) {
				check = sscanf(params[11].c_str(), "%d", &tmI); if (check == 1) maxCracks = tmI;
			}
			if (params.size() > 12) {
				check = sscanf(params[12].c_str(), "%d", &tmI); if (check == 1) method = (FRACTURE_METHOD)tmI;
			}
			if (params.size() > 13) {
				check = sscanf(params[13].c_str(), "%lf", &tmp); if (check == 1) smplDens = tmp;
			}

			if (minVoxelSize > voxelSize) minVoxelSize = voxelSize; // store the smallest voxel-size in the scene
			printf("\n%% ... building breakable rigid body:"
				"\n%% ... crack mesh size %.3lg / voxel size %.3lg = resolution ratio %.1lf"
				"\n%% ... Young's modulus %.3lg, Poisson's ratio %.3lf, density %.3lg"
				"\n%% ... tensile strength %.3lg, tensile toughness %.3lg, compressive factor %.3lf",
				meshSize, voxelSize, meshSize / voxelSize, E, nu, rho, Sc, Kc, cf
			);/**/
			breakableRB = new FractureRB(rb, true);
			breakableRB->setOutputDir(outDir);
			if (useEstSIFs >= 0) breakableRB->setEstSIFsThreshold(useEstSIFs);
			breakableRB->initFractureSim(
				meshSize, voxelSize, E, nu, rho, Sc, Kc, cf, remesh, noSI, material, maxCracks, method, smplDens
			);
		}

		return 0;
	}
}