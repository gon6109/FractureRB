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

using namespace vcg;
using namespace tri;

using namespace std;

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

	FractureDataGenerator::FractureDataGenerator(bool useDefaultSolver)
	{
		std::random_device seedGen;
		random = mt19937(seedGen());

		broadphase = new btDbvtBroadphase();
		collisionConfiguration = new btDefaultCollisionConfiguration();
		dispatcher = new btCollisionDispatcher(collisionConfiguration);
		if (useDefaultSolver)
			solver = new btSequentialImpulseConstraintSolver();
		else
			solver = new btMLCPSolver(new btDantzigSolver());
		dynamicsWorld = new btDiscreteDynamicsWorld(dispatcher, broadphase, solver, collisionConfiguration);
		btGImpactCollisionAlgorithm::registerAlgorithm(dispatcher);
	}

	FractureDataGenerator::~FractureDataGenerator()
	{
		for (int i = 0; i < dynamicsWorld->getNumCollisionObjects(); ++i) {
			string* name = (string*)(dynamicsWorld->getCollisionObjectArray()[i]->getUserPointer());
			if (name) {
				//				printf("%d (%d) -> %s\n",i,dynamicsWorld->getCollisionObjectArray()[i]->getUserIndex(),name->c_str() );
				delete name;
			}
			dynamicsWorld->getCollisionObjectArray()[i]->setUserPointer(NULL);
		}

		delete breakableRB;

		if (dynamicsWorld->getNumCollisionObjects())
			printf("\n%% deleting %d remaining objects from dynamics world", dynamicsWorld->getNumCollisionObjects());
		for (int i = dynamicsWorld->getNumCollisionObjects() - 1; i >= 0; --i) {
			btCollisionObject* obj = dynamicsWorld->getCollisionObjectArray()[i];
			if (obj) {
				//				printf("remove\n");
				dynamicsWorld->removeCollisionObject(obj);
				btRigidBody* rb = dynamic_cast<btRigidBody*>(obj);
				//				printf("motion state\n");
				if (rb && rb->getMotionState()) delete rb->getMotionState();
				//				printf("collision shape\n");
				if (obj->getCollisionShape()) delete obj->getCollisionShape();
				//				printf("del obj\n");
				delete obj;
			}
		}

		ColliderData::clearStoredData(false); // we've already delete the colliison shape for all the rigid bodies in the loop above

//		printf("del dynamics world and stuff\n");
		delete dynamicsWorld;
		delete solver;
		delete dispatcher;
		delete collisionConfiguration;
		delete broadphase;
	}

	int FractureDataGenerator::init(
		std::string inputFile,
		std::string paramFile,
		double splitImpulse,
		int useEstSIFs)
	{
		double minVoxelSize = DBL_MAX;

		vector<string> params;
		loadParamFile(paramFile, params);

		vector<btScalar> vertices;
		vector<int> indices;
		loadObj(inputFile, vertices, indices);

		btTriangleIndexVertexArray* triArray
			= new btTriangleIndexVertexArray(indices.size() / 3, indices.data(), 3 * sizeof(int), vertices.size(), vertices.data(), 3 * sizeof(btScalar));
		btGImpactMeshShape* shape = new btGImpactMeshShape(triArray);
		shape->updateBound();

		btRigidBody::btRigidBodyConstructionInfo rbci(
			0,
			new btDefaultMotionState(),
			shape,
			btVector3()
		);
		btRigidBody* rb = new btRigidBody(rbci);
		dynamicsWorld->addRigidBody(rb);

		if (params.size() > 9) {
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
			if (params.size() > 10) {
				check = sscanf(params[11].c_str(), "%d", &tmI); if (check == 1) maxCracks = tmI;
			}
			if (params.size() > 11) {
				check = sscanf(params[12].c_str(), "%d", &tmI); if (check == 1) method = (FRACTURE_METHOD)tmI;
			}
			if (params.size() > 12) {
				check = sscanf(params[13].c_str(), "%lf", &tmp); if (check == 1) smplDens = tmp;
			}

			if (minVoxelSize > voxelSize) minVoxelSize = voxelSize; // store the smallest voxel-size in the scene
			printf("\n%% ... building breakable rigid body:"
				"\n%% ... crack mesh size %.3lg / voxel size %.3lg = resolution ratio %.1lf"
				"\n%% ... Young's modulus %.3lg, Poisson's ratio %.3lf, density %.3lg"
				"\n%% ... tensile strength %.3lg, tensile toughness %.3lg, compressive factor %.3lf",
				meshSize, voxelSize, meshSize / voxelSize, E, nu, rho, Sc, Kc, cf
			);/**/
			breakableRB = new FractureRB(rb, false); //FractureRB class should not delete this rb since it's the importer's job to do so
			breakableRB->setOutputDir(outDir);
			if (useEstSIFs >= 0) breakableRB->setEstSIFsThreshold(useEstSIFs);
			breakableRB->initFractureSim(
				meshSize, voxelSize, E, nu, rho, Sc, Kc, cf, remesh, noSI, material, maxCracks, method, smplDens
			);
		}

		return 0;
	}

	int FractureDataGenerator::generate(
		std::ostream* out,
		int iterations,
		int gridNum,
		btVector3 gridMin,
		btVector3 gridMax,
		double dt,
		double impulseThreshold,
		double forceThreshold)
	{
		ostream& stressLog = *out;

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

		for (int i = 0; i < iterations; i++)
		{
			printf("iteration: %d\n", i);

			breakableRB->clearContacts();
			double totalImpulse = 0.0;

			int contactNum = random() % 5 + 1; // [1, 5]
			for (int l = 0; l < contactNum; l++)
			{
				auto& elems = breakableRB->getBEM()->getElems();
				auto& nodes = breakableRB->getBEM()->getNodes();
				int elemId = next(elems.begin(), random() % elems.size())->first;

				Eigen::Vector3d elemP0(nodes[elems[elemId][0]][0], nodes[elems[elemId][0]][1], nodes[elems[elemId][0]][2]);
				Eigen::Vector3d elemP1(nodes[elems[elemId][1]][0], nodes[elems[elemId][1]][1], nodes[elems[elemId][1]][2]);
				Eigen::Vector3d elemP2(nodes[elems[elemId][2]][0], nodes[elems[elemId][2]][1], nodes[elems[elemId][2]][2]);
				Eigen::Vector3d elemNormal = (elemP1 - elemP0).cross(elemP2 - elemP0).normalized();

				double t1 = random() % 100 / 100.0;
				double t2 = random() % (int)(100 - t1 * 100) / 100.0;
				Eigen::Vector3d pos = elemP0 + (elemP1 - elemP0) * t1 + (elemP2 - elemP0) * t2;

				Eigen::Vector3d velocity;
				do
				{
					velocity = Eigen::Vector3d(
						random() % 100 / 100.0,
						random() % 100 / 100.0,
						random() % 100 / 100.0
					);
					velocity.normalize();
					velocity *= random() % 500 / 10.0 + 0.01; // velocity [0.001, 50] m/s
				} while (velocity.dot(elemNormal) > 0);
				Eigen::Vector3d dir = velocity * (random() % 80) / 100.0;

				unsigned int contactElem;
				double r = std::abs(breakableRB->getCurvatureAndElementID(pos, contactElem));
				double e = 1.0 / breakableRB->getEovNuSq();
				double m = breakableRB->getRB()->getInvMass();

				double mass = (random() % 950 + 50); // mass [50, 1000] kg
				m += 1.0 / mass;

				r = 1.0 / r;
				e = 1.0 / e;
				m = 1.0 / m;

				double v = std::abs((-velocity).dot(dir.normalized()));
				double t_c = std::pow(m * m / (e * e * r * v), 0.2) * 2.8682656991953313822367;

				double m1 = breakableRB->getRB()->getMass(), m2 = mass;
				double j = std::abs(((1 + breakableRB->getRB()->getRestitution()) * (m1 * m2)
					/ (m1 + m2 + m1 * m2 * (random() % 100 / 1000.0)) * -velocity).dot(dir));

				dir *= -1;
				breakableRB->addContact(contactElem, dir, j, t_c);
				totalImpulse += j;
			}

			printf("fracture sim\n");
			int check = breakableRB->runFractureSim(dt, 0);

			vect3d_map sigma, tau;
			vect3d_map unused;
			breakableRB->getBEM()->computeInteriorStresses(sigma, tau, unused, stressEvalPos);

			printf("contact points: %d, total impulse: %lf, total force: %lf\n", contactNum, totalImpulse, breakableRB->getTotalContactForce());
			string crackObjFileName = "";
			if ((totalImpulse >= impulseThreshold || breakableRB->getTotalContactForce() >= forceThreshold)
				&& breakableRB->getTotalContactForce() > 0.0)
			{
				node_map nodes;
				elem_map elems;
				breakableRB->getBEM()->getLevelSet().mesh(nodes, elems);

				crackObjFileName = to_string(i) + ".obj";
				saveObj(crackObjFileName, nodes, elems);
				printf("save crack .obj file: %s\n", crackObjFileName.c_str());
			}

			for (int l = 0; l < sigma.size(); l++)
			{
				for (int m = 0; m < 3; m++)
				{
					stressLog << sigma[l][m] << ", ";
				}

				for (int m = 0; m < 3; m++)
				{
					stressLog << tau[l][m] << ", ";
				}
			}

			stressLog << crackObjFileName << endl;
		}
		return 0;
	}

	void FractureDataGenerator::setOutputDir(std::string outDir)
	{
	}

	int FractureDataGenerator::loadObj(std::string filename, std::vector<btScalar>& vertices, std::vector<int>& indices)
	{
		localReaderVCG::MyMesh mesh;
		vertices.clear();
		indices.clear();
		int r = io::Importer<localReaderVCG::MyMesh>::Open(mesh, filename.c_str());

		id_map vertexId; unsigned int k = NODE_BASE_INDEX;
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
		std::map<unsigned int, localReaderVCG::MyVertex> vertex_map;
		localReaderVCG::MyMesh::VertexIterator vi = vcg::tri::Allocator<localReaderVCG::MyMesh>::AddVertices(mesh, nodes.size());
		for (auto node : nodes)
		{
			vi->P()[0] = node.second[0];
			vi->P()[1] = node.second[1];
			vi->P()[2] = node.second[2];
			vertex_map[node.first] = *vi;
			vi++;
		}

		localReaderVCG::MyMesh::FaceIterator fi = vcg::tri::Allocator<localReaderVCG::MyMesh>::AddFaces(mesh, elems.size());
		for (auto elem : elems)
		{
			fi->V(0) = &(vertex_map[elem.second[0]]);
			fi->V(1) = &(vertex_map[elem.second[1]]);
			fi->V(2) = &(vertex_map[elem.second[2]]);
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

				return 0;
			}
		}
		in.close();

		return 0;
	}
}