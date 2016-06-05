// File: cleScene.h
// Author: Court Cutting
// Date: December 11, 2010
// Purpose: This class handles interface between Sifakis-Teran corotated linear elasticity physics library
//		and glslTriangle surgical scene.

#ifndef __cleScene__
#define __cleScene__

#include "config.h"
#include "GraphicsUtils/boundingBox.h"
#include <vector>
#include <string>
#include "incision.h"
#include "GraphicsUtils/Vec3f.h"
#include <utility>
#include "GraphicsUtils/skinGraphics.h"

#ifdef WIN32
#include <libHTTPViewer/SERVER.h>
#else
#include <HTTPViewer/SERVER.h>
#endif
#if !defined(NO_PHYSICS) && !defined(USE_RPC_INTERFACE)
#include <EngineFrontend/CLElib.h>
#endif

#if !defined(NO_PHYSICS) &&  defined(USE_RPC_INTERFACE)
#include <cstdint>
#include <arpa/inet.h>
#include <thrift/protocol/TBinaryProtocol.h>
#include <thrift/server/TSimpleServer.h>
#include <thrift/transport/TServerSocket.h>
#include <thrift/transport/TBufferTransports.h>
#include "CLElib_RPC.h"
#include <boost/shared_ptr.hpp>
using namespace ::apache::thrift;
using namespace ::apache::thrift::protocol;
using namespace ::apache::thrift::transport;
using namespace ::apache::thrift::server;
using namespace boost;
using namespace  ::CLElib_RPC_nm;
#endif


// FORWARD DECLARATIONS
class staticTriangle;
class wxGraphics;
class surgicalActions;

class cleScene
{
public:
	cleScene(const Configuration* config);
#if !defined(NO_PHYSICS) && defined(USE_RPC_INTERFACE)
    void InitializeCLERPC( const std::string& server, const std::string& port);
#endif

	~cleScene(void);
	void setWxGraphics(wxGraphics *wxG) {_wxg=wxG;}
	void setSurgicalActions(surgicalActions *sa) {_surgAct=sa;}
	bool loadScene(const std::string& scene_data);
	inline skinGraphics* getElasticSkinGraphics() {return _sg;}
	bool exciseSubobject(int triangle);	// removes entire subobject connected to triangle
//	float* getMaterialCoordinate(int vertexNumber)	{return _materialCoords[vertexNumber]._v; }	// careful no error checking
//	int closestMaterialVertex(float (&position)[3]);
//	bool getMaterialSpatialMapping(const char *objectHit, bool materialToSpatial, const float (&xyzFrom)[3], const float (&nrmFrom)[3], float (&xyzTo)[3], float (&nrmTo)[3]);
	int addHook(int triangle, float (&uv)[2], float(&pos)[3]);
	void deleteHook(const int hook_id);
	void setHookPosition(const float(&location)[3],const int hook_id);
    int addSuture(int triangle0, float (&uv0)[2], int triangle1, float (&uv1)[2]);
    void deleteSuture(const int suture_id);

	bool makeIncision(trianglesUVW* tri,int nPoints,float (*positions)[3],float (*normals)[3],int TstartTriangleEdge,float TstartParam,bool Tout);
	// next routine is for debug - nuke later
	void cleUpdate();
	void getBoundingBox(float (&box)[6]);
	void scaleScene(float scaleFactor);
	void saveModifiedScene();
	void createFixedNodeFile();
	incision* getIncisionPointer() {return &_incis;}
	void setPhysicsPause(bool pause) { _pausePhysics=pause; }

    void reloadEngine();

private:
    const Configuration* _config;
	bool _pausePhysics;
	bool _physicsRunning;
	void createAndStartNewPhysicsModel();
	void setFixedGeometry();
	void recurseTriangleRemoval( int triangle);
	wxGraphics *_wxg;
	std::string _smdDirectory;
	skinGraphics *_sg;	// dynamic triangulated skin object
	trianglesUVW *_tuvw;
	std::list<staticTriangle*> _triList;
	surgicalActions *_surgAct;
	incision _incis;
	std::string _triName;
	// next variables are to CLElib
	std::vector<float> _fixedVerts,_fixedmeshVerts,_collisionVerts; // _cleVerts,
	std::vector<int> _fixedPoints,_fixedSegments,_fixedTriangles,_fixedmeshTriangles,_collisionTris;
    typedef std::pair< std::vector<float>, std::vector<float> > REGION;
    std::vector< REGION > _fixedRegions;
    std::vector< std::vector<float> > _muscleVerts;
    std::vector< std::vector<int> > _muscleTris;
    std::vector< std::vector<float> > _muscleFibers;
    std::vector< float > _muscleMaxStress;
    float _hookStiffness, _sutureStiffness;
    float _poissons_ratio, _youngs_modulus;
    int _refinement;
    float _refinement_ratio;
    int _frame; // current physics frame;

	float _cleNodeSpacing;
#if !defined(NO_PHYSICS) && !defined(USE_RPC_INTERFACE)
	CLElib _cle;
#endif
#if !defined(NO_PHYSICS) && defined(USE_RPC_INTERFACE)
    bool rpc_initialized;
    CLE_Session RPC_Session;//(session);
    shared_ptr<TSocket> socket;//(new TSocket("localhost", 9090));
    shared_ptr<TBufferedTransport> transport;//(new TBufferedTransport(socket));
    shared_ptr<TProtocol> protocol;//(new TBinaryProtocol(transport));
    shared_ptr<CLElib_RPCClient> client;//(protocol);
#endif   

	boundingBox<float> _bBox;

};
#endif	// #ifndef __cleScene__
