// File: cleScene.h
// Author: Court Cutting
// Date: December 11, 2010
// Purpose: This class handles interface between Sifakis-Teran corotated linear elasticity physics library
//		and surgical scene.
#include <exception>
#include "config.h"
#include "cleScene.h"
#include <map>
#include <algorithm>
#include <string.h>
#include <cassert>
#include "GraphicsUtils/staticTriangle.h"
#include "GraphicsUtils/Vec3f.h"
#include "surgicalActions.h"
#include "GraphicsUtils/lines.h"
#include "GraphicsUtils/wxGraphics.h"
#include <fstream>
#include <iostream>
//#include "json_spirit.h"	// for json_spirit.lib version
#include "json_spirit_reader_template.h"	// these 2 lines for header only version
#include "json_spirit_writer_template.h"

#ifndef JSON_SPIRIT_MVALUE_ENABLED
#define JSON_SPIRIT_MVALUE_ENABLED
#endif

#if !defined(NO_PHYSICS) &&  defined(USE_RPC_INTERFACE)
#include <thrift/protocol/TBinaryProtocol.h>
#include <thrift/protocol/TJSONProtocol.h>
#include <thrift/server/TSimpleServer.h>
#include <thrift/transport/TSocket.h>
#include <thrift/transport/TBufferTransports.h>

using namespace ::apache::thrift;
using namespace ::apache::thrift::protocol;
using namespace ::apache::thrift::transport;
using namespace ::apache::thrift::server;
using namespace boost;
using namespace  ::CLElib_RPC_nm;
#endif


namespace PhysBAM{
	int PhysBAM_number_of_threads=12;
	int PhysBAM_random_seed=0;
}



cleScene::cleScene(const Configuration* config): _config(config), _wxg(NULL),_sg(NULL),_tuvw(NULL),_surgAct(NULL),_physicsRunning(false),_pausePhysics(false)
{
	_bBox.Empty_Box();

    #if !defined(NO_PHYSICS) &&  defined(USE_RPC_INTERFACE)
    rpc_initialized = false;
    #endif

    _frame = -1;

    _hookStiffness = 1e4;
    _sutureStiffness = 1e5;
    _refinement_ratio = 3.0f;
    _poissons_ratio = 0.35f;
    _youngs_modulus = 3e4f;

}



void cleScene::reloadEngine(){
    if( _physicsRunning ){
        std::cout << "Cleaning up old scene..." << std::endl;
        
#if !defined(NO_PHYSICS) && defined(USE_RPC_INTERFACE)
        client->Destroy_Model(RPC_Session);	// no need to destroy a non-existent model
#elif !defined(NO_PHYSICS) && !defined(USE_RPC_INTERFACE)
        _cle.Destroy_Model();
#endif
        _physicsRunning=false;               
            
    }
}


#if !defined(NO_PHYSICS) && defined(USE_RPC_INTERFACE)
void cleScene::InitializeCLERPC( const std::string& server, const std::string& port)
{
    std::cout << "Initializing RPC connection to simulation server..." << std::endl;
    socket = shared_ptr<TSocket>(new TSocket(server,atoi(port.c_str())));
    transport = shared_ptr<TBufferedTransport>(new TBufferedTransport(socket));
    protocol = shared_ptr<TProtocol>(new TJSONProtocol(transport));
    client = shared_ptr<CLElib_RPCClient>(new CLElib_RPCClient(protocol));
    transport->open();
    std::cout << "Transport open, creating session..." << std::endl;
    client->Initialize_Session(RPC_Session);
    std::cout << "Session Ready." << std::endl;

    rpc_initialized = true;

}
#endif




cleScene::~cleScene(void)
{
#if !defined(NO_PHYSICS) && defined(USE_RPC_INTERFACE)
    if( rpc_initialized)
        {
            client->Shutdown_Session(RPC_Session);
            transport->close();
        }
#endif
    if(_sg)
        delete _sg;
}

bool cleScene::loadScene(const std::string& scene_data)
{
    if(_sg) 
        delete _sg;
    _sg = new skinGraphics();

    _hookStiffness = 1e4;
    _sutureStiffness = 1e5;
    _refinement_ratio = 3.0f;
    _poissons_ratio = 0.35f;
    _youngs_modulus = 3e4f;

    _wxg->clear();
	_physicsRunning = false;
	bool fixed;	// vs. dynamic
    std::string path;
	    json_spirit::mValue value;
    json_spirit::read_string( scene_data, value );	// header only version
	const json_spirit::mObject& smdObj = value.get_obj();
	//is.close();
	staticTriangle *tr;
//	std::string Obj(smdDirectory),Shader(smdDirectory),unparsedLine;
//	std::vector<std::string> parsedLine;
	std::map<int,GLuint> txMap;
    std::string nrm,tex;
	// get texture files first
	json_spirit::mObject::const_iterator oit,suboit,suboit2;
	if((oit=smdObj.find("textureFiles"))==smdObj.end())	{
		throw std::runtime_error("No texture files in scene file-");
	}
	else	{
		for(suboit=oit->second.get_obj().begin(); suboit!=oit->second.get_obj().end(); ++suboit)	{
			//path = _smdDirectory + suboit->first;
			int txNum = suboit->second.get_int();
            std::string texPath = suboit->first;
			//txMap.insert(std::make_pair(txNum,txNow));
            if(_surgAct->getServerMode()>0) {
                _surgAct->getServerHandle().RegisterTexture( std::to_string((long long)(txNum)), texPath );
                //tex.close();
            }
		}
	}
	if((oit=smdObj.find("fixedObjects"))!=smdObj.end())	{
		for(suboit=oit->second.get_obj().begin(); suboit!=oit->second.get_obj().end(); ++suboit)	{
			path = suboit->first;
			std::map<int,GLuint>::iterator tit;
			std::string nrm,tex;
			std::string nrmMap,texMap;
			for(suboit2=suboit->second.get_obj().begin(); suboit2!=suboit->second.get_obj().end(); ++suboit2)	{
				if(suboit2->first=="textureMap")	{
					texMap = std::to_string((unsigned long long)suboit2->second.get_int());
				}
				else if(suboit2->first=="normalMap")	{
                    nrmMap = std::to_string((unsigned long long)suboit2->second.get_int());					
				}
				else {
					throw std::runtime_error("Incorrect fixed triangle section in .smd input file-");
				}
			}
			if(texMap.empty() || nrmMap.empty()) {
				throw std::runtime_error("Missing texture or normal map in fixed triangle section in .smd input file-");
			}
            std::string obj_data = _surgAct->loadData( path, "STATIC" );
			if((tr=_wxg->loadStaticObjFile(obj_data.c_str(),NULL,NULL,true))==NULL)
			{
				throw std::runtime_error("Unable to load fixed triangle .obj input file-");
			}
			_triList.push_back(tr);
			path = suboit->first;
			size_t pos = path.rfind(".obj");
			path.erase(pos);
			// this is a staticTriangle, not elastic
			tr->setName(path.c_str());

                     if(_surgAct->getServerMode()>0) {
                         std::vector<int> tris;
                         std::vector<float> verts;
                         std::vector<float> uvs;
                         tr->getSurfaceTriangles(tris);
                         tr->getSurfaceVertices(verts);
                         tr->getSurfaceUVs(uvs);
                         std::vector<double> d_verts( verts.begin(), verts.end());
                         std::vector<double> d_uvs( uvs.begin(), uvs.end());
                         _surgAct->getServerHandle().AddStaticObject( tris,d_verts,d_uvs,texMap,nrmMap);
                    }
		}
	}
	if((oit=smdObj.find("collisionObject"))!=smdObj.end())	{
		path = oit->second.get_str();
		fixed = true;
		staticTriangle tempTr;
		// NATHAN - just need obj reader. Don't want to visualize it as already loaded that fixed .obj file.
        std::string obj_data = _surgAct->loadData( path, "STATIC" );
		if(tempTr.readObjFile(obj_data.c_str(),true)>0)
		{
			throw std::runtime_error("Unable to load collision triangle .obj input file-");
		}
		tempTr.getTriangulatedSurface(_collisionVerts,_collisionTris);
	}
	if((oit=smdObj.find("dynamicObjects"))!=smdObj.end())	{
		// this is an elastic object which needs elastic physics simulation
		for(suboit=oit->second.get_obj().begin(); suboit!=oit->second.get_obj().end(); ++suboit)	{
			path = suboit->first;
			std::map<int,GLuint>::iterator tit;
			std::string nrmMap,texMap;
			for(suboit2=suboit->second.get_obj().begin(); suboit2!=suboit->second.get_obj().end(); ++suboit2)	{
                if(suboit2->first=="textureMap")	{
					tex = std::to_string((unsigned long long)suboit2->second.get_int());
					texMap = std::to_string((unsigned long long)suboit2->second.get_int());
				}
				else if(suboit2->first=="normalMap")	{
					nrm = std::to_string((unsigned long long)suboit2->second.get_int());
                    nrmMap = std::to_string((unsigned long long)suboit2->second.get_int());					
				}				
				else {
					throw std::runtime_error("Incorrect fixed triangle section in .smd input file-");
				}
			}
			if(texMap.empty() || nrmMap.empty()) {
				throw std::runtime_error("Missing texture or normal map in fixed triangle section in .smd input file-");
			}
			fixed = false;
			_tuvw = _sg->getTrianglesUVW();
            std::string obj_data = _surgAct->loadData( path, "STATIC" );
			if(_tuvw->readObjFile(obj_data.c_str())) {
				throw std::runtime_error("Unable to load fixed uvwTriangle .obj input file-"); }
			_sg->setWxGraphics(_wxg);
			_sg->setNewTopology();
			_incis.setPreferredEdgeLength(_sg->getMeanEdgeTriangleLength()); // only make this sg call after setting new topology
			//_sg->setTopTextureFilesCreateProgram(texMap.c_str(),nrmMap.c_str());
			path = suboit->first;
			size_t pos = path.rfind(".obj");
			path.erase(pos);
			_sg->setName(path.c_str());
		}
	}
	if((oit=smdObj.find("muscle_layers"))!=smdObj.end())	{
        std::cout << "Loading muscle layers..." << std::endl;
		// this is an muscle glslTriangle which needs elastic physics simulation
        typedef json_spirit::mObject::const_iterator json_iterator;
        typedef json_spirit::mArray json_array;

		for(json_iterator muscle=oit->second.get_obj().begin(); muscle!=oit->second.get_obj().end(); ++muscle)	{
            std::cout << "Loading muscle object: " << muscle->first << std::endl;
			path = muscle->first;
            int muscle_id;
            float muscle_maxstress;
            std::vector<float> muscle_fiber; muscle_fiber.resize(3);
            const json_array& muscle_data = muscle->second.get_array();
            muscle_id = muscle_data[0].get_int();
            for( int v = 0; v < 3; v++ )
                muscle_fiber[v] = (float)muscle_data[1].get_array()[v].get_real();
            muscle_maxstress = (float)muscle_data[2].get_real();

            std::vector< float > temp_vertices;
            std::vector< int > temp_triangles;
            staticTriangle tempTr;
            // NATHAN - just need obj reader. Don't want to visualize it as already loaded that fixed .obj file.
            std::string obj_data = _surgAct->loadData( path, "STATIC" );
            if(tempTr.readObjFile(obj_data.c_str(),true)>0)
                {
                    throw std::runtime_error("Unable to load muscle triangle .obj input file-");
                }
            tempTr.getTriangulatedSurface(temp_vertices, temp_triangles);
            
            _muscleVerts.push_back( temp_vertices );
            _muscleTris.push_back( temp_triangles );
            _muscleFibers.push_back( muscle_fiber );
            _muscleMaxStress.push_back( muscle_maxstress );
                        
        }
    }
/*
	if((oit=smdObj.find("fixedWireframe"))!=smdObj.end())	{
		// COURT - you haven't debugged this JSON version yet.
		json_spirit::mArray::const_iterator ait;
		for(ait=oit->second.get_array().begin(); ait!=oit->second.get_array().end(); ++ait)	{
			path = _smdDirectory + ait->get_str();
			fixed = true;
			staticTriangle trTemp(false);
			if(trTemp.readObjFile(path.c_str(),true))	{
				throw std::runtime_error("Unable to load triangle wireframe .obj input file-");
			}
			std::vector<GLuint> lines;
			trTemp.makeLineList(lines);
			_wxg->getLines()->addLines(trTemp.getPositionsArray(),lines);
			path = ait->get_str();
			size_t pos = path.rfind(".obj");
			path.erase(pos);
			_wxg->getLines()->setName(path.c_str());	
			float color[4]={1.0f,1.0f,1.0f,1.0f};
			_wxg->getLines()->setColor(color);
		}
	}
*/
	if((oit=smdObj.find("incisionSideTexture"))!=smdObj.end())	{
		throw std::runtime_error("There is no incisionSideTexture field in a UVW scene file-");
	}
	if((oit=smdObj.find("physicsNodeSpacing"))!=smdObj.end())	{
		_cleNodeSpacing = (float)oit->second.get_real();
	}
	if((oit=smdObj.find("hookStiffness"))!=smdObj.end())	{
		_hookStiffness = (float)oit->second.get_real();
	}
	if((oit=smdObj.find("sutureStiffness"))!=smdObj.end())	{
		_sutureStiffness = (float)oit->second.get_real();
	}
	if((oit=smdObj.find("poissonsRatio"))!=smdObj.end())	{
		_poissons_ratio = (float)oit->second.get_real();
	}
	if((oit=smdObj.find("youngsModulus"))!=smdObj.end())	{
		_youngs_modulus = (float)oit->second.get_real();
	}
	if((oit=smdObj.find("refinementRatio"))!=smdObj.end())	{
		_refinement_ratio = (float)oit->second.get_real();
	}
	if((oit=smdObj.find("incisionWidth"))!=smdObj.end())	{
		float incisWidth = (float)oit->second.get_real();
        _refinement = std::ceil( _cleNodeSpacing / (incisWidth / _refinement_ratio) );
		_incis.setIncisionWidth(incisWidth);	// incision width also set by this parameter
		_incis.setWxGraphics(_wxg);
	}
	if((oit=smdObj.find("fixedGeometry"))!=smdObj.end())	{
		for(suboit=oit->second.get_obj().begin(); suboit!=oit->second.get_obj().end(); ++suboit)	{
			json_spirit::mArray::const_iterator ait;
			if(suboit->first=="fixedTriangles")	{
				const json_spirit::mArray& fArr = suboit->second.get_array();
				for(ait=fArr.begin(); ait!=fArr.end(); ++ait)
					_fixedTriangles.push_back(ait->get_int());
			}
			if(suboit->first=="lineSegments")	{
				const json_spirit::mArray& fArr = suboit->second.get_array();
				for(ait=fArr.begin(); ait!=fArr.end(); ++ait)
					_fixedSegments.push_back(ait->get_int());
			}
			if(suboit->first=="fixedVertices")	{
				// COURT - haven't debugged this section yet.
				const json_spirit::mArray& fArr = suboit->second.get_array();
				for(ait=fArr.begin(); ait!=fArr.end(); ++ait)
					_fixedPoints.push_back(ait->get_int());
			}
			if(suboit->first=="regions")	{
				const json_spirit::mArray& regionArray = suboit->second.get_array();
                for(ait=regionArray.begin(); ait!=regionArray.end(); ++ait){
                    const json_spirit::mArray& subregion = ait->get_array();
                    REGION region;
                    std::vector<float>& min_corner = region.first;
                    std::vector<float>& max_corner = region.second;
                    for(int i = 0; i < 3; i ++){
                        min_corner.push_back( (float)subregion[0].get_array()[i].get_real() );
                        max_corner.push_back( (float)subregion[1].get_array()[i].get_real() );
                    }
                    _fixedRegions.push_back( region );
                }
			}
			if(suboit->first=="mesh")	{
                path = suboit->second.get_str();
                fixed = true;
                staticTriangle tempTr(false);
                std::string obj_data = _surgAct->loadData( path, "STATIC" );
                if(tempTr.readObjFile(obj_data.c_str(),true)>0)
                    {
                        throw std::runtime_error("Unable to load fixedmesh triangle .obj input file-");
                    }
                tempTr.getTriangulatedSurface(_fixedmeshVerts,_fixedmeshTriangles);
			}
		}
	}
//	std::vector<int> tris(_tuvw->getTriangleArray()->begin(),_tuvw->getTriangleArray()->end());
	//_wxg->frameScene(true);

#ifndef NO_PHYSICS
	setFixedGeometry();	// will be needed later when physics is started.
#endif

    if(_surgAct->getServerMode()>0) { // NATHAN - notice we aren't sending normals or tangents anymore. Computed in javascript client.
        std::vector<int> tris;
        std::vector<float> verts;
        std::vector<float> uvs;
        _sg->getJavascriptData(tris,verts,uvs);
        std::vector<double> d_verts( verts.begin(), verts.end());
        std::vector<double> d_uvs( uvs.begin(), uvs.end());
        _surgAct->getServerHandle().UpdateDynamicData( tris, d_verts, d_uvs);
        _surgAct->getServerHandle().SetTextureName(tex);
        _surgAct->getServerHandle().SetNormalName(nrm);
    }


	return true;
// NATHAN - physics is no longer started here. For the purposes of our demos, it only starts
// when the first hook or suture is added.  In the future physics will be recreated and started
// after incisions and excisions, which are the only surgical elements that change topology.
}

void cleScene::setFixedGeometry()
{ // remember this specification no longer necessarily linked to original model vertices and triangles.
	// next section prepares fixed geometry from input file, which now is related to original model vertices and
	// triangles. This may be dropped at a later date and become totally independent.
	// NATHAN - as discussed this is a messy, probably temporary fix.
	// This should only be done once after initial scene load.
	std::map<int,int> posMap;
	_fixedVerts.clear();
	auto fn = [this,&posMap] (int v) -> int {	std::pair<std::map<int,int>::iterator,bool> pr;
		pr=posMap.insert(std::make_pair(v,(int)(posMap.size())));
		if(pr.second)	{
			float *pos=_tuvw->vertexCoordinate(v);
			for(int j=0; j<3; ++j)
				_fixedVerts.push_back(pos[j]);
		}
		return pr.first->second;
	};
	std::transform(_fixedPoints.begin(),_fixedPoints.end(),_fixedPoints.begin(),fn);
	std::transform(_fixedSegments.begin(),_fixedSegments.end(),_fixedSegments.begin(),fn);
	std::transform(_fixedTriangles.begin(),_fixedTriangles.end(),_fixedTriangles.begin(),fn);
	// Nathan - uncomment next lines for a graphics debug.
/*	std::vector<GLfloat> gPoints;
	gPoints.assign((_fixedVerts.size()/3)<<2,1.0f);
	for(int i=0; i<_fixedVerts.size()/3; ++i) {
		for(int j=0; j<3; ++j)
			gPoints[(i<<2)+j] = _fixedVerts[i*3+j];
	}
	std::vector<GLuint> gLines;
	for(int i=0; i<_fixedSegments.size(); i+=2)	{
		gLines.push_back(_fixedSegments[i]);
		gLines.push_back(_fixedSegments[i+1]);
		gLines.push_back(0xffffffff);
	}
	for(int i=0; i<_fixedTriangles.size(); i+=3)	{
		gLines.push_back(_fixedTriangles[i]);
		gLines.push_back(_fixedTriangles[i+1]);
		gLines.push_back(_fixedTriangles[i+2]);
		gLines.push_back(0xffffffff);
	}
	_wxg->getLines()->addLines(gPoints,gLines); */
}

void cleScene::cleUpdate()
{
    if(_pausePhysics)
        return;

    #if !defined(NO_PHYSICS) && defined(USE_RPC_INTERFACE)
    if( rpc_initialized == false)  return;
    #endif

	// EFTYCHIOS - In non-multithread version frame should be advanced and graphics updated before each redraw.
	// If multithreaded, I am looking for you to set a flag to indicate new CLE vertex positions
	if(_tuvw==NULL)	return;

#ifdef NO_PHYSICS

/*	// oscillates the 0 xcoord y components
  	unsigned int i,j,nVerts = _tr->uniqueCoordinateNumber();
	float vtx[3],ymin=0.0f,ymax=0.0f;
	const float *vptr;
	for(i=0; i<nVerts; ++i)
	{
		vptr = _tr->getVertexCoordinate(i);
		vtx[0]=vptr[0]; vtx[1]=vptr[1]; vtx[2]=vptr[2];
		if(vtx[0]==0.0f)	{
			vtx[1] += _incr;
			ymin = std::min(ymin,vtx[1]);
			ymax = std::max(ymax,vtx[1]);
			_tr->setVertexCoordinate(i,vtx);
		}
	}
	if(ymax>0.4f)
		_incr = -0.001f;
	if(ymin<-0.4f)
		_incr = 0.001f; */

	return;
#else
	if(!_physicsRunning)
		return;

    bool update_available=false;
    int nTverts;
    float *vPos = _tuvw->getPositionArray(nTverts);
    nTverts *= 3;
#if defined(USE_RPC_INTERFACE)
    Vertex_State vstate;
    int frame = _frame;
	client->Advance_One_Time_Step(RPC_Session);
    client->Get_FrameVertices(vstate, RPC_Session, frame);
    if( vstate.vertices.size() == nTverts ){
        update_available = true;
        for(int i=0; i<nTverts; i++)
            vPos[i] = vstate.vertices[i];
        _frame = vstate.frame;
    }
#else
	_cle.Advance_One_Time_Step();
    std::vector<double> d_cleVerts;
    int frame = _frame;
    int fields;
	_cle.Get_Vertex_Data(d_cleVerts, V_DATA_TYPE::POSITION, frame);
    if( d_cleVerts.size() > 0 ){
        update_available = true;
        fields = d_cleVerts.size() / nTverts;
        std::cout << "Updated " << fields << " vertex data fields." << std::endl;
        for(int i=0; i<nTverts; i++)
            vPos[i] = d_cleVerts[i];
        _frame = frame;
    }
#endif
    if( update_available ){
        //std::cout << "Current Physics Frame : " << _frame << std::endl;
/*        int i,numVerts;
          float *posArr = _tuvw->getPositionArray(numVerts);
          for(i=0; i<numVerts; ++i) {
          posArr[i<<2] = _cleVerts[(i<<1)+i];
          posArr[(i<<2)+1] = _cleVerts[(i<<1)+i+1];
          posArr[(i<<2)+2] = _cleVerts[(i<<1)+i+2];
          } */
        _sg->updatePositionsAndNormals();
        if(_surgAct->getServerMode()>0) {
            std::vector<float> scratch_remapped;
            std::vector<float> scratch_raw;
            std::vector<double> d_stress;
            std::vector<double> d_strain;
            std::vector<double> d_position;
            
            std::vector<float> data_texture;
            
            for(int i=0; i<nTverts; i++)
                scratch_raw.push_back(d_cleVerts[i]); // Position is the first third of the vertex data
            _sg->remapJavascriptVertexData(scratch_raw,scratch_remapped);
            d_position.assign(scratch_remapped.begin(),scratch_remapped.end());            
            
            if(fields>=2){
                scratch_raw.clear();
                scratch_remapped.clear();
                
                for(int i=nTverts; i<nTverts*2; i++){
                    scratch_raw.push_back(d_cleVerts[i]); // Stress is the second third of the vertex data
                }
                _sg->remapJavascriptVertexData(scratch_raw,scratch_remapped);          
                d_stress.assign(scratch_remapped.begin(),scratch_remapped.end());
            }
            
            if( fields >= 3){
                scratch_raw.clear();
                scratch_remapped.clear();
                
                for(int i=nTverts*2; i<nTverts*3; i++){
                    scratch_raw.push_back(d_cleVerts[i]); // Strain is the last third of the vertex data
                }
                _sg->remapJavascriptVertexData(scratch_raw,scratch_remapped);          
                d_strain.assign(scratch_remapped.begin(),scratch_remapped.end());
            }
            
            _cle.GetTextureStress(data_texture);            
            _surgAct->getServerHandle().UpdateVertexData( d_position, d_stress, d_strain, data_texture);
        }
    }        
#endif
}

int cleScene::addHook(int triangle, float (&uv)[2], float(&pos)[3] )
{
	if(!_physicsRunning)
		createAndStartNewPhysicsModel();
	int ret=0;
	float vCoord[3];
	_tuvw->getBarycentricPosition(triangle,uv,vCoord);
    if(pos){
        pos[0] = vCoord[0]; pos[1] = vCoord[1]; pos[2] = vCoord[2]; }
#ifndef NO_PHYSICS
#if defined(USE_RPC_INTERFACE)
    std::vector<double> dweights(uv, uv+2);
    std::vector<double> dloc(vCoord, vCoord+3);

	ret = client->Add_Hook(RPC_Session,triangle,dweights);
	client->Move_Hook(RPC_Session,ret,dloc);
#else
	ret = _cle.Add_Hook(triangle,uv);
	_cle.Move_Hook(ret,(float(&)[3])vCoord);
#endif
#endif
	return ret;
}

void cleScene::deleteHook(const int hook_id){
#if !defined(NO_PHYSICS) && !defined(USE_RPC_INTERFACE)
    _cle.Delete_Hook(hook_id); _cle.Advance_One_Time_Step();
#elif !defined(NO_PHYSICS) && defined(USE_RPC_INTERFACE)
    client->Delete_Hook(RPC_Session,hook_id); client->Advance_One_Time_Step(RPC_Session);
#endif
	// NATHAN - what follows is real bad. We have allowed physics nodes to be different than user interface nodes.
	// These must be made unique and consistent between the two environments. Problem began back with Joey Teran and initial interface design.
}

void cleScene::setHookPosition(const float(&location)[3],const int hook_id){
#if !defined(NO_PHYSICS) && !defined(USE_RPC_INTERFACE)
    _cle.Move_Hook(hook_id,location); 
#elif !defined(NO_PHYSICS) && defined(USE_RPC_INTERFACE)
        double dlocation[3];
        dlocation[0] = location[0];
        dlocation[1] = location[1];
        dlocation[2] = location[2];
        const std::vector<double> vlocation(dlocation, dlocation+3);
        client->Move_Hook(RPC_Session,hook_id,vlocation);
#endif
	// NATHAN - what follows is real bad. We have allowed physics nodes to be different than user interface nodes.
	// These must be made unique and consistent between the two environments. Problem began back with Joey Teran and initial interface design.
}


int cleScene::addSuture(int triangle0, float (&uv0)[2], int triangle1, float (&uv1)[2])
{	// can only add on a triangle edge
	if(!_physicsRunning)
		createAndStartNewPhysicsModel();

	int ret=-1;
/*
	float weights0[2],weights1[2];
	if(edge0<1)	{
		weights0[0]=param0;
		weights0[1]=0.0f;	}
	else if(edge0<2)	{
		weights0[0]=1.0f-param0;
		weights0[1]=param0;	}
	else	{
		weights0[0]=0.0f;
		weights0[1]=1.0f-param0; 	}
	if(edge1<1)	{
		weights1[0]=param1;
		weights1[1]=0.0f;	}
	else if(edge1<2)	{
		weights1[0]=1.0f-param1;
		weights1[1]=param1;	}
	else	{
		weights1[0]=0.0f;
		weights1[1]=1.0f-param1; 	}
*/
#ifndef NO_PHYSICS
#if defined(USE_RPC_INTERFACE)
    std::vector<double> dweights0(uv0, uv0+2);
    std::vector<double> dweights1(uv1, uv1+2);
	ret = client->Add_Suture(RPC_Session,triangle0,dweights0,triangle1,dweights1);
#else
	ret = _cle.Add_Suture(triangle0,uv0,triangle1,uv1);
#endif
#endif
	return ret;
}

void cleScene::deleteSuture(const int suture_id){
#if !defined(NO_PHYSICS) && !defined(USE_RPC_INTERFACE)
    _cle.Delete_Suture(suture_id); _cle.Advance_One_Time_Step();
#elif !defined(NO_PHYSICS) && defined(USE_RPC_INTERFACE)
    client->Delete_Suture(RPC_Session,suture_id); client->Advance_One_Time_Step(RPC_Session);
#endif
}

bool cleScene::makeIncision(trianglesUVW* tri,int nPoints,float (*positions)[3],float (*normals)[3],
                            int TstartTriangleEdge,float TstartParam,bool Tout)
{
	if(!_incis.makeIncision(tri,nPoints,positions,normals,TstartTriangleEdge,TstartParam,Tout))
		return false;
	_sg->setNewTopology();
    if(_surgAct->getServerMode()>0) {
        std::vector<int> tris;
        std::vector<float> verts;
        std::vector<float> uvs;
        _sg->getJavascriptData(tris,verts,uvs);
        std::vector<double> d_verts( verts.begin(), verts.end());
        std::vector<double> d_uvs( uvs.begin(), uvs.end());
        _surgAct->getServerHandle().UpdateDynamicData( tris, d_verts, d_uvs );
        // Use these to clear out manipulators
        _surgAct->getServerHandle().UpdateHooks( std::vector< std::pair< int, std::vector<double> > >() );
        _surgAct->getServerHandle().UpdateSutures( std::vector< std::pair< std::vector<int>, std::vector<double> > > () );
    }
// NATHAN - once we have a robust incision tool, we'll turn the next line back on.  For our initial
// demos the physics won't get turned on until after all incisions and excisions are done and a hook
// or suture has been applied.
//	createAndStartNewPhysicsModel();
#if !defined(NO_PHYSICS) && defined(USE_RPC_INTERFACE)
	if(_physicsRunning)
		client->Destroy_Model(RPC_Session);	// no need to destroy a non-existent model
#elif !defined(NO_PHYSICS) && !defined(USE_RPC_INTERFACE)
	if(_physicsRunning)	// no need to destroy a non-existent model
		_cle.Destroy_Model();
#endif
	_physicsRunning=false;

	return true;
}

bool cleScene::exciseSubobject(int triangle)
{	// removes entire subobject connected to triangle
	_tuvw->findAdjacentTriangles();
	recurseTriangleRemoval(triangle);
	_tuvw->cleanAndPack();
	_tuvw->findAdjacentTriangles();
	_sg->setNewTopology();
	_sg->updatePositionsAndNormals();
    if(_surgAct->getServerMode()>0) {
        std::vector<int> tris;
        std::vector<float> verts;
        std::vector<float> uvs;
        _sg->getJavascriptData(tris,verts,uvs);
        std::vector<double> d_verts( verts.begin(), verts.end());
        std::vector<double> d_uvs( uvs.begin(), uvs.end());
        _surgAct->getServerHandle().UpdateDynamicData( tris, d_verts, d_uvs );
    }
// NATHAN - once we have a robust incision tool, we'll turn the next line back on.  For our initial
// demos the physics won't get turned on until after all incisions and excisions are done and a hook
// or suture has been applied.
//	createAndStartNewPhysicsModel();
	return true;
}

void cleScene::recurseTriangleRemoval(int triangle)
{
	int *tp = _tuvw->triangleVertices(triangle);
	if(*tp < 0)
		return;
	unsigned long *adjs = _tuvw->triAdjs(triangle);
	*tp = -1;
	for(int i=0; i<3; ++i)
		recurseTriangleRemoval(adjs[i]>>2);
}

void cleScene::createAndStartNewPhysicsModel()
{	// destroys previous physics model and creates new one
    if(_surgAct->getServerMode()>0){
        _surgAct->getServerHandle().SetSimulatorState(1);
    }

#ifndef NO_PHYSICS
    std::vector<int> js_tris;
    std::vector<float> js_verts;
    std::vector<float> js_uvs;
    std::vector<int> js_map;
    _sg->getJavascriptData(js_tris,js_verts,js_uvs);
    _sg->getJavascriptMap(js_map);
	std::vector<int> *tris = _tuvw->getTriangleArray();
	std::vector<float> *verts = _tuvw->getPositionArray();
    std::vector<float> *texes = _tuvw->getTextureArray();
    std::vector<int> texture_size;
    texture_size.push_back(256);
    texture_size.push_back(256);
//	_tuvw->getTriangleArray(_cleVerts,tris);
    if(_surgAct->getServerMode()>0) {
        std::vector<double> d_verts( js_verts.begin(), js_verts.end());
        std::vector<double> d_uvs( js_uvs.begin(), js_uvs.end());
        _surgAct->getServerHandle().UpdateDynamicData( js_tris, d_verts, d_uvs );
    }

#if 1
#if defined(USE_RPC_INTERFACE)
	if(_physicsRunning)
		client->Destroy_Model(RPC_Session);	// no need to destroy a non-existent model
    std::vector<double> d_cleVerts(verts->begin(), verts->end());
	client->Set_Hook_Stiffness(RPC_Session,_hookStiffness);
	client->Set_Suture_Stiffness(RPC_Session,_sutureStiffness);
    client->Set_Poissons_Ratio(_poissons_ratio);
    client->Set_Youngs_Modulus(_youngs_modulus);
	client->Create_Model(RPC_Session,d_cleVerts,*tris,_cleNodeSpacing,_refinement);    
#else
	if(_physicsRunning)	// no need to destroy a non-existent model
		_cle.Destroy_Model();
    _cle.Set_Poissons_Ratio(_poissons_ratio);
    _cle.Set_Youngs_Modulus(_youngs_modulus);
	_cle.Set_Hook_Stiffness(_hookStiffness);
	_cle.Set_Suture_Stiffness(_sutureStiffness);
	_cle.Create_Model(*verts,*tris,_cleNodeSpacing,_refinement);
    _cle.SetTextureParameters(texture_size, js_uvs, js_tris, js_map);
#endif
#else
#if defined(USE_RPC_INTERFACE)
	if(_physicsRunning)
		client->Destroy_Model(RPC_Session);	// no need to destroy a non-existent model
    std::vector<double> d_cleVerts(verts->begin(), verts->end());
	client->Set_Hook_Stiffness(RPC_Session,_hookStiffness);
	client->Set_Suture_Stiffness(RPC_Session,_sutureStiffness);
    client->Set_Poissons_Ratio(_poissons_ratio);
    client->Set_Youngs_Modulus(_youngs_modulus);
	client->Create_Model(RPC_Session,d_cleVerts,tris,_cleNodeSpacing,_refinement);
#else
	if(_physicsRunning)	// no need to destroy a non-existent model
		_cle.Destroy_Model();
    _cle.Set_Poissons_Ratio(_poissons_ratio);
    _cle.Set_Youngs_Modulus(_youngs_modulus);
	_cle.Set_Hook_Stiffness(_hookStiffness);
	_cle.Set_Suture_Stiffness(_sutureStiffness);
	_cle.Create_Model(*verts,*tris,_cleNodeSpacing,_refinement);
    _cle.SetTextureParameters(texture_size, js_uvs, js_tris, js_map);
#endif
#endif

	// COURT - below is correct. The Static models are just that: static. They allow us to display
        // the scenery in the outside physbam visualizer. They do not effect the physics simulation.
        std::list<staticTriangle*>::iterator tit;
        for(tit=_triList.begin(); tit!=_triList.end(); ++tit)	{
            std::vector<int> stris;
            std::vector<float> sverts;
            (*tit)->getTriangulatedSurface(sverts,stris);           
#if defined(USE_RPC_INTERFACE)
            std::vector<double> dverts(sverts.begin(), sverts.end());
            client->Add_Static_Model(RPC_Session,dverts, stris);
#else
            _cle.Add_Static_Model(sverts,stris);
#endif
            //_cle.WriteDebug();
        }


	// INDEPENDENT fixed geometry input with scene load. Not necessarily related to an existing
	// vertex, line or triangle in the scene. Already computed, just output it.
#if defined(USE_RPC_INTERFACE)
        {
            std::vector<double> d_verts(_fixedVerts.begin(), _fixedVerts.end());
            client->Set_Fixed_Geometry(RPC_Session,d_verts,_fixedPoints,_fixedSegments,_fixedTriangles);  
        }

    for( int i = 0; i < _fixedRegions.size(); i++ ){
        std::vector<double> d_mincorner( _fixedRegions[i].first.begin(), _fixedRegions[i].first.end() );
        std::vector<double> d_maxcorner( _fixedRegions[i].second.begin(), _fixedRegions[i].second.end() );
        client->Set_Fixed_Volume(RPC_Session, d_mincorner, d_maxcorner );
    }
    
    if( _fixedmeshVerts.size() > 0 ){
        std::vector<double> d_verts(_fixedmeshVerts.begin(), _fixedmeshVerts.end());
        client->Set_Fixed_Triangles(RPC_Session,d_verts,_fixedmeshTriangles);  
    }
#else
    _cle.Set_Fixed_Geometry(_fixedVerts,_fixedPoints,_fixedSegments,_fixedTriangles);

    for( int i = 0; i < _fixedRegions.size(); i++ ){
        _cle.Set_Fixed_Volume( _fixedRegions[i].first, _fixedRegions[i].second );
    }

    if( _fixedmeshVerts.size() > 0 ){
        _cle.Set_Fixed_Triangles(_fixedmeshVerts,_fixedmeshTriangles);  
    }

#endif

	// This is for collisions
	if(!_collisionTris.empty()){
		#if defined(USE_RPC_INTERFACE) //TODO: Need to implement this yet!
		std::vector<double> d_verts(_collisionVerts.begin(),_collisionVerts.end());
			client->Set_Collision_Model(RPC_Session,d_verts,_collisionTris);
		#else
			_cle.Set_Collision_Model(_collisionVerts,_collisionTris);
		#endif
	}

    for( int i = 0; i<_muscleVerts.size(); i++)
    {
		#if defined(USE_RPC_INTERFACE) //TODO: Need to implement this yet!
        std::vector<double> d_muscleVerts(_muscleVerts[i].begin(),_muscleVerts[i].end());
		std::vector<double> d_muscleFibers(_muscleFibers[i].begin(),_muscleFibers[i].end());
        client->Add_Muscle_Layer(RPC_Session, d_muscleVerts, _muscleTris[i], d_muscleFibers, _muscleMaxStress[i]  );
        #else
        _cle.Add_Muscle_Layer( _muscleVerts[i], _muscleTris[i], _muscleFibers[i], _muscleMaxStress[i]  );
        #endif
    }

#if defined(USE_RPC_INTERFACE)
    client->Finalize_Initialization(RPC_Session);
#else
	_cle.Finalize_Initialization();
#endif
    //_cle.WriteDebug();
	_physicsRunning=true;
    _frame = -1;
	// NATHAN - I should also write repeat calls from surgical history to reapply hooks and sutures - not done yet.
#endif	// NO_PHYSICS

    if(_surgAct->getServerMode()>0)
        _surgAct->getServerHandle().SetSimulatorState(2);
}

void cleScene::getBoundingBox(float (&box)[6])
{
	_bBox.Empty_Box();
	std::list<staticTriangle*>::iterator git;
	for(git=_triList.begin(); git!=_triList.end(); ++git)	{
		int i,n=(*git)->vertexNumber();
		float xyz[3];
		for(i=0; i<n; ++i)	{
			(*git)->getVertexCoordinate(i,xyz);
			_bBox.Enlarge_To_Include_Point(xyz);
		}
	}
	box[1]=_bBox.xmax; box[0]=_bBox.xmin; 
	box[3]=_bBox.ymax; box[2]=_bBox.ymin; 
	box[5]=_bBox.zmax; box[4]=_bBox.zmin; 
}

void cleScene::scaleScene(float scaleFactor)
{
	assert(false);  // do we need this anymore?
	std::list<staticTriangle*>::iterator git;
	for(git=_triList.begin(); git!=_triList.end(); ++git)	{
		int i,n=(*git)->vertexNumber();
		GLfloat *xyz = (*git)->getPositionArray(n);
		for(i=0; i<n; ++i)	{
			for(int j=0; j<3; ++j)
				xyz[(i<<2)+j] *= scaleFactor;
		}
	}
}

void cleScene::saveModifiedScene()
{
	assert(false); // can we nuke this now?
	std::list<staticTriangle*>::iterator git;
	for(git=_triList.begin(); git!=_triList.end(); ++git)	{
		std::string fname,path(_smdDirectory);
		(*git)->getName(fname);
		size_t pos = fname.rfind(".obj");
		if(pos<1000)
			fname.erase(pos,std::string::npos);
		path.append(fname);
		path.append("_mod.obj");
//		(*git)->writeObjFile(path.c_str());
	}
}

void cleScene::createFixedNodeFile()
{
	std::string path(_smdDirectory);
	if(path.empty())
		return;
	path.append("fixedNodes.txt");
	std::ofstream fin(path.c_str());
    if(!fin.is_open())
        return;
	char s[400];
	// next section finds Aaron's "corner" polylines
	float bbox[6];
	getBoundingBox(bbox);
	std::vector<int> extremeV[4];
	int nv;
	const float *pos = _tuvw->getPositionArray(nv);
	 nv *= 3;
	for(int i=0; i<nv; i+=3)	{
		if(pos[i]<=bbox[0] && pos[i+2]<=bbox[4])
			extremeV[0].push_back(i/3);
		if(pos[i]<=bbox[0] && pos[i+2]>=bbox[5])
			extremeV[1].push_back(i/3);
		if(pos[i]>=bbox[1] && pos[i+2]<=bbox[4])
			extremeV[2].push_back(i/3);
		if(pos[i]>=bbox[1] && pos[i+2]>=bbox[5])
			extremeV[3].push_back(i/3);
	}
/*	for(int i=0; i<4; ++i)	{  // not necessary any longer as no difference between vertex num and its position array index.
		std::transform(extremeV[i].begin(),extremeV[i].end(),extremeV[i].begin(),[this](int p) -> int{
			int j;
			for(j=0; j<_tr->vertexNumber(); ++j)	{
				if(_tr->getPositionArrayIndex(j)>>2==p)
					break;
			}
			return j;}
		);
	} */
	for(int i=0; i<4; ++i)	{
		std::sort(extremeV[i].begin(),extremeV[i].end(),[pos](int m,int n){return pos[(m*3)+1]<pos[(n*3)+1];});
		sprintf(s,"d l %d\n",(int)(extremeV[i].size()));
		fin.write(s,strlen(s));
		std::for_each(extremeV[i].begin(),extremeV[i].end(),[&fin,&s](int n){sprintf(s,"%d\n",n); fin.write(s,strlen(s));});
	}
	fin.close();
}

