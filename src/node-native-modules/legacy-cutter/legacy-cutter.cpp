#include <set>
#include <iostream>
#include <vector>
#include <array>

#include "legacy-cutter.h"
#include <GraphicsUtils/skinGraphics.h>
#include <GraphicsUtils/staticTriangle.h>

#include "incision.h"

namespace {
template<class T>
std::vector<T> load_std_vector( const Nan::FunctionCallbackInfo<v8::Value>& info, int slot ){
    std::vector<T> T_vec;
    Nan::HandleScope scope;
    v8::Handle<v8::Value> val;
    
    if (info[slot]->IsArray()) {
        v8::Handle<v8::Array> jsArray = v8::Handle<v8::Array>::Cast(info[slot]);
        for (unsigned int i = 0; i < jsArray->Length(); i++) {
            val = jsArray->Get(i);
            T_vec.push_back(val->NumberValue());
        }
    }
    else{       
        v8::Local<v8::Object> BufObj    = info[slot]->ToObject();
        char*         BufData   = node::Buffer::Data(BufObj);
        size_t        BufLength = node::Buffer::Length(BufObj);
        T* T_Data = (T*)(BufData);
        T_vec.reserve( (BufLength / sizeof(T)) );
        for( int i = 0; i < (BufLength / sizeof(T)); i++)
            T_vec.push_back( T_Data[i] );
    }
    return T_vec;
}

template<class T, int size>
std::array<T,size> load_std_array( const Nan::FunctionCallbackInfo<v8::Value>& info, int slot ){
    std::array<T,size> T_arr; 
    Nan::HandleScope scope;
    v8::Handle<v8::Value> val;
    
    if (info[slot]->IsArray()) {
        v8::Handle<v8::Array> jsArray = v8::Handle<v8::Array>::Cast(info[slot]);
        if( jsArray->Length() != size )
            Nan::ThrowRangeError("Array must be of fixed length.");
        for (unsigned int i = 0; i < jsArray->Length(); i++) {
            val = jsArray->Get(i);
            T_arr[i] = val->NumberValue();
        }
    }
    else {
        v8::Local<v8::Object> BufObj    = info[slot]->ToObject();
        char*         BufData   = node::Buffer::Data(BufObj);
        size_t        BufLength = node::Buffer::Length(BufObj);
        T* T_Data = (T*)(BufData);
        if( (BufLength / sizeof(T)) != size )
            Nan::ThrowRangeError("Array must be of fixed length.");
        for( int i = 0; i < (BufLength / sizeof(T)); i++)
            T_arr[i] = T_Data[i];
    }
    return T_arr;
}

template<class T>
void save_std_vector( const Nan::FunctionCallbackInfo<v8::Value>& info, int slot, const std::vector<T>& T_vec ){
    v8::Local<v8::Object> BufObj    = info[slot]->ToObject();
    char*         BufData   = node::Buffer::Data(BufObj);
    size_t        BufLength = node::Buffer::Length(BufObj);
    T* T_Data = (T*)(BufData);
    if( T_vec.size() > (BufLength/sizeof(T)) )
        Nan::ThrowRangeError("Insufficient space allocated for array copy.");
    for( int i = 0; i < T_vec.size(); i++)
        T_Data[i] = T_vec[i];   
}
}

void InitAll(v8::Local<v8::Object> exports) {
  Legacy_Cutter::Init(exports);
}

NODE_MODULE(Legacy_Cutter, InitAll)

Nan::Persistent<v8::Function> Legacy_Cutter::constructor;

Legacy_Cutter::Legacy_Cutter()  {
    std::cout << "Legacy_Cutter::Constructor" << std::endl;
    _incis = new incision();
    _sg = new skinGraphics();
}

Legacy_Cutter::~Legacy_Cutter() {
    std::cout << "Legacy_Cutter::Destructor" << std::endl;
    delete _incis;
    delete _sg;
} 


void Legacy_Cutter::Init(v8::Local<v8::Object> exports) {
  Nan::HandleScope scope;

  // Prepare constructor template
  v8::Local<v8::FunctionTemplate> tpl = Nan::New<v8::FunctionTemplate>(New);
  tpl->SetClassName(Nan::New("Legacy_Cutter").ToLocalChecked());
  tpl->InstanceTemplate()->SetInternalFieldCount(1);

  // Prototype
  Nan::SetPrototypeMethod(tpl,  "ParseFile",  __ParseFile);
  Nan::SetPrototypeMethod(tpl,  "ParseStaticFile",  __ParseFileStatic);
  Nan::SetPrototypeMethod(tpl,  "Incise",  __makeIncision);
  Nan::SetPrototypeMethod(tpl,  "Excise",  __makeExcision);
  Nan::SetPrototypeMethod(tpl,  "GetJS_Vertex", __GetJavascriptVertex);
  Nan::SetPrototypeMethod(tpl,  "GetJS_Topology", __GetJavascriptTopology);
  Nan::SetPrototypeMethod(tpl,  "GetJS_UV", __GetJavascriptUV);
  Nan::SetPrototypeMethod(tpl,  "GetRaw_Vertex",__GetRawVertex);
  Nan::SetPrototypeMethod(tpl,  "GetRaw_Topology",__GetRawTopology);

  constructor.Reset(tpl->GetFunction());
  exports->Set(Nan::New("Legacy_Cutter").ToLocalChecked(), tpl->GetFunction());
}

void Legacy_Cutter::New(const Nan::FunctionCallbackInfo<v8::Value>& info) {
  if (info.IsConstructCall()) {
    // Invoked as constructor: `new Legacy_Cutter(...)`
    //double value = info[0]->IsUndefined() ? 0 : info[0]->NumberValue();
    Legacy_Cutter* obj = new Legacy_Cutter();
    obj->Wrap(info.This());
    info.GetReturnValue().Set(info.This());
  } else {
    // Invoked as plain function `Legacy_Cutter(...)`, turn into construct call.
    const int argc = 1;
    v8::Local<v8::Value> argv[argc] = { info[0] };
    v8::Local<v8::Function> cons = Nan::New<v8::Function>(constructor);
    info.GetReturnValue().Set(cons->NewInstance(argc, argv));
  }
}


void Legacy_Cutter::__ParseFile( const Nan::FunctionCallbackInfo<v8::Value>& info ){
    Legacy_Cutter* obj = ObjectWrap::Unwrap<Legacy_Cutter>(info.Holder());   
    std::string tempString(*v8::String::Utf8Value(info[0]));
    trianglesUVW *_tuvw = obj->_sg->getTrianglesUVW();
    if(_tuvw->readObjFile(tempString.c_str()))
        Nan::ThrowError("Unable to load fixed uvwTriangle .obj input file-");
    obj->_sg->setNewTopology();
    obj->_incis->setIncisionWidth(0.0036);
    obj->_incis->setPreferredEdgeLength(obj->_sg->getMeanEdgeTriangleLength());    
}

void Legacy_Cutter::__ParseFileStatic( const Nan::FunctionCallbackInfo<v8::Value>& info ){
    Nan:: HandleScope scope;
    v8::Local<v8::Object> results = Nan::New<v8::Object>();

    Legacy_Cutter* obj = ObjectWrap::Unwrap<Legacy_Cutter>(info.Holder());   
    std::string tempString(*v8::String::Utf8Value(info[0]));
    staticTriangle *tr = new staticTriangle();
    if(tr->readObjFile(tempString.c_str()))
        Nan::ThrowError("Unable to load fixed uvwTriangle .obj input file-");

    std::vector<int> tris;
    std::vector<float> verts;
    std::vector<float> uvs;
    tr->getSurfaceTriangles(tris);
    tr->getSurfaceVertices(verts);
    tr->getSurfaceUVs(uvs);

    delete tr;

    results->Set( Nan::New("vertices").ToLocalChecked(),
                  Nan::CopyBuffer(reinterpret_cast<const char*>(verts.data()),
                                  verts.size()*sizeof(float)).ToLocalChecked() );

    results->Set( Nan::New("topology").ToLocalChecked(),
                  Nan::CopyBuffer(reinterpret_cast<const char*>(tris.data()),
                                  tris.size()*sizeof(int)).ToLocalChecked() );

    results->Set( Nan::New("uv").ToLocalChecked(),
                  Nan::CopyBuffer(reinterpret_cast<const char*>(uvs.data()),
                                  uvs.size()*sizeof(float)).ToLocalChecked() );

    info.GetReturnValue().Set( results );
}


void Legacy_Cutter::__makeIncision( const Nan::FunctionCallbackInfo<v8::Value>& info ){
    Legacy_Cutter* obj = ObjectWrap::Unwrap<Legacy_Cutter>(info.Holder());   

    typedef float T;

    std::vector<int> path_triangles = load_std_vector<int>(info,0);;
    std::vector<float> path_uvs = load_std_vector<T>(info,1);
    std::vector<float> path_positions = load_std_vector<T>(info,2);
    std::vector<float> path_normals = load_std_vector<T>(info,3); 
    bool edge_start = info[4]->BooleanValue();
    bool edge_end = info[5]->BooleanValue();


    std::cout << "Path: " << std::endl;
    for( int i = 0 ; i < path_triangles.size(); i++ ){
        std::cout << "Node " << i;
        std::cout << " Triangle: " << path_triangles.at(i);
        std::cout << " UV: " << path_uvs.at(i*2) << ", " << path_uvs.at(i*2+1);
        std::cout << " Position: " << path_positions.at(i*3) << ", " << path_positions.at(i*3+1) << ", " << path_positions.at(i*3+2);
        std::cout << " Position: " << path_normals.at(i*3) << ", " << path_normals.at(i*3+1) << ", " << path_normals.at(i*3+2);
        std::cout << std::endl;
    }
    std::cout << "T_In: " << edge_start << std::endl;
    std::cout << "T_Out: " << edge_end << std::endl;

    obj->CreateIncision( path_triangles, path_uvs, path_positions, path_normals, edge_start, edge_end);
}

void Legacy_Cutter::__makeExcision( const Nan::FunctionCallbackInfo<v8::Value>& info ){
    Legacy_Cutter* obj = ObjectWrap::Unwrap<Legacy_Cutter>(info.Holder());   

    typedef float T;

    int triangle = info[0]->NumberValue();

    obj->CreateExcision( triangle );
}

void Legacy_Cutter::__GetJavascriptVertex( const Nan::FunctionCallbackInfo<v8::Value>& info )
{
    Legacy_Cutter* obj = ObjectWrap::Unwrap<Legacy_Cutter>(info.Holder());   

    std::vector<int> tris;
    std::vector<float> verts;
    std::vector<float> uvs;
    obj->_sg->getJavascriptData(tris,verts,uvs);
    info.GetReturnValue().Set(Nan::CopyBuffer(reinterpret_cast<const char*>(verts.data()),
                                              verts.size()*sizeof(float)).ToLocalChecked());
}

void Legacy_Cutter::__GetJavascriptTopology( const Nan::FunctionCallbackInfo<v8::Value>& info )
{
    Legacy_Cutter* obj = ObjectWrap::Unwrap<Legacy_Cutter>(info.Holder());   

    std::vector<int> tris;
    std::vector<float> verts;
    std::vector<float> uvs;
    obj->_sg->getJavascriptData(tris,verts,uvs);
    //for(int i=0; i< tris.size(); i++)
    //    tris[i]--;
    info.GetReturnValue().Set(Nan::CopyBuffer(reinterpret_cast<const char*>(tris.data()),
                                              tris.size()*sizeof(int)).ToLocalChecked());
}

void Legacy_Cutter::__GetJavascriptUV( const Nan::FunctionCallbackInfo<v8::Value>& info )
{
    Legacy_Cutter* obj = ObjectWrap::Unwrap<Legacy_Cutter>(info.Holder());   

    std::vector<int> tris;
    std::vector<float> verts;
    std::vector<float> uvs;
    obj->_sg->getJavascriptData(tris,verts,uvs);
    info.GetReturnValue().Set(Nan::CopyBuffer(reinterpret_cast<const char*>(uvs.data()),
                                              uvs.size()*sizeof(float)).ToLocalChecked());
}


void Legacy_Cutter::__GetRawVertex( const Nan::FunctionCallbackInfo<v8::Value>& info )
{
    Legacy_Cutter* obj = ObjectWrap::Unwrap<Legacy_Cutter>(info.Holder());   
    trianglesUVW *tri = obj->_sg->getTrianglesUVW();
    std::vector<float> *verts = tri->getPositionArray();
    info.GetReturnValue().Set(Nan::CopyBuffer(reinterpret_cast<const char*>(verts->data()),
                                              verts->size()*sizeof(float)).ToLocalChecked());
}

void Legacy_Cutter::__GetRawTopology( const Nan::FunctionCallbackInfo<v8::Value>& info )
{
    Legacy_Cutter* obj = ObjectWrap::Unwrap<Legacy_Cutter>(info.Holder());   
    trianglesUVW *tri = obj->_sg->getTrianglesUVW();
    std::vector<int> *tris = tri->getTriangleArray();
    info.GetReturnValue().Set(Nan::CopyBuffer(reinterpret_cast<const char*>(tris->data()),
                                              tris->size()*sizeof(int)).ToLocalChecked());
}


bool Legacy_Cutter::texturePickCode(const int triangle, const float (&uv)[2], float (&uvw)[3], float &triangleDuv){
    // This is required by Court's old incision glue code. We may not need this as we refactor his stuff away.   

    // texture seam triangles will have a large deltaUV in cylindrical or spherical texture mapping
	// return true=user selected a top or bottom triangle, false if an edge triangle was selected
	float *tx[3],mm[6]={1e30f,-1e30f,1e30f,-1e30f,1e30f,-1e30f},p=1.0f-uv[0]-uv[1];
	trianglesUVW *tri = _sg->getTrianglesUVW();
	int *tr = tri->triVerts(triangle);
	if(*tr<0)
		return false;
	for(int j,i=0; i<3; ++i) {
		tx[i] = tri->vertexTexture(tr[i]);
		for(j=0; j<3; ++j) {
			if(mm[j<<1]>tx[i][j])
				mm[j<<1]=tx[i][j];
			if(mm[(j<<1)+1]<tx[i][j])
				mm[(j<<1)+1]=tx[i][j];
		}
	}
	for(int i=0; i<3; ++i)
		uvw[i] = p*tx[0][i] + uv[0]*tx[1][i] + uv[1]*tx[2][i];
	triangleDuv = mm[1]-mm[0] + mm[3]-mm[2];
	if(mm[5]-mm[4]>1e-5f)
		return false;
	if(fabs(1.0f-mm[4])<1e-5f || mm[5]<1e-5f)
		return true;
	return false;
}


void Legacy_Cutter::CreateIncision( const std::vector<int>& path_triangles, const std::vector<float>& path_uvs, 
                                    const std::vector<float>& path_positions, const std::vector<float>& path_normals, 
                                    const bool edge_start, const bool edge_end ){
    // This is Courts old incision tool glue code. It will probably need to be rewritten as we transition to a new
    // incision / cutting code.
    
    bool Tout=false,nukeThis=false;
    int n = path_triangles.size();
    trianglesUVW *tri= _sg->getTrianglesUVW();
    std::vector<float> hitV,hitP;
    std::vector<int> hitTr;
    int TstartTriangleEdge=0x0003;
    float TstartParam=-1.0f;
    float uv[2],uvw[3],triangleDuv,minParam=1.0e15f;
    for(int j,i=0; i<n; ++i)	{ // nH,
        uv[0] = path_uvs[i<<1];
        uv[1] = path_uvs[(i<<1)+1];
        if(!texturePickCode(path_triangles.at(i),uv,uvw,triangleDuv)) // user picked an edge triangle
            Nan::ThrowError("Incision tool error: User picked an edge triangle");         

        if((i==0 && edge_start) )	{	// ||(i==n-1 && edge_end) removed because in T-out self intersection there is no getOppositeIncisionEdge() yet
            if(path_uvs[(i<<1)+1]==0.0f)	{
                j = 0;
                uv[0] = path_uvs[i<<1];
            }
            else if(path_uvs[i<<1]==0.0f)	{
                j = 2;
                uv[0] = 1.0f - path_uvs[(i<<1)+1];
            }
            else	{
                j = 1;
                uv[0] = path_uvs[(i<<1)+1];
            }
            TstartTriangleEdge = path_triangles[i]<<2;
            TstartTriangleEdge |= j;
            TstartParam = uv[0];
        }
    }

    if(!_incis->makeIncision(tri,
                            n,
                            (float (*)[3])&(path_positions[0]),
                            (float (*)[3])&(path_normals[0]),
                            TstartTriangleEdge,
                            TstartParam,
                            edge_end))
        Nan::ThrowError("Incision tool error. Not recoverable");
    
    _sg->setNewTopology();        
}

void Legacy_Cutter::CreateExcision( int triangle ){
    trianglesUVW *_tuvw = _sg->getTrianglesUVW();
    
	_tuvw->findAdjacentTriangles();
	recurseTriangleRemoval(triangle);
	_tuvw->cleanAndPack();
	_tuvw->findAdjacentTriangles();
	_sg->setNewTopology();
	_sg->updatePositionsAndNormals();
}


void Legacy_Cutter::recurseTriangleRemoval(int triangle)
{
    trianglesUVW *_tuvw = _sg->getTrianglesUVW();
    
	int *tp = _tuvw->triangleVertices(triangle);
	if(*tp < 0)
		return;
	unsigned long *adjs = _tuvw->triAdjs(triangle);
	*tp = -1;
	for(int i=0; i<3; ++i)
		recurseTriangleRemoval(adjs[i]>>2);
}
