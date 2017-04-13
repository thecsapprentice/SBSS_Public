#ifndef __LEGACY_CUTTER_H__
#define __LEGACY_CUTTER_H__

#include <nan.h>
#include <vector>
#include <memory>

class skinGraphics;
class incision;

//using namespace PhysBAM;
void InitAll(v8::Local<v8::Object> exports);

class Legacy_Cutter : public Nan::ObjectWrap {
 public:
  static void Init(v8::Local<v8::Object> exports);
 
 private:

 explicit Legacy_Cutter();
  ~Legacy_Cutter();
  
  // Helper methods
  bool texturePickCode(const int triangle, const float (&uv)[2], float (&uvw)[3], float &triangleDuv);
  void CreateIncision( const std::vector<int>& path_triangles, const std::vector<float>& path_uvs, 
                       const std::vector<float>& path_positions, const std::vector<float>& path_normals, 
                       const bool edge_start, const bool edge_end );
  void CreateExcision( int triangle );
  void recurseTriangleRemoval(int triangle);
  
  // Wrapper methods
  static void New(const Nan::FunctionCallbackInfo<v8::Value>& info);
  static Nan::Persistent<v8::Function> constructor;
  
  static void __ParseFile( const Nan::FunctionCallbackInfo<v8::Value>& info );
  static void __ParseFileStatic( const Nan::FunctionCallbackInfo<v8::Value>& info );
  static void __makeIncision( const Nan::FunctionCallbackInfo<v8::Value>& info );
  static void __makeExcision( const Nan::FunctionCallbackInfo<v8::Value>& info );
  
  static void __GetJavascriptVertex( const Nan::FunctionCallbackInfo<v8::Value>& info );
  static void __GetJavascriptTopology( const Nan::FunctionCallbackInfo<v8::Value>& info );
  static void __GetJavascriptUV( const Nan::FunctionCallbackInfo<v8::Value>& info );

  static void __GetRawVertex( const Nan::FunctionCallbackInfo<v8::Value>& info );
  static void __GetRawTopology( const Nan::FunctionCallbackInfo<v8::Value>& info );

  static void __RemapVertexData( const Nan::FunctionCallbackInfo<v8::Value>& info );


  
  // Members
  std::unique_ptr<incision> _incis;
  std::unique_ptr<skinGraphics> _sg; 
     
};

#endif

