#include "CLEjs.h"
#include <iostream>
#include <vector>
#include <array>


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
        std::cout.flush();
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
            Nan::ThrowRangeError("Buffer must be of fixed length.");
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
  CLEjs::Init(exports);
}

NODE_MODULE(CLEjs, InitAll)

Nan::Persistent<v8::Function> CLEjs::constructor;

CLEjs::CLEjs()  {
    std::cout << "CLEjs::Constructor" << std::endl;
}

CLEjs::~CLEjs() {
    std::cout << "CLEjs::Destructor" << std::endl;
} 

void CLEjs::Init(v8::Local<v8::Object> exports) {
  Nan::HandleScope scope;

  // Prepare constructor template
  v8::Local<v8::FunctionTemplate> tpl = Nan::New<v8::FunctionTemplate>(New);
  tpl->SetClassName(Nan::New("CLEjs").ToLocalChecked());
  tpl->InstanceTemplate()->SetInternalFieldCount(1);

  // Prototype
  Nan::SetPrototypeMethod(tpl,  "Set_Youngs_Modulus",  __Set_Youngs_Modulus);
  Nan::SetPrototypeMethod(tpl,  "Set_Poissons_Ratio",  __Set_Poissons_Ratio);
  Nan::SetPrototypeMethod(tpl,  "Create_Model",  __Create_Model);
  Nan::SetPrototypeMethod(tpl,  "Add_Static_Model",  __Add_Static_Model);
  Nan::SetPrototypeMethod(tpl,  "Add_Muscle_Layer",  __Add_Muscle_Layer);
  Nan::SetPrototypeMethod(tpl,  "Set_Collision_Model",  __Set_Collision_Model);
  Nan::SetPrototypeMethod(tpl,  "Set_Collision_Proxy",  __Set_Collision_Proxy);
  Nan::SetPrototypeMethod(tpl,  "Destroy_Model",  __Destroy_Model);
  Nan::SetPrototypeMethod(tpl,  "SetTextureParameters",  __SetTextureParameters);
  Nan::SetPrototypeMethod(tpl,  "GetTextureStress",  __GetTextureStress);
  Nan::SetPrototypeMethod(tpl,  "Set_Fixed_Geometry",  __Set_Fixed_Geometry);
  Nan::SetPrototypeMethod(tpl,  "Set_Fixed_Points",  __Set_Fixed_Points);
  Nan::SetPrototypeMethod(tpl,  "Set_Fixed_Segments",  __Set_Fixed_Segments);
  Nan::SetPrototypeMethod(tpl,  "Set_Fixed_Triangles",  __Set_Fixed_Triangles);
  Nan::SetPrototypeMethod(tpl,  "Set_Fixed_Volume",  __Set_Fixed_Volume);
  Nan::SetPrototypeMethod(tpl,  "Finalize_Initialization",  __Finalize_Initialization);
  Nan::SetPrototypeMethod(tpl,  "Set_Damping",  __Set_Damping);
  Nan::SetPrototypeMethod(tpl,  "Set_Density",  __Set_Density);
  Nan::SetPrototypeMethod(tpl,  "Set_Time_Step",  __Set_Time_Step);
  Nan::SetPrototypeMethod(tpl,  "Set_Hook_Stiffness",  __Set_Hook_Stiffness);
  Nan::SetPrototypeMethod(tpl,  "Set_Suture_Stiffness",  __Set_Suture_Stiffness);
  Nan::SetPrototypeMethod(tpl,  "Add_Hook",  __Add_Hook);
  Nan::SetPrototypeMethod(tpl,  "Move_Hook",  __Move_Hook);
  Nan::SetPrototypeMethod(tpl,  "Delete_Hook",  __Delete_Hook);
  Nan::SetPrototypeMethod(tpl,  "Get_Hook_Position",  __Get_Hook_Position);
  Nan::SetPrototypeMethod(tpl,  "Get_Hook_Triangles",  __Get_Hook_Triangles);
  Nan::SetPrototypeMethod(tpl,  "Get_Active_Hooks",  __Get_Active_Hooks);
  Nan::SetPrototypeMethod(tpl,  "Add_Suture",  __Add_Suture);
  Nan::SetPrototypeMethod(tpl,  "Delete_Suture",  __Delete_Suture);
  Nan::SetPrototypeMethod(tpl,  "Get_Suture_UV",  __Get_Suture_UV);
  Nan::SetPrototypeMethod(tpl,  "Get_Suture_Triangles",  __Get_Suture_Triangles);
  Nan::SetPrototypeMethod(tpl,  "Get_Active_Sutures",  __Get_Active_Sutures);
  Nan::SetPrototypeMethod(tpl,  "Advance_One_Time_Step",  __Advance_One_Time_Step);
  Nan::SetPrototypeMethod(tpl,  "WaitForSolve",  __WaitForSolve);
  Nan::SetPrototypeMethod(tpl,  "Update_Fine_Displacement",  __Update_Fine_Displacement);
  Nan::SetPrototypeMethod(tpl,  "Get_Vertices",  __Get_Vertices);
  Nan::SetPrototypeMethod(tpl,  "Get_Strain",  __Get_Strain);
  Nan::SetPrototypeMethod(tpl,  "Get_Stress",  __Get_Stress);
  Nan::SetPrototypeMethod(tpl,  "Get_Vertex_Data",  __Get_Vertex_Data);
  Nan::SetPrototypeMethod(tpl,  "Update_Embedded_Surfaces",  __Update_Embedded_Surfaces);
  Nan::SetPrototypeMethod(tpl,  "Generate_UFine_MeshMap",  __Generate_UFine_MeshMap);
  Nan::SetPrototypeMethod(tpl,  "Update_Collisions",  __Update_Collisions);
  Nan::SetPrototypeMethod(tpl,  "WriteDebug",  __WriteDebug);
  Nan::SetPrototypeMethod(tpl,  "Apply_Perturbation",  __Apply_Perturbation);

  Nan::SetAccessor(tpl->InstanceTemplate(), Nan::New("V_DATA_POSITION").ToLocalChecked(), __Getter_V_DATA_POSITION);
  Nan::SetAccessor(tpl->InstanceTemplate(), Nan::New("V_DATA_STRESS").ToLocalChecked(), __Getter_V_DATA_STRESS);
  Nan::SetAccessor(tpl->InstanceTemplate(), Nan::New("V_DATA_STRAIN").ToLocalChecked(), __Getter_V_DATA_STRAIN);

  constructor.Reset(tpl->GetFunction());
  exports->Set(Nan::New("CLEjs").ToLocalChecked(), tpl->GetFunction());
}

void CLEjs::New(const Nan::FunctionCallbackInfo<v8::Value>& info) {
  if (info.IsConstructCall()) {
    // Invoked as constructor: `new CLEjs(...)`
    //double value = info[0]->IsUndefined() ? 0 : info[0]->NumberValue();
    CLEjs* obj = new CLEjs();
    obj->Wrap(info.This());
    info.GetReturnValue().Set(info.This());
  } else {
    // Invoked as plain function `CLEjs(...)`, turn into construct call.
    const int argc = 1;
    v8::Local<v8::Value> argv[argc] = { info[0] };
    v8::Local<v8::Function> cons = Nan::New<v8::Function>(constructor);
    info.GetReturnValue().Set(cons->NewInstance(argc, argv));
  }
}


void CLEjs::__Set_Youngs_Modulus(const Nan::FunctionCallbackInfo<v8::Value>& info) {
    CLEjs* obj = ObjectWrap::Unwrap<CLEjs>(info.Holder());

    double youngs_modulus_in = info[0]->IsUndefined() ? 0.0 : info[0]->NumberValue();
    obj->Set_Youngs_Modulus( youngs_modulus_in );
};

void CLEjs::__Set_Poissons_Ratio(const Nan::FunctionCallbackInfo<v8::Value>& info) {
    CLEjs* obj = ObjectWrap::Unwrap<CLEjs>(info.Holder());

    double poissons_ratio_in = info[0]->IsUndefined() ? 0.0 : info[0]->NumberValue();
    obj->Set_Poissons_Ratio( poissons_ratio_in );
};

void CLEjs::__Create_Model(const Nan::FunctionCallbackInfo<v8::Value>& info) {
    CLEjs* obj = ObjectWrap::Unwrap<CLEjs>(info.Holder());

    if(! info.Length() > 0 )
        Nan::ThrowError("Must pass vertices and triangles");
    
    typedef float T;

    std::vector<T> vertex_in = load_std_vector<T>(info,0);
    std::vector<int> triangle_in = load_std_vector<int>(info,1);
    if( vertex_in.size() % 3 )
        Nan::ThrowRangeError("Vertex input must be an array of triplets.");
    if( triangle_in.size() % 3 )
        Nan::ThrowRangeError("Triangle input must be an array of triplets.");    
    
    double dx_in = info[2]->IsUndefined() ? 1.0 : info[2]->NumberValue();
    double refine_in = info[3]->IsUndefined() ? 1.0 : info[3]->NumberValue();

    std::cout << "Create Model: " << vertex_in.size() / 3 << " Vertices, " << triangle_in.size() / 3 << " Triangles." << std::endl;

    obj->Create_Model( vertex_in, triangle_in, dx_in, refine_in );
};

void CLEjs::__Add_Static_Model(const Nan::FunctionCallbackInfo<v8::Value>& info) {
    CLEjs* obj = ObjectWrap::Unwrap<CLEjs>(info.Holder());

    if(! info.Length() > 0 )
        Nan::ThrowError("Must pass vertices and triangles");

    typedef float T;

    std::vector<T> vertex_in = load_std_vector<T>(info,0);
    std::vector<int> triangle_in = load_std_vector<int>(info,1);
    if( vertex_in.size() % 3 )
        Nan::ThrowRangeError("Vertex input must be an array of triplets.");
    if( triangle_in.size() % 3 )
        Nan::ThrowRangeError("Triangle input must be an array of triplets.");    

    obj->Add_Static_Model( vertex_in, triangle_in );        
};

void CLEjs::__Add_Muscle_Layer(const Nan::FunctionCallbackInfo<v8::Value>& info) {
    CLEjs* obj = ObjectWrap::Unwrap<CLEjs>(info.Holder());

    if(! info.Length() > 0 )
        Nan::ThrowError("Must pass vertices and triangles");

    typedef float T;

    std::vector<T> vertex_in = load_std_vector<T>(info,0);
    std::vector<int> triangle_in = load_std_vector<int>(info,1);
    std::vector<T> fiber_dir_in = load_std_vector<T>(info,1);
    if( vertex_in.size() % 3 )
        Nan::ThrowRangeError("Vertex input must be an array of triplets.");
    if( triangle_in.size() % 3 )
        Nan::ThrowRangeError("Triangle input must be an array of triplets."); 
    if( fiber_dir_in.size() != 3 )
        Nan::ThrowRangeError("Fiber_Dir input must be a triplet.");   
    double maxstress_in = info[3]->IsUndefined() ? 1.0 : info[3]->NumberValue();
    
    obj->Add_Muscle_Layer( vertex_in, triangle_in, fiber_dir_in, maxstress_in ); 
    
};

void CLEjs::__Set_Collision_Model(const Nan::FunctionCallbackInfo<v8::Value>& info) {
    CLEjs* obj = ObjectWrap::Unwrap<CLEjs>(info.Holder());

    if(! info.Length() > 0 )
        Nan::ThrowError("Must pass vertices and triangles");

    typedef float T;

    std::vector<T> vertex_in = load_std_vector<T>(info,0);
    std::vector<int> triangle_in = load_std_vector<int>(info,1);
    if( vertex_in.size() % 3 )
        Nan::ThrowRangeError("Vertex input must be an array of triplets.");
    if( triangle_in.size() % 3 )
        Nan::ThrowRangeError("Triangle input must be an array of triplets.");    

    obj->Set_Collision_Model( vertex_in, triangle_in );     
};

void CLEjs::__Set_Collision_Proxy(const Nan::FunctionCallbackInfo<v8::Value>& info) {
    CLEjs* obj = ObjectWrap::Unwrap<CLEjs>(info.Holder());

    if(! info.Length() > 0 )
        Nan::ThrowError("Must pass vertices and triangles");

    typedef float T;

    std::vector<T> vertex_in = load_std_vector<T>(info,0);
    std::vector<int> triangle_in = load_std_vector<int>(info,1);
    if( vertex_in.size() % 3 )
        Nan::ThrowRangeError("Vertex input must be an array of triplets.");
    if( triangle_in.size() % 3 )
        Nan::ThrowRangeError("Triangle input must be an array of triplets.");    

    obj->Set_Collision_Proxy( vertex_in, triangle_in ); 
};

void CLEjs::__Destroy_Model(const Nan::FunctionCallbackInfo<v8::Value>& info) {
    CLEjs* obj = ObjectWrap::Unwrap<CLEjs>(info.Holder());
    obj->Destroy_Model();
};

void CLEjs::__SetTextureParameters(const Nan::FunctionCallbackInfo<v8::Value>& info) {
    CLEjs* obj = ObjectWrap::Unwrap<CLEjs>(info.Holder());
    typedef float T;
    std::vector<int> dimensions_in = load_std_vector<int>(info,0);
    std::vector<T> texcoords_in = load_std_vector<T>(info,1);
    std::vector<int> textriangles_in = load_std_vector<int>(info,2);
    std::vector<int> coord_map_in = load_std_vector<int>(info,3);
    obj->SetTextureParameters( dimensions_in, texcoords_in, textriangles_in, coord_map_in);    
};

void CLEjs::__GetTextureStress(const Nan::FunctionCallbackInfo<v8::Value>& info) {
    CLEjs* obj = ObjectWrap::Unwrap<CLEjs>(info.Holder());
    typedef float T;
    std::vector<T> texdata;
    obj->GetTextureStress( texdata );       
    info.GetReturnValue().Set(Nan::CopyBuffer(reinterpret_cast<const char*>(texdata.data()),
                                              texdata.size()*sizeof(T)).ToLocalChecked());    
};

void CLEjs::__Set_Fixed_Geometry(const Nan::FunctionCallbackInfo<v8::Value>& info) {
    CLEjs* obj = ObjectWrap::Unwrap<CLEjs>(info.Holder());
    typedef float T;
    std::vector<T> vertices_in = load_std_vector<T>(info,0);
    std::vector<int> points_in = load_std_vector<int>(info,1);
    std::vector<int> segments_in = load_std_vector<int>(info,2);
    std::vector<int> triangles_in = load_std_vector<int>(info,3);
    obj->Set_Fixed_Geometry( vertices_in, points_in, segments_in, triangles_in);
};

void CLEjs::__Set_Fixed_Points(const Nan::FunctionCallbackInfo<v8::Value>& info) {
    CLEjs* obj = ObjectWrap::Unwrap<CLEjs>(info.Holder());
    typedef float T;
    std::vector<T> vertices_in = load_std_vector<T>(info,0);
    obj->Set_Fixed_Points( vertices_in);
};

void CLEjs::__Set_Fixed_Segments(const Nan::FunctionCallbackInfo<v8::Value>& info) {
    CLEjs* obj = ObjectWrap::Unwrap<CLEjs>(info.Holder());
    typedef float T;
    std::vector<T> vertices_in = load_std_vector<T>(info,0);
    std::vector<int> segments_in = load_std_vector<int>(info,1);
    obj->Set_Fixed_Segments( vertices_in, segments_in);
};

void CLEjs::__Set_Fixed_Triangles(const Nan::FunctionCallbackInfo<v8::Value>& info) {
    CLEjs* obj = ObjectWrap::Unwrap<CLEjs>(info.Holder());
    typedef float T;
    std::vector<T> vertices_in = load_std_vector<T>(info,0);
    std::vector<int> triangles_in = load_std_vector<int>(info,1);
    obj->Set_Fixed_Triangles( vertices_in, triangles_in);
};

void CLEjs::__Set_Fixed_Volume(const Nan::FunctionCallbackInfo<v8::Value>& info) {
    CLEjs* obj = ObjectWrap::Unwrap<CLEjs>(info.Holder());
    typedef float T;
    std::vector<T> min_corner_in = load_std_vector<T>(info,0);
    std::vector<T> max_corner_in = load_std_vector<T>(info,1);
    obj->Set_Fixed_Volume( min_corner_in, max_corner_in );
};

void CLEjs::__Finalize_Initialization(const Nan::FunctionCallbackInfo<v8::Value>& info) {
    CLEjs* obj = ObjectWrap::Unwrap<CLEjs>(info.Holder());
    obj->Finalize_Initialization();
};

void CLEjs::__Set_Damping(const Nan::FunctionCallbackInfo<v8::Value>& info) {
    CLEjs* obj = ObjectWrap::Unwrap<CLEjs>(info.Holder());

    double rayleigh_coefficient_in = info[0]->IsUndefined() ? 0.0 : info[0]->NumberValue();
    obj->Set_Damping( rayleigh_coefficient_in );
};

void CLEjs::__Set_Density(const Nan::FunctionCallbackInfo<v8::Value>& info) {
    CLEjs* obj = ObjectWrap::Unwrap<CLEjs>(info.Holder());

    double density_in = info[0]->IsUndefined() ? 1.0 : info[0]->NumberValue();
    obj->Set_Density( density_in );
};

void CLEjs::__Set_Time_Step(const Nan::FunctionCallbackInfo<v8::Value>& info) {
    CLEjs* obj = ObjectWrap::Unwrap<CLEjs>(info.Holder());

    double dt_in = info[0]->IsUndefined() ? 1.0 : info[0]->NumberValue();
    obj->Set_Time_Step( dt_in );
};

void CLEjs::__Set_Hook_Stiffness(const Nan::FunctionCallbackInfo<v8::Value>& info) {
    CLEjs* obj = ObjectWrap::Unwrap<CLEjs>(info.Holder());
    double hook_stiffness_in = info[0]->IsUndefined() ? 1.0 : info[0]->NumberValue();
    obj->Set_Hook_Stiffness( hook_stiffness_in );
};

void CLEjs::__Set_Suture_Stiffness(const Nan::FunctionCallbackInfo<v8::Value>& info) {
    CLEjs* obj = ObjectWrap::Unwrap<CLEjs>(info.Holder());
    double suture_stiffness_in = info[0]->IsUndefined() ? 1.0 : info[0]->NumberValue();
    obj->Set_Suture_Stiffness( suture_stiffness_in );
};

void CLEjs::__Add_Hook(const Nan::FunctionCallbackInfo<v8::Value>& info) {
    CLEjs* obj = ObjectWrap::Unwrap<CLEjs>(info.Holder());
    typedef float T;
    int ret = -1;
    
    if( info.Length() == 1 ){
        std::array<T,3> location = load_std_array<T,3>(info,0);      
        ret = obj->Add_Single_Point_Constraint( location );
    }
    if( info.Length() == 2 ){
        int triangle_id_in = info[0]->NumberValue();
        std::array<T,3> location = load_std_array<T,3>(info,1);      
        ret = obj->Add_Single_Point_Constraint( triangle_id_in, location );
    }
    
    info.GetReturnValue().Set(ret);    
};

void CLEjs::__Move_Hook(const Nan::FunctionCallbackInfo<v8::Value>& info) {
    CLEjs* obj = ObjectWrap::Unwrap<CLEjs>(info.Holder());
    typedef float T;
    int hook_id_in = info[0]->NumberValue();
    std::array<T,3> location_in = load_std_array<T,3>(info,1);      
    obj->Move_Single_Point_Constraint( hook_id_in, location_in );
};

void CLEjs::__Delete_Hook(const Nan::FunctionCallbackInfo<v8::Value>& info) {
    CLEjs* obj = ObjectWrap::Unwrap<CLEjs>(info.Holder());
    int hook_id_in = info[0]->NumberValue();
    obj->Delete_Single_Point_Constraint( hook_id_in );
};

void CLEjs::__Get_Hook_Position(const Nan::FunctionCallbackInfo<v8::Value>& info) {
    CLEjs* obj = ObjectWrap::Unwrap<CLEjs>(info.Holder());
    std::vector<int> hook_ids = load_std_vector<int>(info,0);
    std::vector<std::array<float,3> > locations;
    obj->Get_Single_Point_Constraint_Position( hook_ids, locations );
    info.GetReturnValue().Set(Nan::CopyBuffer(reinterpret_cast<const char*>(locations.data()),
                                              3*locations.size()*sizeof(float)).ToLocalChecked());
};

void CLEjs::__Get_Hook_Triangles(const Nan::FunctionCallbackInfo<v8::Value>& info) {
    CLEjs* obj = ObjectWrap::Unwrap<CLEjs>(info.Holder());
    std::vector<int> hook_ids = load_std_vector<int>(info,0);
    std::vector< int > triangles;
    obj->Get_Single_Point_Constraint_Triangles( hook_ids, triangles );
    info.GetReturnValue().Set(Nan::CopyBuffer(reinterpret_cast<const char*>(triangles.data()),
                                              triangles.size()*sizeof(int)).ToLocalChecked());
};

void CLEjs::__Get_Active_Hooks(const Nan::FunctionCallbackInfo<v8::Value>& info) {
    CLEjs* obj = ObjectWrap::Unwrap<CLEjs>(info.Holder());
    std::vector<int> hook_ids;
    obj->Get_Active_Single_Point_Constraints(hook_ids);
    info.GetReturnValue().Set(Nan::CopyBuffer(reinterpret_cast<const char*>(hook_ids.data()),
                                              hook_ids.size()*sizeof(int)).ToLocalChecked());
}


void CLEjs::__Add_Suture(const Nan::FunctionCallbackInfo<v8::Value>& info) {
    CLEjs* obj = ObjectWrap::Unwrap<CLEjs>(info.Holder());
    typedef float T;
    int ret = -1;
    
    if( info.Length() == 4 ){
        int triangle1_id_in = info[0]->NumberValue();
        std::array<T,3> weights1 = load_std_array<T,3>(info,1);      
        int triangle2_id_in = info[2]->NumberValue();
        std::array<T,3> weights2 = load_std_array<T,3>(info,3);      
        ret = obj->Add_Dual_Point_Constraint( triangle1_id_in, weights1, triangle2_id_in, weights2 );       
    }
    if( info.Length() == 2 ){
        std::array<T,3> location1 = load_std_array<T,3>(info,0);
        std::array<T,3> location2 = load_std_array<T,3>(info,1);      
        ret = obj->Add_Dual_Point_Constraint( location1, location2 );
    }
    
    info.GetReturnValue().Set(ret);
};

void CLEjs::__Delete_Suture(const Nan::FunctionCallbackInfo<v8::Value>& info) {
    CLEjs* obj = ObjectWrap::Unwrap<CLEjs>(info.Holder());
    int suture_id_in = info[0]->NumberValue();
    obj->Delete_Dual_Point_Constraint( suture_id_in );
};

void CLEjs::__Get_Suture_UV(const Nan::FunctionCallbackInfo<v8::Value>& info) {
    CLEjs* obj = ObjectWrap::Unwrap<CLEjs>(info.Holder());
    std::vector<int> suture_ids = load_std_vector<int>(info,0);
    std::vector<std::array<float,6> > uvs;
    obj->Get_Dual_Point_Constraint_UV( suture_ids, uvs );
    info.GetReturnValue().Set(Nan::CopyBuffer(reinterpret_cast<const char*>(uvs.data()),
                                              6*uvs.size()*sizeof(float)).ToLocalChecked());
};

void CLEjs::__Get_Suture_Triangles(const Nan::FunctionCallbackInfo<v8::Value>& info) {
    CLEjs* obj = ObjectWrap::Unwrap<CLEjs>(info.Holder());
    std::vector<int> suture_ids = load_std_vector<int>(info,0);
    std::vector< std::array<int,2> > triangles;
    obj->Get_Dual_Point_Constraint_Triangles( suture_ids, triangles );
    info.GetReturnValue().Set(Nan::CopyBuffer(reinterpret_cast<const char*>(triangles.data()),
                                              triangles.size()*2*sizeof(int)).ToLocalChecked());
};

void CLEjs::__Get_Active_Sutures(const Nan::FunctionCallbackInfo<v8::Value>& info) {
    CLEjs* obj = ObjectWrap::Unwrap<CLEjs>(info.Holder());
    std::vector<int> suture_ids;
    obj->Get_Active_Dual_Point_Constraints(suture_ids);
    info.GetReturnValue().Set(Nan::CopyBuffer(reinterpret_cast<const char*>(suture_ids.data()),
                                              suture_ids.size()*sizeof(int)).ToLocalChecked());
}

void CLEjs::__Advance_One_Time_Step(const Nan::FunctionCallbackInfo<v8::Value>& info) {
    CLEjs* obj = ObjectWrap::Unwrap<CLEjs>(info.Holder());
    obj->Advance_One_Time_Step();
};

void CLEjs::__WaitForSolve(const Nan::FunctionCallbackInfo<v8::Value>& info) {
    CLEjs* obj = ObjectWrap::Unwrap<CLEjs>(info.Holder());
    obj->WaitForSolve();
};

void CLEjs::__Update_Fine_Displacement(const Nan::FunctionCallbackInfo<v8::Value>& info) {
    CLEjs* obj = ObjectWrap::Unwrap<CLEjs>(info.Holder());
    obj->Update_Fine_Displacement();
};

void CLEjs::__Get_Vertices(const Nan::FunctionCallbackInfo<v8::Value>& info) {
    CLEjs* obj = ObjectWrap::Unwrap<CLEjs>(info.Holder());
    typedef float T;

    std::vector<T> vertex_data;
    
    if( info.Length() == 1 ){
        int since_frame = info[0]->NumberValue();
        obj->Get_Vertices( vertex_data, since_frame );
    }
    else{
        obj->Get_Vertices( vertex_data );
    }

    info.GetReturnValue().Set(Nan::CopyBuffer(reinterpret_cast<const char*>(vertex_data.data()),
                                              vertex_data.size()*sizeof(T)).ToLocalChecked());               
};

void CLEjs::__Get_Strain(const Nan::FunctionCallbackInfo<v8::Value>& info) {
    CLEjs* obj = ObjectWrap::Unwrap<CLEjs>(info.Holder());
    typedef float T;
    std::vector<T> strain_data;
    obj->Get_Strain( strain_data );
    info.GetReturnValue().Set(Nan::CopyBuffer(reinterpret_cast<const char*>(strain_data.data()),
                                              strain_data.size()*sizeof(T)).ToLocalChecked());      
};

void CLEjs::__Get_Stress(const Nan::FunctionCallbackInfo<v8::Value>& info) {
    CLEjs* obj = ObjectWrap::Unwrap<CLEjs>(info.Holder());
    typedef float T;
    std::vector<T> stress_data;
    obj->Get_Stress( stress_data );
    info.GetReturnValue().Set(Nan::CopyBuffer(reinterpret_cast<const char*>(stress_data.data()),
                                              stress_data.size()*sizeof(T)).ToLocalChecked());      
};

void CLEjs::__Get_Vertex_Data(const Nan::FunctionCallbackInfo<v8::Value>& info) {
    Nan:: HandleScope scope;
    v8::Local<v8::Object> results = Nan::New<v8::Object>();
    
    CLEjs* obj = ObjectWrap::Unwrap<CLEjs>(info.Holder());
    typedef float T;
    std::vector<T> v_data;
    int type_in = info[0]->NumberValue();
    int since_frame_in = info[1]->NumberValue();
    obj->Get_Vertex_Data( v_data, type_in, since_frame_in );

    results->Set( Nan::New("frame").ToLocalChecked(),
                  Nan::New(since_frame_in) );

    results->Set( Nan::New("vdata").ToLocalChecked(),
                  Nan::CopyBuffer(reinterpret_cast<const char*>(v_data.data()),
                                  v_data.size()*sizeof(T)).ToLocalChecked() );

    
    info.GetReturnValue().Set(results);
};

void CLEjs::
__Getter_V_DATA_POSITION(v8::Local<v8::String> property,
                         const Nan::PropertyCallbackInfo<v8::Value>& info){
    info.GetReturnValue().Set( V_DATA_TYPE::POSITION );
}

void CLEjs::
__Getter_V_DATA_STRESS(v8::Local<v8::String> property,
                       const Nan::PropertyCallbackInfo<v8::Value>& info){
    info.GetReturnValue().Set( V_DATA_TYPE::STRESS );
}

void CLEjs::
__Getter_V_DATA_STRAIN(v8::Local<v8::String> property,
                       const Nan::PropertyCallbackInfo<v8::Value>& info){
    info.GetReturnValue().Set( V_DATA_TYPE::STRAIN );
}

void CLEjs::__Update_Embedded_Surfaces(const Nan::FunctionCallbackInfo<v8::Value>& info) {
    CLEjs* obj = ObjectWrap::Unwrap<CLEjs>(info.Holder());
    obj->Update_Embedded_Surfaces();
};

void CLEjs::__Generate_UFine_MeshMap(const Nan::FunctionCallbackInfo<v8::Value>& info) {
    CLEjs* obj = ObjectWrap::Unwrap<CLEjs>(info.Holder());
    obj->Generate_UFine_MeshMap();
};

void CLEjs::__Update_Collisions(const Nan::FunctionCallbackInfo<v8::Value>& info) {
    CLEjs* obj = ObjectWrap::Unwrap<CLEjs>(info.Holder());
    obj->Update_Collisions();
};

void CLEjs::__WriteDebug(const Nan::FunctionCallbackInfo<v8::Value>& info) {
    CLEjs* obj = ObjectWrap::Unwrap<CLEjs>(info.Holder());
    bool andDie_in = info[0]->BooleanValue();
    obj->WriteDebug( andDie_in );
};

void CLEjs::__Apply_Perturbation(const Nan::FunctionCallbackInfo<v8::Value>& info) {
    CLEjs* obj = ObjectWrap::Unwrap<CLEjs>(info.Holder());
    float perturb_in = info[0]->NumberValue();
    obj->Apply_Perturbation( perturb_in );
};


